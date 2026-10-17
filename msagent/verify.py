"""Criterion-by-criterion verification of one study against the requested experiment.

Every criterion ends up ``confirmed``, ``excluded`` or ``unknown`` with the evidence that
decided it. Unknown is never silently treated as a pass. Counting is done at the right
grain: individuals vs biological samples vs runs (data files) vs fractions/channels.
"""

from __future__ import annotations

import re
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass

from . import vocab
from .criteria import CELL_TYPES, Criteria
from .store import ProjectData

CONFIRMED, EXCLUDED, UNKNOWN = "confirmed", "excluded", "unknown"

TECHNICAL_RE = re.compile(
    r"\b(reference|quality control|qc|pool(ed)?|blank|empty|buffer|standard|carrier|bulk control|plate control|"
    r"negative control|spike[- ]?in|library)\b"
)
ROLE_ATTRS = vocab.DISEASE_ATTRS + ("sample type", "tissue type", "tumor status", "tumour status", "specimen type",
                                    "sample group", "group", "histology", "pathology", "tumor", "tumour")
GROUP_NAMES = {
    "tumour_vs_normal": ("tumour", "normal"),
    "disease_vs_control": ("disease", "control"),
    "treated_vs_control": ("treated", "control"),
    "mutant_vs_wildtype": ("mutant", "wild type"),
}


@dataclass
class Check:
    id: str
    label: str
    status: str
    detail: str
    kind: str = "required"
    source: str = ""
    hint: str = ""

    def to_dict(self):
        return asdict(self)


@dataclass
class Sample:
    name: str
    attrs: dict[str, list[str]]
    runs: list[dict]
    technical: str | None = None
    cell_line: str | None = None
    group: str | None = None
    group_evidence: str | None = None
    individual: str | None = None

    def vals(self, attrs) -> list[tuple[str, str]]:
        out = []
        for a in attrs:
            for v in self.attrs.get(a, []):
                if v and not vocab.is_missing(v):
                    out.append((a, v))
        return out


def _q(v: str) -> str:
    return f"'{v}'" if len(v) < 60 else f"'{v[:57]}...'"


def _top(counter: Counter, n: int = 3) -> str:
    return ", ".join(f"{_q(k)} ({v})" for k, v in counter.most_common(n))


def _sample_views(pd: ProjectData) -> list[Sample]:
    runs = pd.runs_by_sample()
    out = []
    for name, attrs in pd.samples.items():
        s = Sample(name=name, attrs=attrs, runs=runs.get(name, []))
        for a, v in s.vals(("sample type",)):
            if TECHNICAL_RE.search(vocab.norm(v)):
                s.technical = f"sample type = {_q(v)}"
        for a, v in s.vals(vocab.INDIVIDUAL_ATTRS):
            nv = vocab.norm(v)
            if nv in {"pooled", "pool", "mixed"}:
                s.technical = s.technical or f"{a} = {_q(v)}"
            elif s.individual is None:
                s.individual = v
        for a, v in s.vals(("material type",)):
            if vocab.norm(v) in {"synthetic", "recombinant protein"}:
                s.technical = s.technical or f"material type = {_q(v)}"
        s.cell_line = _cell_line_reason(s)
        out.append(s)
    return out


def _cell_line_reason(s: Sample) -> str | None:
    for a, v in s.vals(("cell line", "cellosaurus name", "cellosaurus accession")):
        if vocab.norm(v) not in {"no", "none", "primary", "primary cells"}:
            return f"{a} = {_q(v)}"
    for a, v in s.vals(("material type",)):
        if "cell line" in vocab.norm(v) or "cell culture" in vocab.norm(v):
            return f"{a} = {_q(v)}"
    for a, v in s.vals(vocab.TISSUE_ATTRS):
        if vocab.norm(v) in {"cell culture", "cell suspension culture", "cultured cell", "cell lysate"}:
            return f"{a} = {_q(v)}"
    for a, v in s.vals(("cell type",)):
        if vocab.looks_like_cell_line(v):
            return f"{a} = {_q(v)}"
    # HeLa is very often annotated only through its disease (Cellosaurus: cervical adenocarcinoma)
    dz = [v for _, v in s.vals(("disease",))]
    if any("cervical adenocarcinoma" in vocab.norm(v) for v in dz) and not s.individual:
        tissue = [vocab.norm(v) for _, v in s.vals(vocab.TISSUE_ATTRS)]
        if not tissue or all("cervix" in t or "cervical" in t for t in tissue):
            return "likely HeLa: disease = 'cervical adenocarcinoma' with no patient annotation"
    return None


def _pattern(table: dict, key: str) -> vocab.TermPattern:
    return table.get(key) or vocab.TermPattern(key, [vocab.norm(key)])


def _attr_check(samples: list[Sample], attrs, patterns: list[vocab.TermPattern], cid: str, label: str, what: str):
    """Generic sample-level attribute check. Returns (Check, set of matching sample names)."""
    matched, seen = set(), Counter()
    annotated = 0
    for s in samples:
        vals = s.vals(attrs)
        if vals:
            annotated += 1
        for a, v in vals:
            seen[v] += 1
            if any(p.matches(v) for p in patterns):
                matched.add(s.name)
    src = "SDRF " + "/".join(f"[{a}]" for a in attrs[:2])
    if matched:
        hit = Counter(v for s in samples if s.name in matched for _, v in s.vals(attrs) if any(p.matches(v) for p in patterns))
        return Check(cid, label, CONFIRMED, f"{len(matched)}/{len(samples)} samples: {what} = {_top(hit)}", source=src), matched
    if annotated:
        return Check(cid, label, EXCLUDED, f"{what} annotated for {annotated} samples, none match; found {_top(seen)}", source=src), matched
    return Check(cid, label, UNKNOWN, f"{what} not annotated in SDRF", source=src), matched


def _role_disease(s: Sample, factor_attrs: set[str]) -> tuple[str | None, str | None]:
    attrs = tuple(ROLE_ATTRS) + tuple(sorted(factor_attrs))
    vals = s.vals(attrs)
    for a, v in vals:
        if vocab.is_normal_value(v):
            return "normal", f"{a} = {_q(v)}"
    for a, v in s.vals(vocab.TISSUE_ATTRS):
        if vocab.is_normal_value(v):
            return "normal", f"{a} = {_q(v)}"
    for a, v in vals:
        if vocab.is_tumour_value(v):
            return "tumour", f"{a} = {_q(v)}"
    for a, v in s.vals(vocab.TISSUE_ATTRS):
        if vocab.is_tumour_value(v):
            return "tumour", f"{a} = {_q(v)}"
    sn = vocab.norm(s.name)
    if re.search(r"\b(normal|adjacent|non[- ]?tumou?r\w*|nat|para[- ]?(tumou?r|cancer)\w*|healthy)\b", sn):
        return "normal", f"sample name {_q(s.name)} (inferred)"
    if re.search(r"\b(tumou?r|cancer|carcinoma|hcc)\b", sn):
        return "tumour", f"sample name {_q(s.name)} (inferred)"
    return None, None


_TREAT_FACTOR = re.compile(r"treat|compound|drug|dose|stimul|inhibit|exposure|perturb|agent|condition")
_STRICT_WT_RE = re.compile(r"\b(wild[- ]?type|wildtype|wt|\+/\+)\b")
_GENO_FACTOR = re.compile(r"genotype|genetic|knock|mutation|strain|construct|sirna|shrna")
_CONTROL_EXTRA = {"none", "no treatment", "no", "0", "0 nm", "0 um", "0 mm", "0 μm", "0 ug/ml", "pbs", "water", "saline", "wt", "wild type"}


def _role_treatment(s: Sample, factor_attrs: set[str]):
    attrs = tuple(vocab.TREATMENT_ATTRS) + tuple(a for a in sorted(factor_attrs) if _TREAT_FACTOR.search(a))
    vals = s.vals(attrs)
    if not vals:
        return None, None
    for a, v in vals:
        nv = vocab.norm(v)
        if vocab.CONTROL_RE.search(nv) or nv in _CONTROL_EXTRA:
            return "control", f"{a} = {_q(v)}"
    a, v = vals[0]
    return "treated", f"{a} = {_q(v)}"


def _role_genotype(s: Sample, factor_attrs: set[str]):
    attrs = tuple(vocab.GENOTYPE_ATTRS) + tuple(a for a in sorted(factor_attrs) if _GENO_FACTOR.search(a))
    vals = s.vals(attrs)
    for a, v in vals:
        # "control" means wild type only in a genotype column; in phenotype/condition it means healthy/untreated
        pattern = vocab.WILDTYPE_RE if a in ("genotype", "genetic modification") or _GENO_FACTOR.search(a) else _STRICT_WT_RE
        if pattern.search(vocab.norm(v)):
            return "wild type", f"{a} = {_q(v)}"
    for a, v in vals:
        if a in ("genotype", "genetic modification") or vocab.MUTANT_RE.search(vocab.norm(v)):
            return "mutant", f"{a} = {_q(v)}"
    return None, None


_STEM_RE = re.compile(
    r"(tumou?rs?|normal|cancer|carcinoma|adjacent|non[-_ ]?tumou?r\w*|para[-_ ]?\w*|nat|healthy|control|ctrl|"
    r"(?<=[\d_\- ])[tn](?=$|[_\- ])|^[tn](?=[_\- \d])|ca$|nt$)"
)


def _name_stem(name: str) -> str:
    s = _STEM_RE.sub("", vocab.norm(name))
    return re.sub(r"[\s_\-\.]+", "_", s).strip("_")


def verify(pd: ProjectData, c: Criteria) -> dict:
    samples_all = _sample_views(pd)
    bio = [s for s in samples_all if not s.technical]
    checks: list[Check] = []
    pool = list(bio)  # candidates for the cohort; narrowed by each sample-level criterion
    warnings: list[str] = []

    def narrow(check: Check, matched: set[str]):
        nonlocal pool
        checks.append(check)
        if check.status == CONFIRMED:
            pool = [s for s in pool if s.name in matched]

    # ---- biology ------------------------------------------------------------------
    if c.organisms:
        matched = {s.name for s in pool if any(vocab.normalize_organism(v) in c.organisms for _, v in s.vals(("organism",)))}
        seen = Counter(vocab.normalize_organism(v) for s in pool for _, v in s.vals(("organism",)))
        if matched:
            ck = Check("organism", "Organism", CONFIRMED, f"{len(matched)}/{len(pool)} samples: organism = {_top(Counter({o: seen[o] for o in c.organisms if seen[o]}))}", source="SDRF characteristics[organism]")
            others = [o for o in seen if o not in c.organisms]
            if others:
                ck.hint = "also contains: " + ", ".join(others[:3])
        elif seen:
            ck = Check("organism", "Organism", EXCLUDED, f"organisms in study: {_top(seen)}", source="SDRF characteristics[organism]")
        else:
            ck = Check("organism", "Organism", UNKNOWN, "organism not annotated", source="SDRF characteristics[organism]")
        narrow(ck, matched)

    if c.tissues:
        pats = [_pattern(vocab.TISSUES, t) for t in c.tissues]
        ck, matched = _attr_check(pool, vocab.TISSUE_ATTRS, pats, "tissue", "Tissue / biofluid", "organism part")
        narrow(ck, matched)

    if c.cell_types:
        pats = [CELL_TYPES.get(t) or vocab.TermPattern(t, [t]) for t in c.cell_types]
        ck, matched = _attr_check(pool, ("cell type", "organism part"), pats, "cell_type", "Cell type", "cell type")
        narrow(ck, matched)

    if c.cell_lines:
        pats = [vocab.TermPattern(cl, [cl]) for cl in c.cell_lines]
        ck, matched = _attr_check(pool, ("cell line", "cellosaurus name", "cell type"), pats, "cell_line", "Cell line", "cell line")
        if ck.status != CONFIRMED:
            # cell line names often live only in the sample name
            byname = {s.name for s in pool if any(p.matches(s.name) for p in pats)}
            if byname:
                ck = Check("cell_line", "Cell line", CONFIRMED, f"{len(byname)} samples named after {', '.join(c.cell_lines)} (inferred from sample names)", source="SDRF source name")
                matched = byname
        narrow(ck, matched)

    disease_pats = [_pattern(vocab.DISEASES, d) for d in c.diseases]
    groups: dict[str, list[Sample]] = {}
    if c.diseases and c.comparison not in ("tumour_vs_normal", "disease_vs_control"):
        ck, matched = _attr_check(pool, vocab.DISEASE_ATTRS + tuple(sorted(pd.factor_attrs)), disease_pats, "disease", "Disease", "disease")
        narrow(ck, matched)

    # ---- exclusion of cell lines (before the comparison so groups are built from usable samples) ----
    if c.exclude_cell_lines:
        cl = [s for s in pool if s.cell_line]
        prim = [s for s in pool if not s.cell_line]
        tissue_annot = sum(1 for s in prim if s.vals(vocab.TISSUE_ATTRS + ("material type",)))
        if not pool:
            ck = Check("no_cell_lines", "Not cell lines", UNKNOWN, "no samples left to check")
        elif not cl and tissue_annot:
            ck = Check("no_cell_lines", "Not cell lines", CONFIRMED, f"no cell-line annotation; {tissue_annot} samples annotated with a tissue/material", source="SDRF characteristics[cell line]/[material type]/[organism part]")
        elif not cl:
            ck = Check("no_cell_lines", "Not cell lines", UNKNOWN, "no cell-line annotation, but tissue/material type not stated either")
        elif not prim:
            ck = Check("no_cell_lines", "Not cell lines", EXCLUDED, f"all {len(cl)} samples are cell-line material: {cl[0].cell_line}")
        else:
            ck = Check("no_cell_lines", "Not cell lines", CONFIRMED, f"{len(prim)} primary samples kept; {len(cl)} cell-line samples dropped from cohort (e.g. {cl[0].cell_line})")
        checks.append(ck)
        if ck.status != UNKNOWN:
            pool = prim
    if c.only_cell_lines:
        cl = [s for s in pool if s.cell_line]
        if cl:
            ck = Check("only_cell_lines", "Cell-line material", CONFIRMED, f"{len(cl)} cell-line samples, e.g. {cl[0].cell_line}")
            pool = cl
        elif pool:
            ck = Check("only_cell_lines", "Cell-line material", EXCLUDED if any(s.vals(vocab.TISSUE_ATTRS) for s in pool) else UNKNOWN, "no cell-line annotation found")
        else:
            ck = Check("only_cell_lines", "Cell-line material", UNKNOWN, "no samples left to check")
        checks.append(ck)

    # ---- comparison ---------------------------------------------------------------
    if c.comparison:
        a_name, b_name = GROUP_NAMES[c.comparison]
        for s in pool:
            if c.comparison == "tumour_vs_normal":
                role, ev = _role_disease(s, pd.factor_attrs)
                if role == "tumour" and c.diseases and c.diseases != ["cancer"]:
                    sample_dz = s.vals(vocab.DISEASE_ATTRS)
                    if sample_dz and not any(p.matches(v) for _, v in sample_dz for p in disease_pats):
                        # tumour of another type; keep only if the project-level disease matches
                        if not any(p.matches(v) for x in pool for _, v in x.vals(vocab.DISEASE_ATTRS) for p in disease_pats):
                            role = None
                s.group = {"tumour": a_name, "normal": b_name}.get(role or "")
            elif c.comparison == "disease_vs_control":
                role, ev = _role_disease(s, pd.factor_attrs)
                dz_match = any(p.matches(v) for _, v in s.vals(vocab.DISEASE_ATTRS + tuple(sorted(pd.factor_attrs))) for p in disease_pats) if disease_pats else role == "tumour"
                if role == "normal":
                    s.group = b_name
                elif dz_match:
                    s.group = a_name
                    ev = ev or "disease matches"
                    hits = [f"{a} = {_q(v)}" for a, v in s.vals(vocab.DISEASE_ATTRS) if any(p.matches(v) for p in disease_pats)]
                    ev = hits[0] if hits else ev
            elif c.comparison == "treated_vs_control":
                role, ev = _role_treatment(s, pd.factor_attrs)
                s.group = {"treated": a_name, "control": b_name}.get(role or "")
            else:
                role, ev = _role_genotype(s, pd.factor_attrs)
                s.group = {"mutant": a_name, "wild type": b_name}.get(role or "")
            s.group_evidence = ev if s.group else None
        ga = [s for s in pool if s.group == a_name]
        gb = [s for s in pool if s.group == b_name]
        groups = {a_name: ga, b_name: gb}
        ev_a = Counter(s.group_evidence for s in ga)
        ev_b = Counter(s.group_evidence for s in gb)
        label = {"tumour_vs_normal": "Tumour vs normal", "disease_vs_control": "Disease vs control",
                 "treated_vs_control": "Treated vs control", "mutant_vs_wildtype": "Mutant vs wild type"}[c.comparison]
        if ga and gb:
            ck = Check("comparison", label, CONFIRMED, f"{len(ga)} {a_name} ({_top(ev_a, 1)}) and {len(gb)} {b_name} ({_top(ev_b, 1)})", source="SDRF characteristics / factor values")
            if any("(inferred)" in (s.group_evidence or "") for s in ga + gb):
                ck.hint = "some group assignments inferred from sample names"
            if c.comparison == "tumour_vs_normal" and any("adenoma" in vocab.norm(e or "") and "adenocarcinoma" not in vocab.norm(e or "") for e in ev_a):
                ck.hint = "tumours annotated as adenoma (benign); check this fits the question"
        elif ga or gb:
            present = a_name if ga else b_name
            missing = b_name if ga else a_name
            ck = Check("comparison", label, EXCLUDED, f"only {present} samples ({len(ga or gb)}); no {missing} samples annotated")
        else:
            ck = Check("comparison", label, UNKNOWN, "sample annotations do not state the groups needed for this comparison")
        checks.append(ck)
        if ck.status == CONFIRMED:
            pool = ga + gb
        elif ck.status == EXCLUDED:
            pool = ga or gb

        if c.paired:
            inds_a = {s.individual for s in ga if s.individual}
            inds_b = {s.individual for s in gb if s.individual}
            both = inds_a & inds_b
            if both:
                ck = Check("paired", "Paired samples", CONFIRMED, f"{len(both)} individuals have both {a_name} and {b_name} samples", source=f"SDRF characteristics[{pd.info.get('individual_column') or 'individual'}]")
            elif inds_a or inds_b:
                ck = Check("paired", "Paired samples", EXCLUDED, f"individuals annotated ({len(inds_a | inds_b)}), but none has both {a_name} and {b_name} samples")
            else:
                ck = Check("paired", "Paired samples", UNKNOWN, "individuals not annotated, so pairing cannot be confirmed")
                stems_a = {_name_stem(s.name) for s in ga}
                stems_b = {_name_stem(s.name) for s in gb}
                shared = {x for x in stems_a & stems_b if x}
                if len(shared) >= 2:
                    ex = sorted(shared)[0]
                    ck.hint = f"sample names suggest {len(shared)} pairs (e.g. stem {_q(ex)}); check against the publication"
            checks.append(ck)
    elif c.diseases and c.comparison is None:
        pass

    # ---- counts --------------------------------------------------------------------
    cohort = pool
    inds = {s.individual for s in cohort if s.individual}
    has_ind_col = bool(pd.info.get("individual_column"))
    if c.min_individuals:
        n = c.min_individuals
        comparison_ok = any(x.id == "comparison" and x.status == CONFIRMED for x in checks)
        if c.paired and groups and comparison_ok:
            a_name, b_name = GROUP_NAMES[c.comparison]
            ia = {s.individual for s in groups[a_name] if s.individual}
            ib = {s.individual for s in groups[b_name] if s.individual}
            count, unit = len(ia & ib), "paired individuals"
            upper = min(len(groups[a_name]), len(groups[b_name]))
        else:
            count, unit = len(inds), "individuals"
            upper = len(cohort)
        if inds:
            st = CONFIRMED if count >= n else EXCLUDED
            ck = Check("min_individuals", f">= {n} individuals", st, f"{count} {unit} annotated", source=f"SDRF characteristics[{pd.info.get('individual_column')}]")
        elif upper < n:
            ck = Check("min_individuals", f">= {n} individuals", EXCLUDED, f"only {upper} biological samples in the matching groups, so at most {upper} {unit}")
        else:
            ck = Check("min_individuals", f">= {n} individuals", UNKNOWN, f"individuals not annotated{'' if has_ind_col else ' (no individual column)'}; {upper} biological samples is the upper bound")
        checks.append(ck)
    if c.min_samples:
        st = CONFIRMED if len(cohort) >= c.min_samples else EXCLUDED
        checks.append(Check("min_samples", f">= {c.min_samples} samples", st, f"{len(cohort)} biological samples match (technical/pooled samples excluded)"))

    # ---- technology (judged on the cohort's runs) -----------------------------------------
    cohort_names = {s.name for s in cohort}
    runs = [a for a in pd.assays if a["source_name"] in cohort_names] or pd.assays
    if c.instruments:
        subs = [x for i in c.instruments for x in vocab.INSTRUMENT_QUERY.get(i, [vocab.norm(i)])]
        inst = Counter(a["instrument"] for a in runs if a["instrument"] and not vocab.is_missing(a["instrument"]))
        ok = Counter({k: v for k, v in inst.items() if any(x in vocab.norm(k) for x in subs)})
        if ok:
            checks.append(Check("instrument", "Instrument", CONFIRMED, f"{sum(ok.values())}/{len(runs)} runs on {_top(ok)}", source="SDRF comment[instrument]"))
        elif inst:
            checks.append(Check("instrument", "Instrument", EXCLUDED, f"acquired on {_top(inst)}", source="SDRF comment[instrument]"))
        else:
            checks.append(Check("instrument", "Instrument", UNKNOWN, "instrument not annotated"))
    if c.labels:
        lt = Counter(vocab.label_type(a["label"]) for a in runs if a["label"])
        lt.pop(None, None)
        ok = [x for x in c.labels if lt.get(x)]
        if not ok and "TMT" in c.labels and any(vocab.norm(m["name"]).startswith("tmt") for m in pd.modifications):
            ok = ["TMT"]
            lt["TMT (from searched modification)"] += 1
        labels = {a["label"] for a in runs if a["label"]}
        plex = vocab.tmt_plex(labels)
        if any(vocab.norm(m["name"]) in ("tmtpro", "tmtpro16plex") for m in pd.modifications):
            plex = "TMTpro (16/18-plex)"
        if ok:
            st, det = CONFIRMED, f"labels: {_top(lt)}" + (f"; {plex}, {len(labels)} channels" if plex else "")
            if c.tmt_plex and plex and plex != c.tmt_plex:
                st, det = EXCLUDED, f"{plex}, not {c.tmt_plex}"
            checks.append(Check("label", "Quantification", st, det, source="SDRF comment[label]"))
        elif lt:
            checks.append(Check("label", "Quantification", EXCLUDED, f"labels: {_top(lt)}", source="SDRF comment[label]"))
        else:
            checks.append(Check("label", "Quantification", UNKNOWN, "labelling not annotated"))
    if c.acquisition:
        aq = Counter(vocab.acquisition_type(a["acquisition"]) for a in runs if a["acquisition"])
        aq.pop(None, None)
        ok = Counter({k: v for k, v in aq.items() if k in c.acquisition})
        if ok:
            checks.append(Check("acquisition", "Acquisition", CONFIRMED, f"{sum(ok.values())}/{len(runs)} runs: {_top(ok)}", source="SDRF comment[proteomics data acquisition method]"))
        elif aq:
            checks.append(Check("acquisition", "Acquisition", EXCLUDED, f"acquisition: {_top(aq)}"))
        else:
            checks.append(Check("acquisition", "Acquisition", UNKNOWN, "acquisition method not annotated"))

    if c.ptms:
        checks.append(_ptm_check(pd, cohort or bio, c.ptms))

    # ---- data availability -------------------------------------------------------------
    if c.raw_files:
        fmts = Counter(a["file_format"] for a in runs if a["data_file"])
        raw = sum(v for k, v in fmts.items() if k in vocab.RAW_FORMATS | vocab.OPEN_FORMATS)
        files = {a["data_file"] for a in runs if a["data_file"]}
        uri = sum(1 for a in runs if a["file_uri"])
        in_pride = sum(1 for f in files if f in pd.pride_files)
        where = (f"{in_pride}/{len(files)} found in PRIDE file listing" if pd.pride_files
                 else (f"direct URLs for {uri}/{len(runs)} rows" if uri else "download URLs resolvable via the PRIDE file API"))
        if files and raw == sum(fmts.values()):
            checks.append(Check("raw_files", "Raw files mapped", CONFIRMED, f"{len(files)} files ({_top(fmts, 2)}) mapped to samples; {where}", source="SDRF comment[data file]"))
        elif raw:
            checks.append(Check("raw_files", "Raw files mapped", CONFIRMED, f"{raw}/{sum(fmts.values())} rows map to spectra files ({_top(fmts, 3)}); {where}", source="SDRF comment[data file]",
                                hint="some rows point to archives or tables"))
        elif files:
            checks.append(Check("raw_files", "Raw files mapped", UNKNOWN, f"data files are {_top(fmts, 2)}, not vendor raw/mzML"))
        else:
            checks.append(Check("raw_files", "Raw files mapped", EXCLUDED, "no data files in SDRF"))
    if c.peak_lists:
        fmts = Counter(a["file_format"] for a in runs)
        open_n = sum(v for k, v in fmts.items() if k in vocab.OPEN_FORMATS)
        peak = [f for f in pd.pride_files.values() if (f.get("category") or "").upper() == "PEAK"]
        if open_n:
            checks.append(Check("peak_lists", "Peak lists", CONFIRMED, f"{open_n} rows map to open-format spectra ({_top(Counter({k: v for k, v in fmts.items() if k in vocab.OPEN_FORMATS}))})"))
        elif peak:
            checks.append(Check("peak_lists", "Peak lists", CONFIRMED, f"{len(peak)} PEAK files in PRIDE listing", source="PRIDE file API"))
        elif pd.pride_files:
            checks.append(Check("peak_lists", "Peak lists", EXCLUDED, "PRIDE file listing has no PEAK files", source="PRIDE file API"))
        else:
            checks.append(Check("peak_lists", "Peak lists", UNKNOWN, "SDRF maps vendor raw files only; PRIDE file listing not loaded (run `msagent enrich`)"))
    if c.published_after:
        year = None
        src = ""
        if pd.pride and pd.pride.get("publication_date"):
            year, src = pd.pride["publication_date"].year, "PRIDE publicationDate"
        elif pd.info.get("pub_year"):
            year, src = pd.info["pub_year"], "PRIDE FTP path in SDRF file URIs"
        if year is None:
            checks.append(Check("published_after", f"Published >= {c.published_after}", UNKNOWN, "publication date not in SDRF (run `msagent enrich` for PRIDE dates)"))
        else:
            st = CONFIRMED if year >= c.published_after else EXCLUDED
            checks.append(Check("published_after", f"Published >= {c.published_after}", st, f"published {year}", source=src))

    if c.exclude_flagged:
        bad = [k for k in pd.known_issues if k["severity"] in ("critical", "major")]
        other = [k for k in pd.known_issues if k["severity"] not in ("critical", "major")]
        if bad:
            k = bad[0]
            checks.append(Check("known_issues", "No flagged deposit issue", EXCLUDED, f"{k['severity']}: {k['section']} ({k['category']})", kind="safety", source="sdrf-annotated-datasets docs/known-issues.md"))
        else:
            ck = Check("known_issues", "No flagged deposit issue", CONFIRMED, "not listed in curators' known issues", kind="safety")
            if other:
                ck.detail = f"{other[0]['severity']} issue noted: {other[0]['section']}"
                warnings.append(ck.detail)
            checks.append(ck)

    # ---- summary -----------------------------------------------------------------------
    required = [x for x in checks if x.kind != "ranking"]
    n_exc = sum(1 for x in required if x.status == EXCLUDED)
    n_unk = sum(1 for x in required if x.status == UNKNOWN)
    status = "excluded" if n_exc else ("needs review" if n_unk else "qualifies")
    if not cohort and status != "excluded":
        status = "excluded"
    cohort_runs = [a for a in pd.assays if a["source_name"] in cohort_names]
    files = {a["data_file"] for a in cohort_runs if a["data_file"]}
    group_counts = {g: len(v) for g, v in groups.items()} if groups else {}
    technical = [s for s in samples_all if s.technical]
    counts = {
        "samples_total": len(samples_all),
        "technical_samples": len(technical),
        "cohort_samples": len(cohort),
        "individuals": len(inds) if inds else None,
        "individuals_annotated": bool(inds),
        "runs": len(files),
        "assay_rows": len(cohort_runs),
        "max_fraction": max((int(a["fraction"]) for a in cohort_runs if a["fraction"] and a["fraction"].isdigit() and int(a["fraction"]) <= 1000), default=None),
        "label_channels": len({a["label"] for a in cohort_runs if a["label"] and vocab.label_type(a["label"]) not in (None, "label-free")}) or None,
        "groups": group_counts,
    }
    if technical:
        warnings.append(f"{len(technical)} technical/pooled samples excluded from the cohort (e.g. {technical[0].technical})")
    return {
        "status": status,
        "checks": [x.to_dict() for x in checks],
        "n_confirmed": sum(1 for x in required if x.status == CONFIRMED),
        "n_unknown": n_unk,
        "n_excluded": n_exc,
        "counts": counts,
        "design": design_summary(counts),
        "cohort": [
            {"sample": s.name, "group": s.group, "group_evidence": s.group_evidence, "individual": s.individual}
            for s in cohort
        ],
        "warnings": warnings,
    }


def _ptm_check(pd: ProjectData, samples: list[Sample], ptms: list[str]) -> Check:
    want_any = "any" in ptms
    enriched = Counter()
    for s in samples:
        for a, v in s.vals(("enrichment process", "enrichment", "antibody enrichment")):
            cls = vocab.enrichment_ptm(v)
            if cls and (want_any or cls in ptms):
                enriched[(cls, v)] += 1
    searched = {m["ptm_class"]: m for m in pd.modifications if m["ptm_class"] and (want_any or m["ptm_class"] in ptms)}
    protocol_hit = None
    if pd.pride:
        text = vocab.norm(" ".join(str(pd.pride.get(k) or "") for k in ("sample_protocol", "description", "title")))
        for kw, cls in [("tio2", "phosphorylation"), ("imac", "phosphorylation"), ("phosphopeptide enrichment", "phosphorylation"),
                        ("phosphoenrichment", "phosphorylation"), ("k-ε-gg", "ubiquitination"), ("diglycine", "ubiquitination"),
                        ("ptmscan", None), ("glycopeptide enrichment", "glycosylation"), ("acetyl-lysine", "acetylation")]:
            if kw in text and (want_any or cls in ptms or cls is None):
                protocol_hit = kw
                break
    label = "PTM focus"
    if enriched:
        (cls, v), n = enriched.most_common(1)[0]
        return Check("ptm", label, CONFIRMED, f"{cls}: enrichment documented ({_q(v)}, {sum(enriched.values())} samples)", source="SDRF characteristics[enrichment process]")
    if protocol_hit:
        return Check("ptm", label, CONFIRMED, f"sample processing protocol mentions {_q(protocol_hit)}", source="PRIDE sampleProcessingProtocol")
    if searched:
        m = next(iter(searched.values()))
        return Check("ptm", label, UNKNOWN, f"{m['ptm_class']} searched as {m['mod_type'] or 'a'} modification ({m['name']}), but no enrichment step annotated",
                     source="SDRF comment[modification parameters]", hint="could be a global proteome search with PTMs allowed")
    if pd.modifications:
        names = ", ".join(sorted({m["name"] for m in pd.modifications})[:5])
        return Check("ptm", label, EXCLUDED, f"searched modifications: {names}", source="SDRF comment[modification parameters]")
    return Check("ptm", label, UNKNOWN, "modification parameters not annotated")


def design_summary(counts: dict) -> str:
    parts = []
    g = counts.get("groups") or {}
    if g and any(g.values()):
        parts.append(" + ".join(f"{v} {k}" for k, v in g.items()) + " samples")
    else:
        parts.append(f"{counts['cohort_samples']} samples")
    if counts.get("individuals"):
        parts.append(f"from {counts['individuals']} individuals")
    else:
        parts.append("(individuals not annotated)")
    run = f"-> {counts['runs']} data files"
    extra = []
    if counts.get("label_channels"):
        extra.append(f"{counts['label_channels']} label channels")
    if counts.get("max_fraction") and counts["max_fraction"] > 1:
        extra.append(f"up to {counts['max_fraction']} fractions")
    if extra:
        run += " (" + ", ".join(extra) + ")"
    parts.append(run)
    return " ".join(parts)


def classify_samples(pd: ProjectData) -> list[Sample]:
    """Sample views with technical/cell-line flags and every role we can infer (for drill-down)."""
    views = _sample_views(pd)
    for s in views:
        roles = {}
        for name, fn in (("disease_role", _role_disease), ("treatment_role", _role_treatment), ("genotype_role", _role_genotype)):
            r, ev = fn(s, pd.factor_attrs)
            if r:
                roles[name] = (r, ev)
        s.roles = roles  # type: ignore[attr-defined]
    return views


def group_lookup(samples: list[Sample]) -> dict[str, dict]:
    out: dict[str, dict] = defaultdict(dict)
    for s in samples:
        for k, (r, ev) in getattr(s, "roles", {}).items():
            out[s.name][k] = r
    return out
