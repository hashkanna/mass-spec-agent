"""Turn an experimental question into explicit, inspectable criteria.

The parser is deliberately rule-based: every criterion it emits can be shown to the user,
edited, and re-run. An LLM client (via MCP) can also fill ``Criteria`` directly.
"""

from __future__ import annotations

import re
from dataclasses import asdict, dataclass, field
from datetime import date

from . import vocab

NUMBER_WORDS = {
    "one": 1, "two": 2, "three": 3, "four": 4, "five": 5, "six": 6, "seven": 7, "eight": 8, "nine": 9,
    "ten": 10, "eleven": 11, "twelve": 12, "fifteen": 15, "twenty": 20, "thirty": 30, "forty": 40,
    "fifty": 50, "hundred": 100, "a hundred": 100,
}
_NUM = r"(\d+|" + "|".join(sorted(NUMBER_WORDS, key=len, reverse=True)) + r")"
_PEOPLE = r"(patients|individuals|donors|subjects|participants|people|persons|cases|humans|mice|animals)"

CELL_TYPES = {
    "t cell": vocab.TermPattern("T cell", ["t cell", "t-cell", "t lymphocyte", "cd4", "cd8"]),
    "b cell": vocab.TermPattern("B cell", ["b cell", "b-cell", "b lymphocyte"]),
    "nk cell": vocab.TermPattern("NK cell", ["natural killer", "nk cell"]),
    "macrophage": vocab.TermPattern("macrophage", ["macrophage"]),
    "monocyte": vocab.TermPattern("monocyte", ["monocyte"]),
    "neutrophil": vocab.TermPattern("neutrophil", ["neutrophil"]),
    "dendritic cell": vocab.TermPattern("dendritic cell", ["dendritic"]),
    "neuron": vocab.TermPattern("neuron", ["neuron", "neural cell"]),
    "astrocyte": vocab.TermPattern("astrocyte", ["astrocyte"]),
    "microglia": vocab.TermPattern("microglia", ["microglia"]),
    "hepatocyte": vocab.TermPattern("hepatocyte", ["hepatocyte"]),
    "fibroblast": vocab.TermPattern("fibroblast", ["fibroblast"]),
    "cardiomyocyte": vocab.TermPattern("cardiomyocyte", ["cardiomyocyte", "cardiac muscle cell"]),
    "endothelial cell": vocab.TermPattern("endothelial cell", ["endothelial"]),
    "epithelial cell": vocab.TermPattern("epithelial cell", ["epithelial"]),
    "keratinocyte": vocab.TermPattern("keratinocyte", ["keratinocyte"]),
    "stem cell": vocab.TermPattern("stem cell", ["stem cell", "ipsc", "pluripotent"]),
    "organoid": vocab.TermPattern("organoid", ["organoid"]),
    "platelet": vocab.TermPattern("platelet", ["platelet"]),
    "erythrocyte": vocab.TermPattern("erythrocyte", ["erythrocyte", "red blood cell"]),
    "sperm": vocab.TermPattern("sperm", ["sperm", "spermatozo"]),
    "oocyte": vocab.TermPattern("oocyte", ["oocyte"]),
    "pbmc": vocab.TermPattern("PBMC", ["peripheral blood mononuclear", "pbmc"]),
    "adipocyte": vocab.TermPattern("adipocyte", ["adipocyte"]),
}
CELL_TYPE_ALIASES = {"t cells": "t cell", "b cells": "b cell", "neurons": "neuron", "macrophages": "macrophage",
                     "monocytes": "monocyte", "neutrophils": "neutrophil", "hepatocytes": "hepatocyte",
                     "fibroblasts": "fibroblast", "cardiomyocytes": "cardiomyocyte", "organoids": "organoid",
                     "platelets": "platelet", "stem cells": "stem cell", "ipscs": "stem cell", "ipsc": "stem cell",
                     "pbmcs": "pbmc", "astrocytes": "astrocyte", "keratinocytes": "keratinocyte"}

STOPWORDS = set(
    """a an the and or of in on for with without from to by at as is are be that this these those which what find
    show me give get list want need looking look studies study dataset datasets data project projects experiment
    experiments proteomics proteomic proteome mass spectrometry ms samples sample files file available
    least more than enough each per their its into using use used any all some only also please direct links
    link table clean verified return between across where who whose have has had over under within last
    years year recent processed peak lists list raw paired matched versus vs exclude excluding not no
    patients patient individuals donors subjects participants tissue tissues human mouse vs. after before
    treatment treated control controls healthy normal non tumour tumor cancer cell lines line two three four five
    ten twenty ftp http accessible metadata map mapping enough give direct link links looking disease diseases
    cells cell type types knockout wild modifications post-translational focus focused about related""".split()
)


@dataclass
class Criteria:
    text: str = ""
    organisms: list[str] = field(default_factory=list)
    tissues: list[str] = field(default_factory=list)
    diseases: list[str] = field(default_factory=list)
    cell_types: list[str] = field(default_factory=list)
    cell_lines: list[str] = field(default_factory=list)
    instruments: list[str] = field(default_factory=list)
    labels: list[str] = field(default_factory=list)
    tmt_plex: str | None = None
    acquisition: list[str] = field(default_factory=list)
    ptms: list[str] = field(default_factory=list)
    comparison: str | None = None  # tumour_vs_normal | disease_vs_control | treated_vs_control | mutant_vs_wildtype
    paired: bool = False
    min_individuals: int | None = None
    min_samples: int | None = None
    exclude_cell_lines: bool = False
    only_cell_lines: bool = False
    raw_files: bool = False
    peak_lists: bool = False
    published_after: int | None = None
    exclude_flagged: bool = True
    keywords: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "Criteria":
        known = {f for f in cls.__dataclass_fields__}
        return cls(**{k: v for k, v in d.items() if k in known and v is not None})

    def describe(self) -> list[dict]:
        """Human-readable list of criteria, in the order they are checked."""
        out = []

        def add(cid, label, kind="required"):
            out.append({"id": cid, "label": label, "kind": kind})

        if self.organisms:
            add("organism", "Organism: " + " or ".join(self.organisms))
        if self.tissues:
            add("tissue", "Tissue / biofluid: " + " or ".join(_tissue_label(t) for t in self.tissues))
        if self.diseases and self.comparison not in ("tumour_vs_normal", "disease_vs_control"):
            add("disease", "Disease: " + " or ".join(_disease_label(d) for d in self.diseases))
        if self.cell_types:
            add("cell_type", "Cell type: " + " or ".join(self.cell_types))
        if self.cell_lines:
            add("cell_line", "Cell line: " + " or ".join(self.cell_lines))
        if self.comparison == "tumour_vs_normal":
            dz = " or ".join(_disease_label(d) for d in self.diseases) or "tumour"
            add("comparison", f"Contains {dz} AND non-tumour/normal samples")
        elif self.comparison == "disease_vs_control":
            dz = " or ".join(_disease_label(d) for d in self.diseases) or "disease"
            add("comparison", f"Contains {dz} AND healthy/control samples")
        elif self.comparison == "treated_vs_control":
            add("comparison", "Contains treated AND control/untreated samples")
        elif self.comparison == "mutant_vs_wildtype":
            add("comparison", "Contains mutant/knockout AND wild-type samples")
        if self.paired:
            add("paired", "Pairing documented (same individual in both groups)")
        if self.min_individuals:
            add("min_individuals", f"At least {self.min_individuals} independent individuals" + (" with paired samples" if self.paired else ""))
        if self.min_samples:
            add("min_samples", f"At least {self.min_samples} biological samples")
        if self.exclude_cell_lines:
            add("no_cell_lines", "Primary material, not cell lines")
        if self.only_cell_lines:
            add("only_cell_lines", "Cell-line material")
        if self.instruments:
            add("instrument", "Instrument: " + " or ".join(self.instruments))
        if self.labels:
            lab = " or ".join(self.labels) + (f" ({self.tmt_plex})" if self.tmt_plex else "")
            add("label", "Quantification: " + lab)
        if self.acquisition:
            add("acquisition", "Acquisition: " + " or ".join(self.acquisition))
        if self.ptms:
            add("ptm", "PTM focus (enrichment documented): " + " or ".join(self.ptms))
        if self.raw_files:
            add("raw_files", "Vendor raw files mapped to samples")
        if self.peak_lists:
            add("peak_lists", "Processed peak lists (mzML/MGF) available")
        if self.published_after:
            add("published_after", f"Published in or after {self.published_after}")
        if self.exclude_flagged:
            add("known_issues", "No curator-flagged critical/major deposit issue", kind="safety")
        if self.keywords:
            add("keywords", "Topic relevance: " + ", ".join(self.keywords), kind="ranking")
        return out

    def is_empty(self) -> bool:
        d = self.to_dict()
        d.pop("text")
        d.pop("exclude_flagged")
        return not any(v for v in d.values())


def _tissue_label(t: str) -> str:
    return vocab.TISSUES[t].label if t in vocab.TISSUES else t


def _disease_label(d: str) -> str:
    return vocab.DISEASES[d].label if d in vocab.DISEASES else d


def _find_phrases(text: str, phrases: dict[str, str]) -> list[tuple[str, int, int]]:
    """Longest-first, non-overlapping whole-word phrase matches -> (value, start, end)."""
    hits = []
    taken = [False] * len(text)
    for phrase in sorted(phrases, key=len, reverse=True):
        for m in re.finditer(r"(?<![a-z0-9])" + re.escape(phrase) + r"(?![a-z0-9])", text):
            if any(taken[m.start():m.end()]):
                continue
            for i in range(m.start(), m.end()):
                taken[i] = True
            hits.append((phrases[phrase], m.start(), m.end()))
    return sorted(hits, key=lambda h: h[1])


def _num(s: str) -> int:
    return int(s) if s.isdigit() else NUMBER_WORDS[s]


def _uniq(xs):
    seen, out = set(), []
    for x in xs:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


def parse_query(text: str, known_organisms: list[str] | None = None, today: date | None = None) -> Criteria:
    today = today or date.today()
    t = " " + vocab.norm(text).replace("’", "'") + " "
    c = Criteria(text=text)
    consumed: list[tuple[int, int]] = []

    # ---- organism -------------------------------------------------------------------
    org_phr = dict(vocab.ORGANISM_SYNONYMS)
    for o in known_organisms or []:
        o = vocab.norm(o)
        if len(o) > 5 and " " in o:
            org_phr[o] = o
    for v, s, e in _find_phrases(t, org_phr):
        c.organisms.append(v)
        consumed.append((s, e))
    c.organisms = _uniq(c.organisms)

    # ---- comparison & pairing ---------------------------------------------------------
    normal_words = r"(non[- ]?tumou?r(al)?|non[- ]?cancer(ous)?|non[- ]?malignant|(adjacent )?normal|adjacent|healthy|benign|control|nat)"
    if re.search(r"\b(paired|matched)\b", t):
        c.paired = True
    tumour_ctx = re.search(r"tumou?r|cancer|carcinoma|malignan|neoplas|hcc\b|nsclc|pdac", t)
    contrast = re.search(
        r"(with|and|vs\.?|versus|against|plus|paired|matched|or)\s+(paired |matched )?" + normal_words
        + r"|" + normal_words + r"\s+(vs\.?|versus|and)\s+(tumou?r|cancer|disease)|tumou?r[- ]normal",
        t,
    )
    if contrast:
        c.comparison = "tumour_vs_normal" if tumour_ctx else "disease_vs_control"
        consumed.append(contrast.span())
    if re.search(r"\b(treated|treatment|drug|compound|inhibitor|stimulat\w*|perturb\w*|exposure|exposed|vehicle|dmso|untreated)\b", t):
        if not c.comparison:
            c.comparison = "treated_vs_control"
    for m in re.finditer(r"\b(knock[- ]?outs?|knock[- ]?downs?|ko|kd|wild[- ]?type|wt|mutants?|deficient|treated|treatment|drug|untreated|vehicle|controls?)\b", t):
        consumed.append(m.span())
    if re.search(r"\b(knock[- ]?out|knock[- ]?down|ko|kd|wild[- ]?type|wt|mutant|mutants|deficient|null mice)\b", t):
        if not c.comparison or c.comparison == "treated_vs_control":
            c.comparison = "mutant_vs_wildtype"

    # ---- counts ----------------------------------------------------------------------
    for pat in [
        r"(at least|minimum of|min\.?|>=|≥|no fewer than)\s*" + _NUM + r"\s+(?:\w+\s+){0,2}?" + _PEOPLE,
        _NUM + r"\s*\+\s*" + _PEOPLE,
        _NUM + r"\s+or more\s+" + _PEOPLE,
        r"(more than|over|>)\s*" + _NUM + r"\s+(?:\w+\s+){0,2}?" + _PEOPLE,
    ]:
        m = re.search(pat, t)
        if m:
            groups = [g for g in m.groups() if g]
            nums = [g for g in groups if g.isdigit() or g in NUMBER_WORDS]
            n = _num(nums[0])
            if pat.startswith(r"(more than"):
                n += 1
            c.min_individuals = n
            consumed.append(m.span())
            break
    m = re.search(r"(at least|minimum of|>=|≥|more than|over)\s*" + _NUM + r"\s+(biological\s+)?samples", t)
    if m:
        c.min_samples = _num(m.group(2)) + (1 if m.group(1) in ("more than", "over") else 0)

    # ---- cell lines -----------------------------------------------------------------
    if re.search(r"(exclud\w*|no|not|without|except|avoid|rather than)\s+(any\s+)?(immortali[sz]ed\s+)?cell[- ]lines?", t) or re.search(
        r"\b(primary (human )?tissues?|patient[- ]derived tissues?|tissue only|clinical (samples|specimens|tissue))\b", t
    ):
        c.exclude_cell_lines = True
    elif re.search(r"\bcell[- ]lines?\b", t):
        c.only_cell_lines = True
    for cl in vocab.KNOWN_CELL_LINES:
        if re.search(r"(?<![a-z0-9])" + re.escape(cl) + r"(?![a-z0-9])", t):
            c.cell_lines.append(cl)
    c.cell_lines = _uniq(c.cell_lines)

    # ---- tissue / disease / cell type -----------------------------------------------
    tissue_phr = {k: k for k in vocab.TISSUES if k != "tumor"}
    tissue_phr.update(vocab.TISSUE_ALIASES)
    tissue_phr.update({"blood plasma": "plasma", "blood serum": "serum", "plasma proteome": "plasma"})
    for v, s, e in _find_phrases(t, tissue_phr):
        c.tissues.append(v)
        consumed.append((s, e))
    c.tissues = _uniq(c.tissues)

    disease_phr = {k: k for k in vocab.DISEASES}
    disease_phr.update(vocab.DISEASE_ALIASES)
    # tissue + tumour wording ("liver tumour", "breast tumours") -> specific cancer
    for tissue, dz in [("liver", "liver tumour"), ("breast", "breast cancer"), ("lung", "lung cancer"),
                       ("colon", "colorectal cancer"), ("colorectal", "colorectal cancer"), ("prostate", "prostate cancer"),
                       ("ovarian", "ovarian cancer"), ("gastric", "gastric cancer"), ("stomach", "gastric cancer"),
                       ("pancreatic", "pancreatic cancer"), ("kidney", "kidney cancer"), ("renal", "kidney cancer"),
                       ("brain", "glioblastoma"), ("endometrial", "endometrial cancer")]:
        for w in ("tumour", "tumor", "tumours", "tumors", "cancer", "carcinoma"):
            disease_phr[f"{tissue} {w}"] = "liver cancer" if tissue == "liver" and w in ("cancer", "carcinoma") else dz
    hits = _find_phrases(t, disease_phr)
    # don't read the comparison partner ("non-tumour", "normal") as a disease
    for v, s, e in hits:
        prefix = t[max(0, s - 5):s]
        if re.search(r"non[- ]?$", prefix):
            continue
        c.diseases.append(v)
        consumed.append((s, e))
    c.diseases = _uniq(c.diseases)
    if len(c.diseases) > 1 and "cancer" in c.diseases:
        c.diseases.remove("cancer")
    if c.paired and tumour_ctx and not c.comparison:
        c.comparison = "tumour_vs_normal"
    if c.comparison == "tumour_vs_normal" and not c.diseases:
        c.diseases = ["cancer"]

    ct_phr = {k: k for k in CELL_TYPES}
    ct_phr.update(CELL_TYPE_ALIASES)
    for v, s, e in _find_phrases(t, ct_phr):
        if v == "pbmc" and "blood" in c.tissues:
            continue
        c.cell_types.append(v)
        consumed.append((s, e))
    c.cell_types = _uniq(c.cell_types)

    # ---- technology -----------------------------------------------------------------
    for v, s, e in _find_phrases(t, {k: k for k in vocab.INSTRUMENT_QUERY}):
        c.instruments.append(v)
        consumed.append((s, e))
    c.instruments = _uniq(c.instruments)

    for v, s, e in _find_phrases(t, vocab.LABEL_QUERY):
        c.labels.append(v)
        consumed.append((s, e))
    c.labels = _uniq(c.labels)
    if "tmtpro" in t or re.search(r"\b(16|18)[- ]?plex\b", t):
        c.tmt_plex = "TMTpro (16/18-plex)"
        if "TMT" not in c.labels:
            c.labels.append("TMT")
    elif re.search(r"\b(10|11)[- ]?plex\b", t):
        c.tmt_plex = "TMT 10/11-plex"
    elif re.search(r"\b6[- ]?plex\b", t) and "itraq" not in t:
        c.tmt_plex = "TMT 6-plex"

    for v, s, e in _find_phrases(t, vocab.ACQ_QUERY):
        c.acquisition.append(v)
        consumed.append((s, e))
    c.acquisition = _uniq(c.acquisition)

    for key, cls in vocab.PTM_QUERY.items():
        m = re.search(r"(?<![a-z])" + re.escape(key) + r"[a-z\-]*", t)
        if m:
            c.ptms.append(cls)
            consumed.append(m.span())
    m = re.search(r"post[- ]?translational modifications?|\bptms?\b", t)
    if m:
        consumed.append(m.span())
        if not c.ptms:
            c.ptms = ["any"]
    c.ptms = _uniq(c.ptms)

    # ---- files & dates ----------------------------------------------------------------
    if re.search(r"\braw (files?|data)\b|\b(downloadable|reprocess\w*|re-?analy[sz]\w*)\b|\bfiles available\b", t):
        c.raw_files = True
    if re.search(r"peak[- ]?lists?|\bmzml\b|\bmgf\b|processed (spectra|files)", t):
        c.peak_lists = True
    m = re.search(r"(last|past|previous)\s+" + _NUM + r"\s+years?", t)
    if m:
        c.published_after = today.year - _num(m.group(2))
    else:
        m = re.search(r"(since|after|from|in or after)\s+(20\d\d)", t)
        if m:
            c.published_after = int(m.group(2)) + (1 if m.group(1) == "after" else 0)
        elif re.search(r"\b(last|past) year\b", t):
            c.published_after = today.year - 1

    # ---- residual keywords (ranking only) ----------------------------------------------
    mask = list(t)
    for s, e in consumed:
        for i in range(max(0, s), min(len(mask), e)):
            mask[i] = " "
    residual = re.findall(r"[a-z][a-z0-9\-]{2,}", "".join(mask))
    c.keywords = _uniq([w for w in residual if w not in STOPWORDS])[:8]
    return c
