"""Normalisation of SDRF values and natural-language terms onto shared concepts.

Repository metadata uses inconsistent spellings for the same thing ("blood plasma" / "plasma",
"TMT127N" / "TMT 10-plex", "q exactive hf" / "Q Exactive HF"). Everything here is a small,
explicit, inspectable mapping; ontology accessions from the SDRF are kept alongside.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

from .sdrf import is_missing


def norm(s: str | None) -> str:
    if s is None:
        return ""
    s = s.replace("–", "-").replace("—", "-").replace("→", "->").replace("_", " ")
    return " ".join(s.strip().lower().split())


# --------------------------------------------------------------------------------------
# Organisms
# --------------------------------------------------------------------------------------

ORGANISM_SYNONYMS = {
    "human": "homo sapiens",
    "humans": "homo sapiens",
    "patient": "homo sapiens",
    "patients": "homo sapiens",
    "mouse": "mus musculus",
    "mice": "mus musculus",
    "murine": "mus musculus",
    "rat": "rattus norvegicus",
    "rats": "rattus norvegicus",
    "yeast": "saccharomyces cerevisiae",
    "budding yeast": "saccharomyces cerevisiae",
    "fission yeast": "schizosaccharomyces pombe",
    "e. coli": "escherichia coli",
    "e.coli": "escherichia coli",
    "arabidopsis": "arabidopsis thaliana",
    "zebrafish": "danio rerio",
    "drosophila": "drosophila melanogaster",
    "fruit fly": "drosophila melanogaster",
    "c. elegans": "caenorhabditis elegans",
    "worm": "caenorhabditis elegans",
    "pig": "sus scrofa",
    "porcine": "sus scrofa",
    "cow": "bos taurus",
    "bovine": "bos taurus",
    "cattle": "bos taurus",
    "chicken": "gallus gallus",
    "dog": "canis lupus familiaris",
    "canine": "canis lupus familiaris",
    "macaque": "macaca mulatta",
    "rhesus": "macaca mulatta",
    "malaria parasite": "plasmodium falciparum",
    "plasmodium": "plasmodium falciparum",
    "rice": "oryza sativa",
    "maize": "zea mays",
    "tomato": "solanum lycopersicum",
    "xenopus": "xenopus laevis",
    "tuberculosis bacterium": "mycobacterium tuberculosis",
    "gut metagenome": "human gut metagenome",
}


def normalize_organism(value: str) -> str:
    v = norm(value)
    v = re.sub(r"\s*\(.*?\)\s*$", "", v)  # "homo sapiens (human)"
    return ORGANISM_SYNONYMS.get(v, v)


# --------------------------------------------------------------------------------------
# Instruments
# --------------------------------------------------------------------------------------

# (substring in normalised name, family, vendor). First match wins, so order matters.
_INSTRUMENT_FAMILIES = [
    ("astral", "Orbitrap Astral", "Thermo"),
    ("exploris", "Orbitrap Exploris", "Thermo"),
    ("q exactive", "Q Exactive", "Thermo"),
    ("qexactive", "Q Exactive", "Thermo"),
    ("lumos", "Orbitrap Tribrid", "Thermo"),
    ("eclipse", "Orbitrap Tribrid", "Thermo"),
    ("ascend", "Orbitrap Tribrid", "Thermo"),
    ("orbitrap fusion", "Orbitrap Tribrid", "Thermo"),
    ("orbitrap elite", "LTQ Orbitrap", "Thermo"),
    ("orbitrap velos", "LTQ Orbitrap", "Thermo"),
    ("ltq orbitrap", "LTQ Orbitrap", "Thermo"),
    ("orbitrap", "Orbitrap (other)", "Thermo"),
    ("ltq ft", "LTQ FT", "Thermo"),
    ("ltq", "LTQ (ion trap)", "Thermo"),
    ("velos", "LTQ (ion trap)", "Thermo"),
    ("tsq", "TSQ (triple quad)", "Thermo"),
    ("stellar", "Stellar", "Thermo"),
    ("timstof", "timsTOF", "Bruker"),
    ("tims tof", "timsTOF", "Bruker"),
    ("impact", "Bruker QTOF", "Bruker"),
    ("maxis", "Bruker QTOF", "Bruker"),
    ("autoflex", "MALDI-TOF", "Bruker"),
    ("ultraflex", "MALDI-TOF", "Bruker"),
    ("rapiflex", "MALDI-TOF", "Bruker"),
    ("zenotof", "ZenoTOF", "SCIEX"),
    ("tripletof", "TripleTOF", "SCIEX"),
    ("triple tof", "TripleTOF", "SCIEX"),
    ("qtrap", "QTRAP", "SCIEX"),
    ("triple quad", "Triple Quad", "SCIEX"),
    ("api ", "SCIEX (other)", "SCIEX"),
    ("synapt", "Synapt", "Waters"),
    ("xevo", "Xevo", "Waters"),
    ("q-tof", "QTOF (other)", "Other"),
    ("qtof", "QTOF (other)", "Other"),
    ("6550", "Agilent QTOF", "Agilent"),
    ("6545", "Agilent QTOF", "Agilent"),
    ("agilent", "Agilent (other)", "Agilent"),
    ("olink", "Olink (affinity)", "Olink"),
    ("somascan", "SomaScan (affinity)", "SomaLogic"),
    ("somalogic", "SomaScan (affinity)", "SomaLogic"),
]

_INSTRUMENT_CANON = {
    "q exactive": "Q Exactive",
    "q exactive hf": "Q Exactive HF",
    "q exactive hf-x": "Q Exactive HF-X",
    "q exactive plus": "Q Exactive Plus",
    "orbitrap fusion": "Orbitrap Fusion",
    "orbitrap fusion lumos": "Orbitrap Fusion Lumos",
    "orbitrap eclipse": "Orbitrap Eclipse",
    "orbitrap exploris 480": "Orbitrap Exploris 480",
    "orbitrap exploris 240": "Orbitrap Exploris 240",
    "orbitrap exploris 120": "Orbitrap Exploris 120",
    "orbitrap astral": "Orbitrap Astral",
    "orbitrap astral zoom": "Orbitrap Astral Zoom",
    "ltq orbitrap": "LTQ Orbitrap",
    "ltq orbitrap velos": "LTQ Orbitrap Velos",
    "ltq orbitrap elite": "LTQ Orbitrap Elite",
    "orbitrap elite": "Orbitrap Elite",
    "ltq orbitrap xl": "LTQ Orbitrap XL",
    "timstof pro": "timsTOF Pro",
    "timstof pro 2": "timsTOF Pro 2",
    "timstof ht": "timsTOF HT",
    "timstof scp": "timsTOF SCP",
    "timstof ultra": "timsTOF Ultra",
    "timstof ultra 2": "timsTOF Ultra 2",
    "timstof": "timsTOF",
    "tripletof 5600": "TripleTOF 5600",
    "tripletof 5600+": "TripleTOF 5600+",
    "tripletof 6600": "TripleTOF 6600",
    "zenotof 7600": "ZenoTOF 7600",
}


def instrument_info(raw: str | None) -> tuple[str, str, str] | None:
    """Return (model, family, vendor) for an instrument name, or None if missing."""
    if raw is None or is_missing(raw):
        return None
    v = norm(raw)
    v = re.sub(r"^(thermo( fisher)?( scientific)?|bruker( daltonics)?|ab sciex|sciex|waters)\s+", "", v)
    if v in {"instrument model", "mass spectrometer"}:
        return None
    model = _INSTRUMENT_CANON.get(v, raw.strip())
    for key, family, vendor in _INSTRUMENT_FAMILIES:
        if key in v:
            return model, family, vendor
    return model, "Other", "Other"


# Query phrase -> substrings that a matching instrument model/family must contain.
INSTRUMENT_QUERY = {
    "orbitrap astral": ["astral"],
    "astral": ["astral"],
    "exploris": ["exploris"],
    "q exactive": ["q exactive"],
    "qe hf": ["q exactive hf"],
    "lumos": ["lumos"],
    "eclipse": ["eclipse"],
    "tribrid": ["fusion", "lumos", "eclipse", "ascend"],
    "fusion": ["fusion"],
    "orbitrap": ["orbitrap", "q exactive", "exploris", "astral", "lumos", "eclipse"],
    "timstof": ["timstof"],
    "tims tof": ["timstof"],
    "tripletof": ["tripletof", "triple tof"],
    "zenotof": ["zenotof"],
    "sciex": ["tripletof", "zenotof", "qtrap", "triple quad"],
    "bruker": ["timstof", "impact", "maxis", "flex"],
    "synapt": ["synapt"],
    "waters": ["synapt", "xevo"],
    "maldi": ["flex", "maldi"],
    "olink": ["olink"],
    "somascan": ["somascan", "somalogic"],
}


# --------------------------------------------------------------------------------------
# Labelling / quantification
# --------------------------------------------------------------------------------------


def label_type(label: str | None) -> str | None:
    if label is None or is_missing(label):
        return None
    v = norm(label)
    if "label free" in v or v in {"lfq", "label-free", "label free sample", "unlabeled", "unlabelled"}:
        return "label-free"
    if v.startswith("tmt") or "tmt" in v:
        return "TMT"
    if "itraq" in v:
        return "iTRAQ"
    if "silac" in v:
        return "SILAC"
    if "dimethyl" in v:
        return "dimethyl"
    if "18o" in v:
        return "18O"
    if "15n" in v:
        return "15N"
    if "isotope label" in v:
        return "isotope (unspecified)"
    return "other"


_TMT_CHANNEL = re.compile(r"tmt\s*(?:pro)?\s*-?\s*(\d{3})([nc]?)")


def tmt_plex(labels: set[str]) -> str | None:
    """Infer the TMT plex from the set of channel labels used in a project."""
    chans = set()
    for lab in labels:
        m = _TMT_CHANNEL.search(norm(lab))
        if m:
            chans.add(m.group(1) + m.group(2))
    if not chans:
        return None
    hi = max(int(c[:3]) for c in chans)
    if hi >= 134 or len(chans) > 11:
        return "TMTpro (16/18-plex)"
    if hi >= 131 or len(chans) > 6:
        return "TMT 10/11-plex"
    return "TMT 6-plex"


LABEL_QUERY = {
    "tmtpro": "TMT",
    "tmt": "TMT",
    "tandem mass tag": "TMT",
    "itraq": "iTRAQ",
    "silac": "SILAC",
    "label-free": "label-free",
    "label free": "label-free",
    "lfq": "label-free",
    "dimethyl label": "dimethyl",
}


# --------------------------------------------------------------------------------------
# Acquisition
# --------------------------------------------------------------------------------------


def acquisition_type(raw: str | None) -> str | None:
    if raw is None or is_missing(raw):
        return None
    v = norm(raw)
    if "independent" in v or "swath" in v or "diapasef" in v or "hdmse" in v or v == "dia" or "mse" == v:
        return "DIA"
    if "dependent" in v or v == "dda" or "ddapasef" in v:
        return "DDA"
    if "parallel reaction" in v or v == "prm":
        return "PRM"
    if "selected reaction" in v or "multiple reaction" in v or v in {"srm", "mrm"}:
        return "SRM/MRM"
    return "other"


ACQ_QUERY = {
    "data-independent": "DIA",
    "data independent": "DIA",
    "dia": "DIA",
    "swath": "DIA",
    "diapasef": "DIA",
    "data-dependent": "DDA",
    "data dependent": "DDA",
    "dda": "DDA",
    "prm": "PRM",
    "parallel reaction monitoring": "PRM",
    "srm": "SRM/MRM",
    "mrm": "SRM/MRM",
    "targeted": "PRM",
}


# --------------------------------------------------------------------------------------
# Modifications / PTMs
# --------------------------------------------------------------------------------------

ROUTINE_MODS = {
    "carbamidomethyl",
    "oxidation",
    "deamidated",
    "deamidation",
    "gln->pyro-glu",
    "glu->pyro-glu",
    "met-loss",
    "met-loss+acetyl",
    "propionamide",
    "methylthio",
    "carbamyl",
    "ammonia-loss",
    "cation:na",
    "amidated",
}


@dataclass
class Modification:
    name: str
    accession: str | None = None
    targets: str | None = None
    mod_type: str | None = None
    position: str | None = None

    @property
    def ptm_class(self) -> str | None:
        return ptm_class(self.name, self.targets, self.position, self.accession)


def parse_modification(cell_dict: dict[str, str]) -> Modification | None:
    name = cell_dict.get("NT", "").strip()
    if is_missing(name):
        return None
    return Modification(
        name=name,
        accession=cell_dict.get("AC"),
        targets=cell_dict.get("TA"),
        mod_type=cell_dict.get("MT"),
        position=cell_dict.get("PP"),
    )


def ptm_class(name: str, targets: str | None = None, position: str | None = None, accession: str | None = None) -> str | None:
    """Biologically interesting PTM class for a searched modification, or None for routine ones."""
    n = norm(name)
    acc = (accession or "").upper()
    t = (targets or "").upper()
    p = norm(position)
    if "phospho" in n or acc == "UNIMOD:21":
        return "phosphorylation"
    if n in {"gg", "glygly", "gly-gly", "diglycine"} or acc == "UNIMOD:121" or "ubiquit" in n:
        return "ubiquitination"
    if n in ROUTINE_MODS:
        if n.startswith("deamidat") and "R" in t and "N" not in t and "Q" not in t:
            return "citrullination"
        return None
    if n == "acetyl" or acc == "UNIMOD:1":
        if "K" in t and "n-term" not in p:
            return "acetylation"
        return None
    if n in {"methyl", "trimethyl"} or (n == "dimethyl" and "R" in t):
        return "methylation"
    if "hexnac" in n or n.startswith("hex(") or n == "hex" or "glyco" in n or "neuac" in n or "fucose" in n:
        return "glycosylation"
    for key, cls in [
        ("succinyl", "succinylation"),
        ("crotonyl", "crotonylation"),
        ("malonyl", "malonylation"),
        ("lactyl", "lactylation"),
        ("citrullin", "citrullination"),
        ("nitrosyl", "S-nitrosylation"),
        ("nitro", "nitration"),
        ("sulfo", "sulfation"),
        ("adp-ribosyl", "ADP-ribosylation"),
        ("palmitoyl", "lipidation"),
        ("myristoyl", "lipidation"),
        ("farnesyl", "lipidation"),
        ("sumo", "SUMOylation"),
        ("hydroxyl", "hydroxylation"),
        ("formyl", "formylation"),
        ("crotonaldehyde", "carbonylation"),
        ("hne", "carbonylation"),
    ]:
        if key in n:
            return cls
    return None


def enrichment_ptm(value: str) -> str | None:
    v = norm(value)
    if "phospho" in v or "imac" in v or "tio2" in v or "titanium dioxide" in v or "fe-nta" in v:
        return "phosphorylation"
    if "glyco" in v or "lectin" in v or "hilic" in v:
        return "glycosylation"
    if "ubiquit" in v or "diglycine" in v or "k-gg" in v or "k-e-gg" in v or "ggk" in v:
        return "ubiquitination"
    if "acetyl" in v:
        return "acetylation"
    if "succinyl" in v:
        return "succinylation"
    return None


PTM_QUERY = {
    "phosphoproteom": "phosphorylation",
    "phosphopeptide": "phosphorylation",
    "phosphoryl": "phosphorylation",
    "phospho": "phosphorylation",
    "ubiquitylom": "ubiquitination",
    "ubiquitinom": "ubiquitination",
    "ubiquitin": "ubiquitination",
    "digly": "ubiquitination",
    "k-gg": "ubiquitination",
    "acetylom": "acetylation",
    "acetylation": "acetylation",
    "glycoproteom": "glycosylation",
    "glycopeptide": "glycosylation",
    "glycosyl": "glycosylation",
    "glycan": "glycosylation",
    "methylation": "methylation",
    "methylom": "methylation",
    "succinyl": "succinylation",
    "citrullin": "citrullination",
    "lactyl": "lactylation",
    "sumo": "SUMOylation",
    "nitrosyl": "S-nitrosylation",
}
ANY_PTM_PHRASES = ["post-translational modification", "post translational modification", "ptms", "ptm "]


# --------------------------------------------------------------------------------------
# Tissues / organism parts
# --------------------------------------------------------------------------------------


@dataclass
class TermPattern:
    """A concept matched by substrings, with exclusions. Shown to the user as-is."""

    label: str
    include: list[str]
    exclude: list[str] = field(default_factory=list)

    def matches(self, value: str) -> bool:
        v = norm(value)
        if not v or is_missing(v):
            return False
        if any(x in v for x in self.exclude):
            return False
        return any(_contains(v, x) for x in self.include)


def _contains(v: str, needle: str) -> bool:
    # short needles must match whole words ("csf", "hcc"), longer ones may be prefixes ("hepat")
    if len(needle) <= 4:
        return re.search(r"(?<![a-z0-9])" + re.escape(needle) + r"(?![a-z0-9])", v) is not None
    return needle in v


TISSUES: dict[str, TermPattern] = {
    "plasma": TermPattern("blood plasma", ["plasma"], ["plasma membrane", "seminal plasma", "plasma cell", "blood plasma cell"]),
    "serum": TermPattern("blood serum", ["serum"], ["serum-free", "serum free", "fetal bovine serum"]),
    "blood": TermPattern("blood (incl. plasma/serum/PBMC)", ["blood", "plasma", "serum", "pbmc", "peripheral blood"], ["plasma membrane", "seminal plasma", "plasma cell"]),
    "cerebrospinal fluid": TermPattern("cerebrospinal fluid", ["cerebrospinal", "csf"]),
    "csf": TermPattern("cerebrospinal fluid", ["cerebrospinal", "csf"]),
    "urine": TermPattern("urine", ["urine", "urinary"]),
    "saliva": TermPattern("saliva", ["saliva", "salivary"]),
    "seminal": TermPattern("seminal plasma / semen", ["seminal", "semen", "sperm"]),
    "milk": TermPattern("milk", ["milk"]),
    "tear": TermPattern("tear fluid", ["tear"]),
    "liver": TermPattern("liver", ["liver", "hepat"], ["hepatocyte-like cell line"]),
    "brain": TermPattern(
        "brain (any region)",
        ["brain", "cerebr", "cortex", "hippocamp", "cerebell", "striatum", "thalamus", "amygdala", "substantia nigra",
         "midbrain", "putamen", "caudate", "frontal lobe", "temporal lobe", "prefrontal", "neocortex", "white matter",
         "grey matter", "gray matter", "spinal cord"],
        ["adrenal cortex", "renal cortex", "kidney cortex", "cortex of kidney", "bone cortex"],
    ),
    "kidney": TermPattern("kidney", ["kidney", "renal", "nephr", "glomerul"], ["adrenal"]),
    "heart": TermPattern("heart", ["heart", "cardiac", "ventric", "atrium", "atrial", "myocard"], ["cerebral ventricle"]),
    "lung": TermPattern("lung", ["lung", "pulmonary", "bronch", "alveol"]),
    "breast": TermPattern("breast / mammary", ["breast", "mammary"]),
    "colon": TermPattern("colon / colorectum", ["colon", "colorectal", "rectum", "rectal", "cecum", "caecum", "sigmoid"]),
    "intestine": TermPattern("intestine", ["intestin", "colon", "ileum", "jejunum", "duodenum", "gut", "cecum", "caecum"]),
    "stomach": TermPattern("stomach", ["stomach", "gastric"]),
    "pancreas": TermPattern("pancreas", ["pancrea", "islet"]),
    "muscle": TermPattern("skeletal muscle", ["muscle", "gastrocnemius", "quadriceps", "soleus", "tibialis", "vastus", "myotube"], ["smooth muscle", "heart muscle", "cardiac muscle"]),
    "skin": TermPattern("skin", ["skin", "epiderm", "dermis", "dermal"]),
    "bone marrow": TermPattern("bone marrow", ["bone marrow"]),
    "spleen": TermPattern("spleen", ["spleen", "splen"]),
    "lymph node": TermPattern("lymph node", ["lymph node"]),
    "thymus": TermPattern("thymus", ["thymus"]),
    "ovary": TermPattern("ovary", ["ovary", "ovarian"]),
    "testis": TermPattern("testis", ["testis", "testes", "testicular"]),
    "prostate": TermPattern("prostate", ["prostat"]),
    "uterus": TermPattern("uterus / endometrium", ["uter", "endometri", "myometri"]),
    "cervix": TermPattern("uterine cervix", ["cervix", "cervical"], ["cervical spinal", "cervical lymph"]),
    "placenta": TermPattern("placenta", ["placent"]),
    "adipose": TermPattern("adipose tissue", ["adipose", "fat"]),
    "bone": TermPattern("bone", ["bone", "osteo"], ["bone marrow"]),
    "cartilage": TermPattern("cartilage", ["cartilage", "chondro"]),
    "retina": TermPattern("retina / eye", ["retina", "eye", "lens", "cornea"]),
    "thyroid": TermPattern("thyroid", ["thyroid"]),
    "esophagus": TermPattern("esophagus", ["esophag", "oesophag"]),
    "bladder": TermPattern("urinary bladder", ["bladder"]),
    "feces": TermPattern("feces / stool", ["feces", "faeces", "fecal", "faecal", "stool"]),
    "platelet": TermPattern("platelets", ["platelet"]),
    "tumor": TermPattern("tumour tissue", ["tumor", "tumour"]),
    "leaf": TermPattern("leaf", ["leaf", "leaves"]),
    "root": TermPattern("root", ["root"]),
    "seed": TermPattern("seed", ["seed"]),
}
TISSUE_ALIASES = {
    "hepatic": "liver",
    "cardiac": "heart",
    "myocardium": "heart",
    "renal": "kidney",
    "pulmonary": "lung",
    "mammary": "breast",
    "colorectal": "colon",
    "gut": "intestine",
    "gastric": "stomach",
    "skeletal muscle": "muscle",
    "pbmc": "blood",
    "pbmcs": "blood",
    "whole blood": "blood",
    "spinal fluid": "csf",
    "stool": "feces",
    "faeces": "feces",
    "fecal": "feces",
    "endometrium": "uterus",
    "cortex": "brain",
    "hippocampus": "brain",
    "cerebellum": "brain",
    "neuronal tissue": "brain",
    "eye": "retina",
}
TISSUE_ATTRS = ("organism part", "tissue", "sampling site", "host body site", "body site", "sample source")


# --------------------------------------------------------------------------------------
# Disease
# --------------------------------------------------------------------------------------

_NORMAL_RE = re.compile(
    r"(?<!ab)\b(normal|healthy|disease[- ]?free|non[- ]?tumou?r(al)?|non[- ]?cancer(ous)?|non[- ]?malignant|"
    r"tumou?r[- ]adjacent|adjacent|para[- ]?(cancer(ous)?|tumou?r(al)?|carcinoma)|nat|control|"
    r"no disease|not diseased|non[- ]?diseased|non[- ]?failing|donor)\b"
)
_TUMOUR_RE = re.compile(
    r"(cancer|carcinoma|tumou?r|neoplas|sarcoma|leuka?emia|lymphoma|melanoma|glioma|glioblastoma|blastoma|"
    r"adenoma|myeloma|mesothelioma|metasta|malignan|hcc|nsclc|pdac|ccrcc|tnbc)"
)


def is_normal_value(value: str) -> bool:
    v = norm(value)
    if not v or is_missing(v):
        return False
    return bool(_NORMAL_RE.search(v))


def is_tumour_value(value: str) -> bool:
    v = norm(value)
    if not v or is_missing(v) or is_normal_value(v):
        return False
    return bool(_TUMOUR_RE.search(v))


DISEASES: dict[str, TermPattern] = {
    "hepatocellular carcinoma": TermPattern("hepatocellular carcinoma", ["hepatocellular carcinoma", "hepatocellular cancer", "hcc", "hepatoma", "liver cancer", "liver carcinoma"]),
    "liver cancer": TermPattern("liver cancer", ["hepatocellular carcinoma", "hepatocellular cancer", "hcc", "hepatoma", "liver cancer", "liver carcinoma", "cholangiocarcinoma", "hepatoblastoma"]),
    "liver tumour": TermPattern("liver tumour (malignant or benign)", ["hepatocellular", "hcc", "hepatoma", "liver cancer", "liver carcinoma", "cholangiocarcinoma", "hepatoblastoma", "liver tumo", "hepatic adenoma", "liver adenoma"]),
    "breast cancer": TermPattern("breast cancer", ["breast cancer", "breast carcinoma", "breast adenocarcinoma", "ductal carcinoma", "lobular carcinoma", "mammary carcinoma", "triple-negative breast", "tnbc", "breast tumo"]),
    "colorectal cancer": TermPattern("colorectal cancer", ["colorectal", "colon cancer", "colon carcinoma", "colon adenocarcinoma", "rectal cancer", "rectal adenocarcinoma", "rectum adenocarcinoma", "crc"]),
    "lung cancer": TermPattern("lung cancer", ["lung cancer", "lung adenocarcinoma", "lung carcinoma", "nsclc", "non-small cell", "small cell lung", "lung squamous", "lung non-small"]),
    "prostate cancer": TermPattern("prostate cancer", ["prostate cancer", "prostate carcinoma", "prostate adenocarcinoma", "prostatic"]),
    "ovarian cancer": TermPattern("ovarian cancer", ["ovarian cancer", "ovarian carcinoma", "ovarian adenocarcinoma", "ovarian high grade", "ovarian serous", "neoplasm of ovary", "ovary carcinoma", "hgsoc"]),
    "gastric cancer": TermPattern("gastric cancer", ["gastric cancer", "gastric carcinoma", "gastric adenocarcinoma", "stomach carcinoma", "stomach cancer", "stomach adenocarcinoma"]),
    "pancreatic cancer": TermPattern("pancreatic cancer", ["pancreatic cancer", "pancreatic ductal", "pancreatic adenocarcinoma", "pancreatic carcinoma", "pdac"]),
    "kidney cancer": TermPattern("kidney cancer", ["renal cell carcinoma", "kidney cancer", "ccrcc", "renal carcinoma", "kidney carcinoma"]),
    "glioblastoma": TermPattern("glioblastoma / glioma", ["glioblastoma", "glioma", "gbm", "astrocytoma"]),
    "leukemia": TermPattern("leukaemia", ["leukemia", "leukaemia", "aml", "cml", "acute lymphoblastic"]),
    "lymphoma": TermPattern("lymphoma", ["lymphoma", "dlbcl"]),
    "myeloma": TermPattern("multiple myeloma", ["myeloma"]),
    "melanoma": TermPattern("melanoma", ["melanoma"]),
    "endometrial cancer": TermPattern("endometrial cancer", ["endometrial carcinoma", "endometrial cancer", "endometrial adenocarcinoma", "uterine corpus"]),
    "head and neck cancer": TermPattern("head and neck cancer", ["head and neck", "hnscc", "oral squamous", "nasopharyngeal"]),
    "alzheimer": TermPattern("Alzheimer's disease", ["alzheimer"]),
    "parkinson": TermPattern("Parkinson's disease", ["parkinson"]),
    "huntington": TermPattern("Huntington's disease", ["huntington"]),
    "als": TermPattern("amyotrophic lateral sclerosis", ["amyotrophic", "als"]),
    "multiple sclerosis": TermPattern("multiple sclerosis", ["multiple sclerosis"]),
    "covid": TermPattern("COVID-19", ["covid", "sars-cov-2", "sars-cov2"]),
    "sepsis": TermPattern("sepsis", ["sepsis", "septic"]),
    "diabetes": TermPattern("diabetes", ["diabet"]),
    "type 2 diabetes": TermPattern("type 2 diabetes", ["type 2 diabetes", "type ii diabetes", "t2d"]),
    "type 1 diabetes": TermPattern("type 1 diabetes", ["type 1 diabetes", "type i diabetes", "t1d"]),
    "obesity": TermPattern("obesity", ["obes"]),
    "heart failure": TermPattern("heart failure", ["heart failure", "hfpef", "hfref", "cardiomyopathy"]),
    "cardiomyopathy": TermPattern("cardiomyopathy", ["cardiomyopathy"]),
    "ibd": TermPattern("inflammatory bowel disease", ["crohn", "ulcerative colitis", "inflammatory bowel"]),
    "crohn": TermPattern("Crohn's disease", ["crohn"]),
    "ulcerative colitis": TermPattern("ulcerative colitis", ["ulcerative colitis"]),
    "nafld": TermPattern("fatty liver disease (NAFLD/NASH/MASLD)", ["nafld", "nash", "masld", "mash", "fatty liver", "steatohepatitis", "steatosis"]),
    "cirrhosis": TermPattern("liver cirrhosis", ["cirrhosis"]),
    "fibrosis": TermPattern("fibrosis", ["fibrosis"]),
    "rheumatoid arthritis": TermPattern("rheumatoid arthritis", ["rheumatoid"]),
    "osteoarthritis": TermPattern("osteoarthritis", ["osteoarthritis"]),
    "lupus": TermPattern("lupus", ["lupus"]),
    "preeclampsia": TermPattern("preeclampsia", ["eclampsia"]),
    "tuberculosis": TermPattern("tuberculosis", ["tuberculosis"]),
    "malaria": TermPattern("malaria", ["malaria", "plasmodium"]),
    "muscular dystrophy": TermPattern("muscular dystrophy", ["muscular dystrophy", "duchenne"]),
    "cancer": TermPattern("any cancer", ["cancer", "carcinoma", "tumor", "tumour", "neoplas", "sarcoma", "leukemia", "leukaemia", "lymphoma", "melanoma", "glioma", "glioblastoma", "blastoma", "myeloma", "adenoma", "malignan"], ["non-tumor", "non-tumour", "non tumor", "non tumour"]),
}
DISEASE_ALIASES = {
    "hcc": "hepatocellular carcinoma",
    "liver tumour": "liver tumour",
    "liver tumor": "liver tumour",
    "liver tumours": "liver tumour",
    "liver tumors": "liver tumour",
    "hepatic cancer": "liver cancer",
    "crc": "colorectal cancer",
    "colon cancer": "colorectal cancer",
    "nsclc": "lung cancer",
    "lung adenocarcinoma": "lung cancer",
    "pdac": "pancreatic cancer",
    "rcc": "kidney cancer",
    "renal cancer": "kidney cancer",
    "gbm": "glioblastoma",
    "glioma": "glioblastoma",
    "aml": "leukemia",
    "leukaemia": "leukemia",
    "alzheimer's": "alzheimer",
    "alzheimers": "alzheimer",
    "parkinson's": "parkinson",
    "parkinsons": "parkinson",
    "covid-19": "covid",
    "sars-cov-2": "covid",
    "t2d": "type 2 diabetes",
    "nash": "nafld",
    "masld": "nafld",
    "fatty liver": "nafld",
    "inflammatory bowel disease": "ibd",
    "dlbcl": "lymphoma",
    "tumour": "cancer",
    "tumor": "cancer",
    "tumours": "cancer",
    "tumors": "cancer",
    "carcinoma": "cancer",
    "oncology": "cancer",
}
DISEASE_ATTRS = ("disease", "disease state", "condition", "phenotype", "host disease status", "diagnosis")


# --------------------------------------------------------------------------------------
# Cell lines
# --------------------------------------------------------------------------------------

KNOWN_CELL_LINES = [
    "hela", "hek293", "hek 293", "293t", "hek-293", "k562", "k-562", "jurkat", "mcf-7", "mcf7", "a549", "a-549",
    "u2os", "u-2 os", "hepg2", "hep-g2", "huh7", "huh-7", "hct116", "hct 116", "sh-sy5y", "thp-1", "u937",
    "lncap", "pc3", "pc-3", "du145", "mda-mb-231", "mda-mb-468", "t47d", "skbr3", "sk-br-3", "caco-2", "caco2",
    "ht-29", "ht29", "sw480", "sw620", "rko", "hl-60", "molm-13", "mv4-11", "raji", "ramos", "nb4", "kasumi",
    "imr-90", "bj fibroblast", "rpe-1", "rpe1", "nih3t3", "nih-3t3", "3t3", "c2c12", "raw 264.7", "raw264.7",
    "cho", "cos-7", "vero", "mcf10a", "mcf-10a", "u-87", "u87", "hap1", "ipsc", "ips cell", "oci-ly", "su-dhl",
]
CELL_LINE_MARKERS = ("cell line", "cell culture", "cultured cell", "cell suspension culture", "cellosaurus", "immortali")


def looks_like_cell_line(value: str) -> bool:
    v = norm(value)
    if not v or is_missing(v):
        return False
    if any(m in v for m in CELL_LINE_MARKERS):
        return True
    return any(_contains(v, c) for c in KNOWN_CELL_LINES)


# --------------------------------------------------------------------------------------
# Comparison groups (used when checking whether a study contains a required contrast)
# --------------------------------------------------------------------------------------

CONTROL_RE = re.compile(
    r"\b(control|ctrl|untreated|non[- ]?treated|vehicle|dmso|mock|placebo|no drug|unstimulated|baseline|"
    r"sham|naive|na[iï]ve|scrambled?|empty vector|non[- ]?targeting|sicontrol|si-?ctrl|parental)\b"
)
WILDTYPE_RE = re.compile(r"\b(wild[- ]?type|wildtype|wt|parental|control|\+/\+)\b")
MUTANT_RE = re.compile(
    r"(knock[- ]?out|knock[- ]?down|\bko\b|\bkd\b|mutant|mutation|deletion|delta|-/-|\+/-|overexpress|sirna|shrna|"
    r"crispr|transgenic|null|deficient|\bcre\b|fl/fl|g\d+[a-z]\b|\b[a-z]\d+[a-z]\b|∆|δ)"
)
TREATMENT_ATTRS = ("treatment", "compound", "drug", "stimulation", "spiked compound", "dose", "exposure", "condition", "growth condition", "experimental condition", "perturbation")
GENOTYPE_ATTRS = ("genotype", "genetic modification", "strain or breed", "phenotype", "condition", "experimental condition")
INDIVIDUAL_ATTRS = ("individual", "participant identifier", "host subject id", "patient", "patient id", "donor", "subject", "subject id", "case id")


def ext_format(data_file: str) -> str:
    f = data_file.lower().strip()
    for suffix, fmt in [
        (".raw", "Thermo RAW"),
        (".wiff2", "SCIEX WIFF"),
        (".wiff", "SCIEX WIFF"),
        (".wiff.scan", "SCIEX WIFF"),
        (".d.zip", "Bruker/Agilent .d"),
        (".d.tar", "Bruker/Agilent .d"),
        (".d", "Bruker/Agilent .d"),
        (".baf", "Bruker BAF"),
        (".mzml.gz", "mzML"),
        (".mzml", "mzML"),
        (".mzxml", "mzXML"),
        (".mgf", "MGF peak list"),
        (".raw.zip", "Thermo RAW"),
        (".yep", "Bruker"),
        (".scan", "SCIEX WIFF"),
    ]:
        if f.endswith(suffix):
            return fmt
    if f.endswith((".zip", ".gz", ".7z", ".rar", ".tar", ".tgz")):
        return "archive"
    if f.endswith((".csv", ".txt", ".tsv", ".xlsx", ".xls", ".parquet")):
        return "table (not spectra)"
    return "other"


RAW_FORMATS = {"Thermo RAW", "SCIEX WIFF", "Bruker/Agilent .d", "Bruker BAF", "Bruker"}
OPEN_FORMATS = {"mzML", "mzXML", "MGF peak list"}
