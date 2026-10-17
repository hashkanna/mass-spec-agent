from datetime import date

from msagent import vocab
from msagent.criteria import parse_query
from msagent.sdrf import parse_cell, parse_header


def test_parse_cell_and_header():
    assert parse_cell("NT=Phospho;AC=UNIMOD:21;TA=S,T;MT=Variable") == {"NT": "Phospho", "AC": "UNIMOD:21", "TA": "S,T", "MT": "Variable"}
    assert parse_cell("liver") == {"NT": "liver"}
    assert parse_header("characteristics [Organism Part]") == ("characteristics", "organism part")
    assert parse_header("factor value[disease]") == ("factor", "disease")
    assert parse_header("source name") == ("core", "source name")


def test_instruments_labels_acquisition():
    assert vocab.instrument_info("orbitrap astral zoom")[1:] == ("Orbitrap Astral", "Thermo")
    assert vocab.instrument_info("q exactive hf")[0] == "Q Exactive HF"
    assert vocab.instrument_info("timsTOF Pro 2")[1] == "timsTOF"
    assert vocab.instrument_info("not available") is None
    assert vocab.label_type("TMT127N") == "TMT"
    assert vocab.label_type("label free sample") == "label-free"
    assert vocab.tmt_plex({"TMT126", "TMT134N"}) == "TMTpro (16/18-plex)"
    assert vocab.tmt_plex({"TMT126", "TMT131"}) == "TMT 10/11-plex"
    assert vocab.acquisition_type("SWATH MS") == "DIA"
    assert vocab.acquisition_type("Data-Dependent Acquisition") == "DDA"


def test_ptm_classes():
    assert vocab.ptm_class("Phospho", "S,T") == "phosphorylation"
    assert vocab.ptm_class("GG", "K") == "ubiquitination"
    assert vocab.ptm_class("Acetyl", "K") == "acetylation"
    assert vocab.ptm_class("Acetyl", None, "Protein N-term") is None
    assert vocab.ptm_class("Carbamidomethyl", "C") is None
    assert vocab.enrichment_ptm("enrichment of phosphorylated Protein") == "phosphorylation"


def test_tissue_and_disease_patterns():
    plasma = vocab.TISSUES["plasma"]
    assert plasma.matches("blood plasma")
    assert not plasma.matches("plasma membrane")
    assert not plasma.matches("seminal plasma")
    assert vocab.TISSUES["brain"].matches("hippocampus")
    assert not vocab.TISSUES["brain"].matches("adrenal cortex")
    assert vocab.is_normal_value("non-tumoral liver")
    assert vocab.is_normal_value("normal")
    assert not vocab.is_normal_value("abnormal karyotype")
    assert vocab.is_tumour_value("hepatocellular carcinoma")
    assert not vocab.is_tumour_value("non-tumor")
    assert vocab.DISEASES["liver cancer"].matches("hepatocellular carcinoma")
    assert not vocab.DISEASES["liver cancer"].matches("hepatocellular adenoma")
    assert vocab.DISEASES["liver tumour"].matches("hepatocellular adenoma")
    assert vocab.looks_like_cell_line("HEK293T cell")
    assert not vocab.looks_like_cell_line("hepatocyte")


def test_parse_hero_question():
    c = parse_query(
        "Find human liver tumour proteomics studies with paired non-tumour tissue, at least ten patients, "
        "accessible raw files. Exclude cell lines."
    )
    assert c.organisms == ["homo sapiens"]
    assert c.tissues == ["liver"]
    assert c.diseases == ["liver tumour"]
    assert c.comparison == "tumour_vs_normal"
    assert c.paired and c.min_individuals == 10
    assert c.exclude_cell_lines and c.raw_files
    ids = [x["id"] for x in c.describe()]
    assert "disease" not in ids and "comparison" in ids


def test_parse_misc_questions():
    c = parse_query("Orbitrap Astral human plasma post-translational modifications from the last two years, processed peak lists",
                    today=date(2026, 10, 17))
    assert c.instruments == ["orbitrap astral"] and c.tissues == ["plasma"]
    assert c.ptms == ["any"] and c.peak_lists and c.published_after == 2024
    c = parse_query("Alzheimer's CSF vs healthy controls, 20+ patients, DIA")
    assert c.diseases == ["alzheimer"] and c.tissues == ["csf"]
    assert c.comparison == "disease_vs_control" and c.min_individuals == 20 and c.acquisition == ["DIA"]
    c = parse_query("mouse brain knockout vs wild type")
    assert c.organisms == ["mus musculus"] and c.comparison == "mutant_vs_wildtype"
    c = parse_query("human plasma TMTpro phospho after drug treatment")
    assert c.labels == ["TMT"] and c.tmt_plex == "TMTpro (16/18-plex)" and c.ptms == ["phosphorylation"]
    assert c.comparison == "treated_vs_control"
    assert "rug" not in c.keywords
