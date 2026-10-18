"""A small synthetic SDRF corpus (clearly fake accessions PXD99xxxx) built once per test session."""

from __future__ import annotations

import pytest

HEADER = [
    "source name", "characteristics[organism]", "characteristics[organism part]", "characteristics[disease]",
    "characteristics[individual]", "characteristics[cell line]", "characteristics[treatment]",
    "characteristics[genotype]", "characteristics[enrichment process]", "characteristics[biological replicate]",
    "assay name", "comment[data file]", "comment[file uri]", "comment[fraction identifier]", "comment[technical replicate]",
    "comment[label]", "comment[instrument]", "comment[proteomics data acquisition method]",
    "comment[modification parameters]", "comment[modification parameters]", "factor value[disease]",
]
NA = "not applicable"
LF = "label free sample"
QE = "AC=MS:1002523;NT=Q Exactive HF"
DDA = "NT=data-dependent acquisition;AC=PRIDE:0000627"
DIA = "NT=data-independent acquisition;AC=PRIDE:0000450"
CAM = "NT=Carbamidomethyl;AC=UNIMOD:4;TA=C;MT=Fixed"
OX = "NT=Oxidation;AC=UNIMOD:35;TA=M;MT=Variable"
PHOS = "NT=Phospho;AC=UNIMOD:21;TA=S,T,Y;MT=Variable"


def row(src, org, part, dz, ind, cl, treat, geno, enr, rep, assay, f, uri, frac, label, inst, acq, m1, m2, fdz):
    return [src, org, part, dz, ind, cl, treat, geno, enr, rep, assay, f, uri, frac, "1", label, inst, acq, m1, m2, fdz]


def corpus_rows():
    data = {}
    # 1: paired HCC tumour / adjacent normal, 3 patients, direct URIs
    rows = []
    for p in (1, 2, 3):
        for kind, dz in (("T", "hepatocellular carcinoma"), ("N", "normal")):
            f = f"P{p}_{kind}.raw"
            rows.append(row(f"patient{p}_{kind}", "homo sapiens", "liver", dz, f"patient {p}", NA, NA, NA, NA, "1",
                            f"run_P{p}{kind}", f, f"https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/05/PXD990001/{f}",
                            "1", LF, QE, DDA, CAM, OX, dz))
    # patient 4 has a tumour sample only: must not enter a paired cohort
    rows.append(row("patient4_T", "homo sapiens", "liver", "hepatocellular carcinoma", "patient 4", NA, NA, NA, NA, "1",
                    "run_P4T", "P4_T.raw", "https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/05/PXD990001/P4_T.raw",
                    "1", LF, QE, DDA, CAM, OX, "hepatocellular carcinoma"))
    data["PXD990001"] = rows
    # 2: liver tumours only, no normal tissue
    data["PXD990002"] = [
        row(f"tumour{i}", "homo sapiens", "liver", "hepatocellular carcinoma", f"donor {i}", NA, NA, NA, NA, "1",
            f"t{i}", f"t{i}.raw", "", "1", LF, QE, DDA, CAM, OX, "hepatocellular carcinoma")
        for i in range(1, 5)
    ]
    # 3: HeLa, TMT 10-plex, phospho enrichment, DMSO vs drug
    tmt = ["TMT126", "TMT127N", "TMT127C", "TMT128N", "TMT128C", "TMT129N", "TMT129C", "TMT130N", "TMT130C", "TMT131"]
    rows = []
    for i, lab in enumerate(tmt):
        treat = "DMSO" if i < 5 else "trametinib 100 nM"
        for frac in (1, 2):
            rows.append(row(f"hela_{i}", "homo sapiens", "cell culture", "cervical adenocarcinoma", NA, "HeLa", treat, NA,
                            "enrichment of phosphorylated Protein", str(i + 1), f"plex1_f{frac}_{lab}", f"plex1_f{frac}.raw", "",
                            str(frac), lab, "NT=Orbitrap Fusion Lumos;AC=MS:1002732", DDA, CAM, PHOS, "cervical adenocarcinoma"))
    data["PXD990003"] = rows
    # 4: mouse brain, knockout vs wild type, DIA on timsTOF
    rows = []
    for i in range(1, 7):
        geno = "wild type" if i <= 3 else "Snca knockout"
        rows.append(row(f"mouse{i}", "Mus musculus", "hippocampus", "normal", NA, NA, NA, geno, NA, str(i), f"m{i}",
                        f"m{i}.d", "", "1", LF, "NT=timsTOF Pro 2;AC=MS:1003230", DIA, CAM, OX, "normal"))
    data["PXD990004"] = rows
    # 5: timsTOF runs are DDA, DIA runs are on a Q Exactive: no single run is timsTOF + DIA
    data["PXD990005"] = [
        row("plasma1", "homo sapiens", "blood plasma", "normal", "donor 1", NA, NA, NA, NA, "1", "a1", "a1.d", "", "1", LF,
            "NT=timsTOF Pro;AC=MS:1003005", DDA, CAM, OX, "normal"),
        row("plasma2", "homo sapiens", "blood plasma", "normal", "donor 2", NA, NA, NA, NA, "1", "a2", "a2.raw", "", "1", LF,
            QE, DIA, CAM, OX, "normal"),
    ]
    return data


@pytest.fixture(scope="session")
def corpus(tmp_path_factory):
    root = tmp_path_factory.mktemp("corpus")
    for acc, rows in corpus_rows().items():
        d = root / "datasets" / acc
        d.mkdir(parents=True)
        with open(d / f"{acc}.sdrf.tsv", "w") as fh:
            fh.write("\t".join(HEADER) + "\n")
            for r in rows:
                fh.write("\t".join(r) + "\n")
    (root / "docs").mkdir()
    (root / "docs" / "known-issues.md").write_text(
        "# Known issues\n\n## Critical: data file is not instrument data\n\nCategory: `data-file-not-spectra`.\n\n"
        "| Dataset | file |\n|---|---|\n| PXD990002 | `t1.raw` |\n"
    )
    return root


@pytest.fixture(scope="session")
def store(corpus, tmp_path_factory):
    from msagent import index, ingest
    from msagent.store import Store

    out = tmp_path_factory.mktemp("build")
    db = out / "test.duckdb"
    ingest.build(sdrf_dir=corpus, db_path=db, verbose=False)
    index.build(db_path=db, index_dir=out / "index", verbose=False)
    return Store(db_path=db, index_dir=out / "index")
