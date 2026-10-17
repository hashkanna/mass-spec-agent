"""PRIDE enrichment parsing, against a record shaped like the v3 API (the live API is not called)."""

import duckdb

from msagent.ingest import SCHEMA
from msagent.pride_api import parse_project

RECORD = {
    "project": {
        "accession": "PXD990001",
        "title": "Paired HCC proteomes",
        "projectDescription": "Tumour and adjacent liver.",
        "sampleProcessingProtocol": "Phosphopeptides were enriched with TiO2.",
        "dataProcessingProtocol": "MaxQuant.",
        "keywords": ["HCC", "liver"],
        "submissionDate": "2021-03-01",
        "publicationDate": "2021-05-20",
        "references": [{"referenceLine": "x", "pubmedID": 123, "doi": "10.1/abc"}],
    },
    "files": [
        {"fileName": "P1_T.raw", "fileCategory": {"value": "RAW"}, "fileSizeBytes": 123,
         "publicFileLocations": [{"name": "FTP Protocol", "value": "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/05/PXD990001/P1_T.raw"}]},
        {"fileName": "P1_T.mgf", "fileCategory": {"value": "PEAK"}, "fileSizeBytes": 45, "publicFileLocations": []},
    ],
}


def test_parse_and_store(tmp_path):
    prow, frows = parse_project(RECORD)
    assert prow[0] == "PXD990001" and prow[5] == ["HCC", "liver"] and prow[8] == ["123"]
    assert frows[0][4] == "https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/05/PXD990001/P1_T.raw"
    assert frows[1][2] == "PEAK" and frows[1][4] is None
    con = duckdb.connect(str(tmp_path / "t.duckdb"))
    con.execute(SCHEMA)
    con.execute("INSERT INTO pride_projects VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", list(prow))
    con.executemany("INSERT INTO pride_files VALUES (?, ?, ?, ?, ?, ?)", frows)
    assert con.execute("SELECT year(publication_date) FROM pride_projects").fetchone()[0] == 2021
    assert con.execute("SELECT count(*) FROM pride_files WHERE category = 'PEAK'").fetchone()[0] == 1
