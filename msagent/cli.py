"""Command-line interface: `msagent <command>`."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

from .config import DB_PATH, INDEX_DIR, SDRF_DIR, SDRF_REPO_URL

ICON = {"confirmed": "+", "excluded": "x", "unknown": "?"}


def _store():
    from .store import Store

    return Store()


def cmd_fetch(args):
    if (SDRF_DIR / ".git").exists():
        print(f"updating {SDRF_DIR}")
        subprocess.run(["git", "-C", str(SDRF_DIR), "pull", "--ff-only", "--depth", "1"], check=False)
    else:
        SDRF_DIR.parent.mkdir(parents=True, exist_ok=True)
        print(f"cloning {SDRF_REPO_URL} -> {SDRF_DIR}")
        subprocess.run(["git", "clone", "--depth", "1", SDRF_REPO_URL, str(SDRF_DIR)], check=True)


def cmd_ingest(args):
    from .ingest import build

    stats = build(limit=args.limit)
    print(json.dumps(stats, indent=1))


def cmd_index(args):
    from .index import build

    print(json.dumps(build(), indent=1))


def cmd_build(args):
    if not (SDRF_DIR / "datasets").exists():
        cmd_fetch(args)
    cmd_ingest(args)
    cmd_index(args)
    if not args.no_map:
        try:
            cmd_map(args)
        except ImportError:
            print("skipping project map (pip install 'msagent[map]')")


def cmd_enrich(args):
    from .pride_api import PrideUnavailable, enrich

    try:
        print(json.dumps(enrich(limit=args.limit), indent=1))
    except PrideUnavailable as exc:
        sys.exit(str(exc))
    print("re-run `msagent index` so titles/descriptions are searchable")


def cmd_criteria(args):
    from .criteria import parse_query

    c = parse_query(args.question)
    if args.json:
        print(json.dumps(c.to_dict(), indent=1))
        return
    for item in c.describe():
        print(f"  [{item['kind']}] {item['label']}")


def print_results(res: dict, show_checks: bool = True):
    print(f"\nQuestion: {res['query']}")
    print("Criteria:")
    for item in res["criteria_list"]:
        print(f"  - {item['label']}" + ("" if item["kind"] == "required" else f"  ({item['kind']})"))
    s = res["summary"]
    print(f"\n{s['facet_candidates']} candidates from {s['projects_indexed']} studies -> verified {s['verified']}: "
          f"{s['qualifies']} qualify, {s['needs_review']} need review, {s['excluded']} excluded  ({s['seconds']}s)\n")
    for r in res["results"]:
        marks = "".join(ICON[c["status"]] for c in r["checks"])
        print(f"{r['accession']:<14} {r['status'].upper():<13} [{marks}]  {r['design']}")
        if r.get("title"):
            print(f"{'':14} {r['title'][:100]}")
        if show_checks:
            for c in r["checks"]:
                if c["status"] != "confirmed" or c["id"] in ("comparison", "paired", "min_individuals", "ptm"):
                    hint = f"  [{c['hint']}]" if c.get("hint") else ""
                    print(f"{'':16}{ICON[c['status']]} {c['label']}: {c['detail']}{hint}")
        print(f"{'':16}{r['links'].get('pride', '')}")


def cmd_search(args):
    from .cohort import find_cohorts

    res = find_cohorts(_store(), args.question, limit=args.limit)
    if args.json:
        print(json.dumps(res, indent=1, default=str))
    else:
        print_results(res, show_checks=not args.brief)


def cmd_project(args):
    from .project import project_summary

    print(json.dumps(project_summary(_store(), args.accession.upper()), indent=1, default=str))


def cmd_files(args):
    from .project import sample_file_map

    res = sample_file_map(_store(), args.accession.upper(), args.where)
    if args.json:
        print(json.dumps(res, indent=1, default=str))
        return
    if "error" in res:
        sys.exit(res["error"])
    print(res["answer"])
    for n in res["notes"]:
        print(f"  note: {n}")
    for s in res["samples"]:
        roles = ", ".join(f"{v['role']} ({v['evidence']})" for v in s["roles"].values())
        flag = f"  [technical: {s['technical']}]" if s["technical"] else ""
        print(f"- {s['sample']}" + (f"  individual={s['individual']}" if s["individual"] else "") + (f"  {roles}" if roles else "") + flag)
        for f in s["files"]:
            extra = " ".join(x for x in [f"fraction {f['fraction']}" if f["fraction"] not in (None, "1") else "",
                                         f["label"] if f["label"] and "free" not in f["label"].lower() else ""] if x)
            print(f"    {f['data_file']}  {extra}  {f['url'] or ''}")


def cmd_manifest(args):
    from .manifest import build_manifest, to_json, to_tsv

    man = build_manifest(_store(), [a.upper() for a in args.accessions], args.query)
    out = to_json(man) if args.format == "json" else to_tsv(man)
    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        Path(args.output).write_text(out)
        print(f"wrote {man['n_rows']} rows ({man['n_files']} files) to {args.output}")
    else:
        sys.stdout.write(out)


def cmd_map(args):
    from .projectmap import build

    print(json.dumps(build(), indent=1))


def cmd_serve(args):
    import uvicorn

    uvicorn.run("msagent.web.app:app", host=args.host, port=args.port, reload=False)


def cmd_mcp(args):
    from .mcp_server import main

    main()


def cmd_bench(args):
    from .bench import run

    run(args.queries, args.output, pride_live=args.pride_live)


def cmd_ask(args):
    from .agent import ask

    print(ask(args.question))


def cmd_status(args):
    st = _store()
    print(json.dumps({"db": str(DB_PATH), "index": str(INDEX_DIR), **st.meta,
                      "text_index": st.text.meta if st.text else None}, indent=1))


def main(argv=None):
    p = argparse.ArgumentParser(prog="msagent", description="Mass-Spec Agent: verified cohort search over public proteomics metadata")
    sub = p.add_subparsers(dest="cmd", required=True)

    sub.add_parser("fetch", help="clone/update the curated SDRF corpus").set_defaults(fn=cmd_fetch)
    sp = sub.add_parser("ingest", help="build the DuckDB metadata graph from SDRF files")
    sp.add_argument("--limit", type=int)
    sp.set_defaults(fn=cmd_ingest)
    sub.add_parser("index", help="build keyword + embedding indexes").set_defaults(fn=cmd_index)
    sp = sub.add_parser("build", help="fetch (if needed) + ingest + index + map")
    sp.add_argument("--limit", type=int)
    sp.add_argument("--no-map", action="store_true")
    sp.set_defaults(fn=cmd_build)
    sp = sub.add_parser("enrich", help="add PRIDE titles, descriptions, dates and file listings (needs www.ebi.ac.uk)")
    sp.add_argument("--limit", type=int)
    sp.set_defaults(fn=cmd_enrich)
    sub.add_parser("map", help="compute the 2-D project map for the web UI").set_defaults(fn=cmd_map)
    sub.add_parser("status", help="show index provenance").set_defaults(fn=cmd_status)

    sp = sub.add_parser("criteria", help="show how a question is interpreted")
    sp.add_argument("question")
    sp.add_argument("--json", action="store_true")
    sp.set_defaults(fn=cmd_criteria)
    sp = sub.add_parser("search", help="find and verify studies for an experimental question")
    sp.add_argument("question")
    sp.add_argument("--limit", type=int, default=10)
    sp.add_argument("--json", action="store_true")
    sp.add_argument("--brief", action="store_true")
    sp.set_defaults(fn=cmd_search)
    sp = sub.add_parser("project", help="experimental design of one study")
    sp.add_argument("accession")
    sp.set_defaults(fn=cmd_project)
    sp = sub.add_parser("files", help="sample -> file map of one study, optionally filtered")
    sp.add_argument("accession")
    sp.add_argument("--where", help="e.g. 'wild-type controls', 'tumour', 'patient 3'")
    sp.add_argument("--json", action="store_true")
    sp.set_defaults(fn=cmd_files)
    sp = sub.add_parser("manifest", help="export a cohort manifest (TSV/JSON)")
    sp.add_argument("accessions", nargs="+")
    sp.add_argument("--query", "-q", default="", help="question used to select samples within each study")
    sp.add_argument("--format", choices=["tsv", "json"], default="tsv")
    sp.add_argument("--output", "-o")
    sp.set_defaults(fn=cmd_manifest)
    sp = sub.add_parser("serve", help="run the web UI")
    sp.add_argument("--host", default="127.0.0.1")
    sp.add_argument("--port", type=int, default=8000)
    sp.set_defaults(fn=cmd_serve)
    sub.add_parser("mcp", help="run the MCP server on stdio").set_defaults(fn=cmd_mcp)
    sp = sub.add_parser("bench", help="run the retrieval benchmark")
    sp.add_argument("--queries", default=str(Path(__file__).resolve().parent.parent / "benchmarks" / "queries.yaml"))
    sp.add_argument("--output", default=str(Path(__file__).resolve().parent.parent / "benchmarks" / "results.json"))
    sp.add_argument("--pride-live", action="store_true", help="also score PRIDE's own keyword search (needs www.ebi.ac.uk)")
    sp.set_defaults(fn=cmd_bench)
    sp = sub.add_parser("ask", help="answer a question with Claude driving the tools (needs Anthropic credentials)")
    sp.add_argument("question")
    sp.set_defaults(fn=cmd_ask)

    args = p.parse_args(argv)
    args.fn(args)


if __name__ == "__main__":
    main()
