"""`msagent ask`: a Claude tool-use loop over the same tools the MCP server exposes.

Needs Anthropic credentials (ANTHROPIC_API_KEY or an `ant auth login` profile).
"""

from __future__ import annotations

import asyncio
import json

MODEL = "claude-opus-5-5"

SYSTEM = """You help proteomics researchers assemble cohorts of public mass-spectrometry data.
Use the tools to interpret the question into explicit criteria, find and verify candidate studies,
and inspect sample-to-file mappings. In your answer:
- lead with the studies that qualify, then those needing review, then notable exclusions with the reason;
- for each study give individuals vs samples vs data files, and quote the evidence for key criteria;
- state plainly which criteria are unknown (never present an unknown as confirmed);
- end with the next step (e.g. export_manifest for the chosen accessions).
Be concise; use a compact table where it helps."""


def _tools():
    from . import mcp_server

    listed = asyncio.run(mcp_server.mcp.list_tools())
    specs = [{"name": t.name, "description": t.description or "", "input_schema": t.inputSchema} for t in listed]
    funcs = {t.name: getattr(mcp_server, t.name) for t in listed}
    return specs, funcs


def ask(question: str, max_turns: int = 12, verbose: bool = True) -> str:
    import anthropic

    client = anthropic.Anthropic()
    tools, funcs = _tools()
    messages = [{"role": "user", "content": question}]
    for _ in range(max_turns):
        response = client.beta.messages.create(
            model=MODEL,
            max_tokens=16000,
            system=SYSTEM,
            tools=tools,
            messages=messages,
            output_config={"effort": "medium"},
            # on a safety decline, retry server-side on a fallback model
            betas=["server-side-fallback-2026-07-01"],
            fallbacks="default",
        )
        if response.stop_reason == "refusal":
            return "The model declined this request."
        if response.stop_reason != "tool_use":
            return "".join(b.text for b in response.content if b.type == "text")
        messages.append({"role": "assistant", "content": response.content})
        results = []
        for block in response.content:
            if block.type != "tool_use":
                continue
            if verbose:
                print(f"  -> {block.name}({json.dumps(block.input)[:160]})", flush=True)
            try:
                out = funcs[block.name](**block.input)
                results.append({"type": "tool_result", "tool_use_id": block.id, "content": json.dumps(out, default=str)[:60000]})
            except Exception as exc:  # report to the model rather than crash the loop
                results.append({"type": "tool_result", "tool_use_id": block.id, "content": f"error: {exc}", "is_error": True})
        messages.append({"role": "user", "content": results})
    return "Stopped after the maximum number of tool-use turns."
