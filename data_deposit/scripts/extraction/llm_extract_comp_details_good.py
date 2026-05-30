#!/usr/bin/env python3
# llm_extract_comp_details.py
"""
Batch‑extract computational‑chemistry details from text files (converted PDFs)
flagged by detect_cc_metadata.py.  Uses an LLM provider of your choice and
emits one JSON per source file.

Directory layout (default):
.
├── chem_metadata_flagged.txt   # produced earlier
└── comp_details/               # JSON output goes here
"""

from pathlib import Path
import textwrap
import sys
import json
import requests
from requests.exceptions import Timeout, RequestException

# ─── USER SETTINGS ─────────────────────────────────────────────────────────────
PROVIDER         = "openai"         # "openai" | "gemini" | "anthropic" | "deepseek" | "local"
MODEL            = "gpt-4o-mini"    # model for the chosen provider

#FLAG_LIST        = Path("chem_metadata_flagged.txt")   # input list
#OUT_DIR          = Path("comp_details")                # output folder
    
# API keys (fill only the ones you need; prefer env vars over committing keys)
# Example: OPENAI_KEY = os.environ.get("OPENAI_API_KEY", "")
OPENAI_KEY       = "sk-proj-REPLACE_ME"
GEMINI_KEY       = "GOAI-REPLACE_ME"
ANTHROPIC_KEY    = "sk-ant-REPLACE_ME"
DEEPSEEK_KEY     = "dpk-REPLACE_ME"
DEEPSEEK_BASEURL = "https://api.deepseek.com/v1"

# Local endpoint (Ollama / llama.cpp compatible)
LOCAL_URL        = "http://localhost:11434/api/generate"
LOCAL_STREAM     = False           # keep False – easier JSON handling

# Model‑specific rough character limits (heuristic)
MODEL_LIMITS = {
    "gpt-4o-mini":               200_000,
    "gpt-4o":                    200_000,
    "claude-3-sonnet-20240229":  120_000,
    "gemini-1.5-flash":          200_000,
    "deepseek-chat":             160_000,
    "qwen-chem:latest":          200_000,
}
DEFAULT_LIMIT   = 200_000
# ───────────────────────────────────────────────────────────────────────────────


def max_chars(model: str) -> int:
    return MODEL_LIMITS.get(model, DEFAULT_LIMIT)


# ─── PROMPT TEMPLATE ───────────────────────────────────────────────────────────
SYSTEM_MSG = (
    "You are a meticulous computational chemist. "
    "Given a plain‑text methods section, you extract structured details."
)
USER_TEMPLATE = textwrap.dedent("""
    **Methods text**
    ----------------
    {text}

    **Task**
    1. Read the text above and identify computational‑chemistry details.
    2. Produce a JSON object with these keys (omit any that are not present):
       - "software":            e.g. "Gaussian 16", "ORCA 5.0.4"
       - "level_of_theory":     full string, e.g. "B3LYP‑D3(BJ)/6‑31G(d,p)", can be one or multiple
       - "basis_set":           e.g. "6‑31G(d,p)", "def2‑TZVP", can be one or multiple
       - "functional":          e.g. "B3LYP‑D3(BJ)", "PBE0", can be one or multiple
       - "method_family":       one or multiple of ["DFT", "HF", "MP2", "CCSD", "CCSD(T)",
                                        "TD‑DFT", "semi‑empirical", ...]
       - "calculation_types":   what was actually computed, choose (can be one or multiple) from
                                ["geometry_optimization",
                                 "transition_state_optimization",
                                 "single_point",
                                 "frequency_analysis"
                                 "vibrational_analysis",
                                 "IRC",
                                 "NMR",
                                 "UV-Vis",
                                 "thermochemistry",
                                 "MD",
                                 "other"]
       - "special_treatments":  solvent model, dispersion, relativistic, etc.
       - "comments":            free‑text notes worth keeping
    3. Output **only** a fenced ```json block, nothing else.
""").strip()

# ─── Provider helpers (mirrors the structure of html_to_biblatex-2.py) ─────────
def ask_openai(prompt: str) -> str:
    from openai import OpenAI
    client = OpenAI(api_key=OPENAI_KEY)
    resp = client.chat.completions.create(
        model=MODEL,
        temperature=0.2,
        messages=[
            {"role": "system", "content": SYSTEM_MSG},
            {"role": "user",   "content": prompt},
        ],
    )
    return resp.choices[0].message.content.strip()


def ask_deepseek(prompt: str) -> str:
    from openai import OpenAI
    client = OpenAI(api_key=DEEPSEEK_KEY, base_url=DEEPSEEK_BASEURL)
    resp = client.chat.completions.create(
        model=MODEL,
        temperature=0.2,
        messages=[
            {"role": "system", "content": SYSTEM_MSG},
            {"role": "user",   "content": prompt},
        ],
    )
    return resp.choices[0].message.content.strip()


def ask_anthropic(prompt: str) -> str:
    import anthropic
    client = anthropic.Anthropic(api_key=ANTHROPIC_KEY)
    resp = client.messages.create(
        model=MODEL,
        max_tokens=4096,
        temperature=0.2,
        system=SYSTEM_MSG,
        messages=[{"role": "user", "content": prompt}],
    )
    return resp.content[0].text.strip()


def ask_gemini(prompt: str) -> str:
    import google.generativeai as genai
    genai.configure(api_key=GEMINI_KEY)
    model = genai.GenerativeModel(MODEL)
    resp  = model.generate_content(prompt)
    return resp.text.strip()


def ask_local(prompt: str) -> str:
    payload = {"model": MODEL, "prompt": prompt, "stream": LOCAL_STREAM}
    r = requests.post(LOCAL_URL, json=payload, timeout=60)
    r.raise_for_status()
    return r.json()["response"].strip()


ASK_FN = {
    "openai":    ask_openai,
    "deepseek":  ask_deepseek,
    "anthropic": ask_anthropic,
    "gemini":    ask_gemini,
    "local":     ask_local,
}.get(PROVIDER)
# ───────────────────────────────────────────────────────────────────────────────


def extract_details(text: str) -> dict[str, str]:
    """Send text to the model and parse the fenced JSON block."""
    prompt   = USER_TEMPLATE.format(text=text[:max_chars(MODEL)])
    raw_resp = ASK_FN(prompt)

    # Expect model to wrap JSON in triple‑backtick ```json ... ```
    try:
        json_block = raw_resp.split("```json")[1].split("```")[0]
    except IndexError:
        raise ValueError("LLM response lacked a ```json block")

    try:
        return json.loads(json_block)
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid JSON from model: {exc}") from exc


def main(root: Path) -> None:
    if ASK_FN is None:
        sys.exit(f"Unsupported PROVIDER: {PROVIDER}")

    OUT_DIR = Path(f"{root}/comp_details")                # output folder
    OUT_DIR.mkdir(exist_ok=True)

    # every pdf has a chem_metadata_flagged.txt, loop over each , write a metadata extraction json for each
    for flag_path in root.glob("comp_details/*chem_metadata_flagged.txt"):

        FLAG_LIST = Path(flag_path)   # input list
        if not FLAG_LIST.is_file():
            sys.exit(f"Flag list not found: {FLAG_LIST}")

        prefix = flag_path.name.rsplit("chem_metadata_flagged.txt", 1)[0]
        out_path = OUT_DIR / f"{prefix}comp_detail.json"
        if out_path.is_file():
            print(f"{out_path} is already performed.")
            continue

        flagged_paths = [
            Path(line.strip()) for line in FLAG_LIST.read_text().splitlines() if line.strip()
        ]

        if not flagged_paths:
            print(f"No files listed in {FLAG_LIST}")
            continue

        # ── NEW: combine every file, one LLM call ────────────────────────────────────
        combined_parts = []
        for txt_path in flagged_paths:
            if not txt_path.is_file():
                print(f"[skip] Missing file: {txt_path}", file=sys.stderr)
                continue
            combined_parts.append(txt_path.read_text(encoding="utf-8", errors="ignore"))

        if not combined_parts:
            sys.exit("No readable files in chem_metadata_flagged.txt")
        print(f"• Combining {len(combined_parts)} files → single LLM call")

        try:
            comp_details = extract_details("\n\n".join(combined_parts))
        except Exception as exc:
            sys.exit(f"error: {exc}")
        out_path.write_text(json.dumps(comp_details, indent=2) + "\n", encoding="utf-8")
        print(out_path.relative_to(OUT_DIR.parent))
        print("\nDone! Extracted details saved in:", OUT_DIR.resolve())


if __name__ == "__main__":
    scan_root = Path(sys.argv[1]).expanduser().resolve() if len(sys.argv) > 1 else Path.cwd()
    main(scan_root)
