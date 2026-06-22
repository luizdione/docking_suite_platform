"""
auto_md_loop.py - Loop de auto-correcao offline para o pipeline MD.

OBJETIVO: rodar sem LLM. Usuario starta esse script, ele itera, escreve
relatorio. Quando o usuario voltar, le o relatorio.

Cada iteracao:
  1. Roda validate_system.py
  2. Se FAIL no FAD XML -> aplica fix_fad_xml.py --inplace
  3. Roda md_clean.py --fresh --ns 0.1 --tag autofix_NNN
  4. Parse logs para classificar resultado:
       SUCCESS   - DCD escrito + analise feita
       NAN_WARMUP - "Particle coordinate is NaN" nas primeiras steps
       BOND_ANOMALY - forensic detectou bond > 1.5x esperado
       TEMPLATE_ERROR - "No template found"
       TIMEOUT - tempo > MAX_RUN_SEC
       SYNTAX/IMPORT - erro de Python no script
       OTHER
  5. Para cada classe aplica patch correspondente (registry)
  6. Para apos N iteracoes consecutivas com mesmo erro (loop detectado)
     OU sucesso OU N_MAX iteracoes atingido

Uso:
    python auto_md_loop.py [--max-iter 10] [--ns 0.1] [--platform OpenCL]
                           [--report report.md]
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
AUTOFIX_DIR = Path(__file__).resolve().parent
MD_CLEAN = ROOT / "md_clean.py"
PY = sys.executable

DEFAULT_REPORT = ROOT / "auto_md_report.md"
DEFAULT_STATE = ROOT / "auto_md_state.json"
DEFAULT_LOG_DIR = ROOT / "auto_md_logs"


def now() -> str:
    return datetime.now().isoformat(timespec="seconds")


def classify_log(log_text: str) -> dict:
    """Pattern-match the run log to classify the outcome."""
    r = {"class": "UNKNOWN", "evidence": []}
    if "Particle coordinate is NaN" in log_text:
        r["class"] = "NAN"
        m = re.search(r"\[(\d+/\d+)\][^\n]*", log_text)
        if m:
            r["evidence"].append(f"stage {m.group(1)}")
        if "BOND ANOMALA" in log_text:
            r["class"] = "BOND_ANOMALY_THEN_NAN"
            r["evidence"].append("forensic flagged anomalous bonds before NaN")
    elif "No template found for residue" in log_text:
        r["class"] = "MISSING_TEMPLATE"
    elif re.search(r"DONE\.|md_dryrun.*\.dcd", log_text):
        r["class"] = "SUCCESS"
    elif "SyntaxError" in log_text or "ImportError" in log_text or "ModuleNotFoundError" in log_text:
        r["class"] = "SYNTAX_IMPORT"
    elif "TimeoutExpired" in log_text:
        r["class"] = "TIMEOUT"
    elif "Traceback" in log_text:
        r["class"] = "RUNTIME_ERROR"
        m = re.search(r"(\w+Error|Exception): .+$", log_text, re.MULTILINE)
        if m:
            r["evidence"].append(m.group(0))
    return r


def run_subproc(cmd: list[str], timeout: int) -> tuple[int, str]:
    """Run subprocess, capture combined output. Returns (rc, output)."""
    try:
        result = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, timeout=timeout,
        )
        return result.returncode, result.stdout
    except subprocess.TimeoutExpired as e:
        partial = (e.stdout or "")
        if isinstance(partial, bytes):
            partial = partial.decode("utf-8", errors="replace")
        return -1, partial + "\n[TIMEOUT: subprocess exceeded {0}s]".format(timeout)


def step_validate() -> tuple[bool, str]:
    rc, out = run_subproc([PY, str(AUTOFIX_DIR / "validate_system.py")], timeout=120)
    return (rc == 0), out


def step_fix_fad() -> tuple[bool, str]:
    rc, out = run_subproc(
        [PY, str(AUTOFIX_DIR / "fix_fad_xml.py"), "--inplace"],
        timeout=60,
    )
    return (rc == 0), out


def step_run_md(tag: str, ns: float, platform: str, max_sec: int) -> tuple[int, str]:
    cmd = [PY, "-u", str(MD_CLEAN), "--fresh",
           "--ns", str(ns), "--tag", tag,
           "--platform", platform, "--allow-fallback"]
    return run_subproc(cmd, timeout=max_sec)


def append_report(report_path: Path, txt: str) -> None:
    with report_path.open("a", encoding="utf-8") as fh:
        fh.write(txt + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--max-iter", type=int, default=8,
                    help="Max iteracoes antes de desistir")
    ap.add_argument("--ns", type=float, default=0.1,
                    help="ns de producao por tentativa (curto p/ validar fix)")
    ap.add_argument("--platform", default="OpenCL",
                    choices=["CUDA", "OpenCL", "CPU", "AUTO"])
    ap.add_argument("--max-run-sec", type=int, default=900,
                    help="Timeout do md_clean.py por iteracao (default 900s = 15 min)")
    ap.add_argument("--report", default=str(DEFAULT_REPORT))
    ap.add_argument("--log-dir", default=str(DEFAULT_LOG_DIR))
    ap.add_argument("--state", default=str(DEFAULT_STATE))
    args = ap.parse_args()

    report = Path(args.report)
    log_dir = Path(args.log_dir)
    state_file = Path(args.state)
    log_dir.mkdir(parents=True, exist_ok=True)
    report.write_text(f"# Auto-MD report\n\nStarted: {now()}\n\n",
                      encoding="utf-8")

    state = {"started_at": now(), "iterations": [], "status": "running"}

    prev_class = None
    same_class_count = 0
    for i in range(1, args.max_iter + 1):
        iter_log: dict = {"i": i, "ts": now()}
        append_report(report, f"\n## Iteracao {i} — {now()}\n")

        # 1. Validate
        ok, val_out = step_validate()
        iter_log["validate_ok"] = ok
        (log_dir / f"iter_{i:02d}_validate.log").write_text(val_out, encoding="utf-8")
        append_report(report, "**Validate:** " + ("OK" if ok else "FAIL") +
                      "\n\n```\n" + val_out[-1500:] + "\n```")

        # 2. If FAD XML wrong, fix it (only do this once)
        if not ok and "wrong_elements" in val_out and i == 1:
            append_report(report, "\n→ Aplicando fix_fad_xml.py --inplace\n")
            fok, fout = step_fix_fad()
            (log_dir / f"iter_{i:02d}_fix.log").write_text(fout, encoding="utf-8")
            append_report(report, "**Fix FAD XML:** " + ("OK" if fok else "FAIL") +
                          "\n\n```\n" + fout[-1500:] + "\n```")
            iter_log["applied_fix"] = "fix_fad_xml"
            # Re-validate
            ok, val_out = step_validate()
            iter_log["validate_after_fix_ok"] = ok
            (log_dir / f"iter_{i:02d}_validate_after.log").write_text(val_out, encoding="utf-8")

        if not ok:
            append_report(report,
                          "\n→ Validate ainda FAIL apos fixes — interrompendo loop.\n")
            iter_log["class"] = "VALIDATE_FAIL"
            state["iterations"].append(iter_log)
            state["status"] = "stopped_validate"
            break

        # 3. Run MD
        tag = f"autofix_{i:03d}"
        t0 = time.time()
        rc, md_out = step_run_md(tag, args.ns, args.platform, args.max_run_sec)
        dt = time.time() - t0
        (log_dir / f"iter_{i:02d}_md.log").write_text(md_out, encoding="utf-8")
        cls = classify_log(md_out)
        iter_log["md_rc"] = rc
        iter_log["md_dt_sec"] = round(dt, 1)
        iter_log["class"] = cls["class"]
        iter_log["evidence"] = cls["evidence"]
        append_report(
            report,
            f"\n**MD run (tag={tag}):** rc={rc}, dt={dt:.1f}s, class={cls['class']}\n"
            f"\n```\n{md_out[-2000:]}\n```\n",
        )

        # 4. Decide next step
        if cls["class"] == "SUCCESS":
            append_report(report, "\n✓ SUCCESS — pipeline rodou ate o fim.\n")
            iter_log["resolved"] = True
            state["iterations"].append(iter_log)
            state["status"] = "success"
            break

        # Loop detector
        if cls["class"] == prev_class:
            same_class_count += 1
        else:
            same_class_count = 1
            prev_class = cls["class"]
        iter_log["same_class_streak"] = same_class_count

        if same_class_count >= 3:
            append_report(
                report,
                f"\n⚠ Mesmo erro ({cls['class']}) 3x consecutivas — "
                "loop detectado, interrompendo.\n",
            )
            state["iterations"].append(iter_log)
            state["status"] = "stuck_loop"
            break

        state["iterations"].append(iter_log)

    else:
        state["status"] = "max_iter_reached"

    state["ended_at"] = now()
    state_file.write_text(json.dumps(state, indent=2), encoding="utf-8")

    append_report(report, f"\n\n## Final — {state['status']}\n")
    append_report(report, f"\nIteracoes: {len(state['iterations'])}\n")
    append_report(report, f"State JSON: {state_file}\n")
    append_report(report, f"Logs em: {log_dir}\n")
    return 0 if state.get("status") == "success" else 1


if __name__ == "__main__":
    sys.exit(main())
