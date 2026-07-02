#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""health.py - auditoria de manutencao e autoteste do bonder.

Reune, num relatorio unico, o estado de sanidade da plataforma:
  - dependencias Python (nucleo x opcionais)
  - ferramentas externas no PATH (SAINTexpress, Rscript, DIA-NN, etc.)
  - pacotes R gold-standard (DEP/limma/MSstats) quando Rscript existe
  - GPU (probe seguro) e VRAM livre
  - configuracao (iscas, contrastes, avisos de validacao)
  - bancos de referencia em cache + FASTA de busca
  - caminhos de dados/corridas (existencia + escrita)
  - corridas registradas e seus status

self_test() roda o pipeline ponta-a-ponta em dados sinteticos, num diretorio
temporario, e reporta o status de cada etapa - a checagem de sanidade dos
scripts e processos, sem depender de GPU nem de R.
"""
from __future__ import annotations

import importlib
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Dict, List

from bonder import __version__
from bonder.config import load_config


# Dependencias Python: (modulo_import, pacote_pip, nucleo?)
PY_DEPS = [
    ("numpy", "numpy", True),
    ("pandas", "pandas", True),
    ("pyarrow", "pyarrow", True),
    ("networkx", "networkx", True),
    ("yaml", "pyyaml", True),
    ("requests", "requests", True),
    ("flask", "flask", True),
    ("scipy", "scipy", False),          # declarado; nao usado no caminho central
    ("matplotlib", "matplotlib", False),  # so figuras (volcano); degrada com aviso
    ("pyteomics", "pyteomics", False),
    ("gprofiler", "gprofiler-official", False),
    ("rq", "rq", False),
    ("redis", "redis", False),
    ("torch", "torch", False),
    ("transformers", "transformers", False),
]

# Ferramentas externas: nome no PATH -> (proposito, obrigatoria?)
EXTERNAL_TOOLS = [
    ("SAINTexpress-int", "especificidade isca-presa (SAINT)", False),
    ("Rscript", "estatistica DEP/limma/MSstats (R)", False),
    ("diann", "busca DIA acelerada por GPU", False),
    ("ms2rescore", "rescoring por deep learning (GPU)", False),
    ("colabfold_batch", "validacao estrutural AF2-Multimer (GPU)", False),
    ("nvidia-smi", "GPU NVIDIA", False),
]

R_PACKAGES = ["DEP", "limma", "MSstats"]


def _pkg_version(pip_name: str) -> str:
    # importlib.metadata evita atributos __version__ deprecados (ex.: Flask 3.1+)
    try:
        from importlib.metadata import version, PackageNotFoundError
        try:
            return version(pip_name)
        except PackageNotFoundError:
            return "?"
    except Exception:
        return "?"


def check_python_deps() -> List[Dict[str, Any]]:
    out = []
    for imp, pip_name, core in PY_DEPS:
        try:
            importlib.import_module(imp)
            out.append({"name": pip_name, "import": imp, "ok": True,
                        "version": _pkg_version(pip_name), "core": core})
        except Exception:
            out.append({"name": pip_name, "import": imp, "ok": False,
                        "version": None, "core": core})
    return out


def check_tools() -> List[Dict[str, Any]]:
    out = []
    for name, purpose, required in EXTERNAL_TOOLS:
        path = shutil.which(name)
        out.append({"name": name, "purpose": purpose, "ok": path is not None,
                    "path": path, "required": required})
    return out


def check_r_packages() -> Dict[str, Any]:
    if shutil.which("Rscript") is None:
        return {"rscript": False, "packages": {p: False for p in R_PACKAGES},
                "note": "Rscript ausente (roda no conteiner de analise)."}
    expr = (
        "cat(paste(sapply(c(%s), function(p) "
        "paste0(p, '=', requireNamespace(p, quietly=TRUE))), collapse=';'))"
        % ",".join("'%s'" % p for p in R_PACKAGES)
    )
    try:
        r = subprocess.run(["Rscript", "-e", expr], capture_output=True, text=True, timeout=60)
        pkgs = {}
        for tok in r.stdout.strip().split(";"):
            if "=" in tok:
                k, v = tok.split("=", 1)
                pkgs[k] = v.strip().upper() == "TRUE"
        return {"rscript": True, "packages": pkgs or {p: False for p in R_PACKAGES}}
    except Exception as e:
        return {"rscript": True, "packages": {p: False for p in R_PACKAGES}, "error": str(e)}


def gpu_status() -> Dict[str, Any]:
    from bonder.gpu.detect import probe
    from bonder.gpu import gpu_gate
    try:
        dec = probe().to_dict()
    except Exception as e:
        dec = {"device": "cpu", "reason": "probe falhou: %s" % e}
    return {"device": dec, "gpu_free_mb": gpu_gate.gpu_free_mb(), "llm_present": gpu_gate.llm_present()}


def check_config() -> Dict[str, Any]:
    try:
        cfg = load_config()
        from bonder.db.baits import validate_baits
        warns = validate_baits(cfg.baits, cfg.config_dir)
        return {
            "ok": True,
            "taxid": cfg.taxid,
            "baits": cfg.bait_ids(),
            "contrasts": [c.get("name") for c in cfg.contrasts()],
            "n_samples": len(cfg.samples.get("samples", [])),
            "warnings": warns,
        }
    except Exception as e:
        return {"ok": False, "error": str(e)}


def check_databases() -> Dict[str, Any]:
    import json
    cfg = load_config()
    reg = cfg.data_root / "registry.json"
    entries = {}
    if reg.exists():
        try:
            entries = json.loads(reg.read_text(encoding="utf-8"))
        except Exception:
            entries = {}
    search_db = cfg.data_root / "search_db.fasta"
    return {
        "data_root": str(cfg.data_root),
        "offline": cfg.offline,
        "cached_sources": sorted(entries.keys()),
        "n_cached": len(entries),
        "search_db_fasta": search_db.exists(),
    }


def check_paths() -> List[Dict[str, Any]]:
    cfg = load_config()
    out = []
    for label, p in [("data_root", cfg.data_root), ("runs_dir", cfg.runs_dir), ("config_dir", cfg.config_dir)]:
        p = Path(p)
        writable = False
        try:
            p.mkdir(parents=True, exist_ok=True)
            probe = p / ".bonder_write_test"
            probe.write_text("ok", encoding="utf-8")
            probe.unlink(missing_ok=True)
            writable = True
        except Exception:
            writable = False
        out.append({"label": label, "path": str(p), "exists": p.exists(), "writable": writable})
    return out


def check_runs() -> Dict[str, Any]:
    from bonder.web import jobs
    cfg = load_config()
    runs = []
    if cfg.runs_dir.exists():
        for d in sorted(cfg.runs_dir.iterdir(), reverse=True):
            if d.is_dir():
                st = jobs.read_status(d.name) or {}
                runs.append({"run_id": d.name, "status": st.get("status", "unknown")})
    by_status: Dict[str, int] = {}
    for r in runs:
        by_status[r["status"]] = by_status.get(r["status"], 0) + 1
    return {"total": len(runs), "by_status": by_status, "recent": runs[:10]}


def full_report() -> Dict[str, Any]:
    """Auditoria completa da plataforma. overall='ok' se o nucleo esta sao."""
    deps = check_python_deps()
    tools = check_tools()
    cfg_report = check_config()
    paths = check_paths()

    core_missing = [d["name"] for d in deps if d["core"] and not d["ok"]]
    paths_ok = all(p["writable"] for p in paths)
    overall = "ok" if (not core_missing and cfg_report.get("ok") and paths_ok) else "error"

    return {
        "overall": overall,
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "system": {
            "bonder_version": __version__,
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "executable": sys.executable,
        },
        "core_missing": core_missing,
        "python_deps": deps,
        "tools": tools,
        "r_packages": check_r_packages(),
        "gpu": gpu_status(),
        "config": cfg_report,
        "databases": check_databases(),
        "paths": paths,
        "runs": check_runs(),
    }


# ---------------------------------------------------------------------------
# Autoteste do pipeline (sanidade de scripts e processos, sem GPU/R)
# ---------------------------------------------------------------------------
def _synthetic_maxquant(path: Path) -> None:
    import pandas as pd
    cols = {
        "Majority protein IDs": ["P16110", "NSs_ICOARACI", "P12345", "Q8R1", "CON__P00761"],
        "Gene names": ["Lgals3", "NSs", "Prey1", "Prey2", "Trypsin"],
        "Peptides": ["8", "12", "5", "6", "3"],
        "Razor + unique peptides": ["6", "10", "4", "5", "2"],
        "Reverse": ["", "", "", "", ""],
        "Potential contaminant": ["", "", "", "", "+"],
    }
    for s, base in [("NSs_IP_1", 5000), ("NSs_IP_2", 5200), ("NSs_IP_3", 4800)]:
        cols["LFQ intensity %s" % s] = [str(base + 100), str(base * 10), str(base), str(base // 2), "300"]
    for s in ["NSs_CTRL_1", "NSs_CTRL_2", "NSs_CTRL_3"]:
        cols["LFQ intensity %s" % s] = ["0", "50", "700", "0", "320"]
    pd.DataFrame(cols).to_csv(path, sep="\t", index=False)


def self_test() -> Dict[str, Any]:
    """Roda o pipeline em dados sinteticos, num diretorio temporario.

    Retorna passou/falhou por artefato + os passos do manifesto. Nao toca no
    diretorio de corridas real. Ferramentas externas ausentes contam como
    'pulado' (esperado), nao como falha.
    """
    import json
    from dataclasses import replace
    started = time.time()
    result: Dict[str, Any] = {"ok": False, "checks": [], "steps": []}
    try:
        cfg = load_config()
        with tempfile.TemporaryDirectory(prefix="bonder_selftest_") as td:
            td = Path(td)
            cfg_tmp = replace(cfg, runs_dir=td / "runs", data_root=td / "data")
            mq = td / "proteinGroups.txt"
            _synthetic_maxquant(mq)
            contrast = next((c for c in cfg_tmp.contrasts() if c.get("name") == "NSs_vs_CTRL"),
                            cfg_tmp.contrasts()[0])

            from bonder.pipeline import run_pipeline
            run_dir = Path(run_pipeline(mq, contrast, "selftest", cfg=cfg_tmp))

            checks = []
            for art in ("interactors.tsv", "manifest.json", "qc.json", "network.json", "results.json"):
                checks.append({"artifact": art, "ok": (run_dir / art).exists()})
            checks.append({"artifact": "saint/interaction.txt", "ok": (run_dir / "saint" / "interaction.txt").exists()})
            checks.append({"artifact": "dep/data.tsv", "ok": (run_dir / "dep" / "data.tsv").exists()})

            manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
            steps = manifest.get("steps", [])
            import pandas as pd
            table = pd.read_csv(run_dir / "interactors.tsv", sep="\t")
            checks.append({"artifact": "bait excluida dos preys",
                           "ok": "NSs_ICOARACI" not in set(table["prey"].astype(str))})
            checks.append({"artifact": "contaminante excluido",
                           "ok": "P00761" not in set(table["prey"].astype(str))})
            checks.append({"artifact": "prey real presente",
                           "ok": "P12345" in set(table["prey"].astype(str))})

            result["checks"] = checks
            result["steps"] = steps
            result["n_candidates"] = int(table.shape[0])
            result["ok"] = all(c["ok"] for c in checks)
    except Exception as e:
        import traceback
        result["error"] = str(e)
        result["traceback"] = traceback.format_exc()
    result["elapsed_s"] = round(time.time() - started, 2)
    return result
