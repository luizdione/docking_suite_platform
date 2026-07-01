#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""tools.py - wrappers das ferramentas de GPU consagradas (revisadas por pares).

Todas rodam no conteiner GPU (Dockerfile.gpu) e passam pelo gpu_gate. Cada wrapper
levanta erro claro se o binario/lib faltar - nada e reimplementado.

  DIA-NN      : busca DIA acelerada por GPU (Demichev et al., Nat Methods 2020)
  MS2Rescore  : rescoring com features de deep learning (AlphaPeptDeep/MS2PIP/DeepLC)
  AlphaFold2-Multimer / ColabFold : validacao estrutural isca-presa dos candidatos
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import List, Optional


def _require(binary: str, purpose: str) -> str:
    exe = shutil.which(binary)
    if exe is None:
        raise FileNotFoundError(
            "%r nao encontrado (%s). Rode no conteiner GPU: 'docker compose run gpu ...'. "
            "Esta etapa e opcional e tem fallback." % (binary, purpose)
        )
    return exe


def run_diann(fasta: Path, raw_files: List[Path], out_dir: Path,
              binary: str = "diann", threads: int = 8, extra: Optional[List[str]] = None) -> Path:
    """Busca DIA com DIA-NN (GPU). Retorna o report.pg_matrix.tsv gerado."""
    exe = _require(binary, "busca DIA por GPU")
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [exe, "--fasta", str(fasta), "--out", str(out_dir / "report.tsv"),
           "--matrices", "--threads", str(threads)]
    for rf in raw_files:
        cmd += ["--f", str(rf)]
    if extra:
        cmd += list(extra)
    subprocess.run(cmd, check=True)
    pg = out_dir / "report.pg_matrix.tsv"
    if not pg.exists():
        raise RuntimeError("DIA-NN terminou sem gerar %s" % pg)
    return pg


def run_ms2rescore(psm_file: Path, spectrum_path: Path, out_dir: Path,
                   binary: str = "ms2rescore", extra: Optional[List[str]] = None) -> Path:
    """Rescoring com MS2Rescore (features AlphaPeptDeep/MS2PIP/DeepLC na GPU)."""
    exe = _require(binary, "rescoring por deep learning")
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [exe, "-p", str(psm_file), "-s", str(spectrum_path), "-o", str(out_dir)]
    if extra:
        cmd += list(extra)
    subprocess.run(cmd, check=True)
    return out_dir


def run_colabfold(fasta: Path, out_dir: Path,
                  binary: str = "colabfold_batch", extra: Optional[List[str]] = None) -> Path:
    """Predicao de complexo isca-presa (AF2-Multimer via ColabFold) para validacao."""
    exe = _require(binary, "validacao estrutural AF2-Multimer")
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [exe, str(fasta), str(out_dir)]
    if extra:
        cmd += list(extra)
    subprocess.run(cmd, check=True)
    return out_dir


def write_complex_fasta(bait_seq: str, prey_seq: str, bait_id: str, prey_id: str, out_path: Path) -> Path:
    """FASTA de complexo (2 cadeias) para AF2-Multimer/ColabFold."""
    out_path = Path(out_path); out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(">%s\n%s\n>%s\n%s\n" % (bait_id, bait_seq, prey_id, prey_seq), encoding="utf-8")
    return out_path
