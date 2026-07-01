#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gpu_gate.py - lease cooperativo da GPU (compartilhamento com o LLM local).

Mesma filosofia do gpu_gate.py da suite docking: antes de carregar um consumidor
de GPU (ESM2, DIA-NN, rescoring), checa via nvidia-smi se ha VRAM livre suficiente
e se nao ha um processo de LLM (Ollama/llama) ocupando a placa. Se nao houver
folga, espera com backoff ou deixa o chamador cair para CPU. Escreve
gpu_lease_status.json para o dashboard.
"""
from __future__ import annotations

import json
import subprocess
import time
from pathlib import Path
from typing import List, Optional, Tuple


def _smi(query: str) -> str:
    r = subprocess.run(
        ["nvidia-smi", query, "--format=csv,noheader,nounits"],
        capture_output=True, text=True, timeout=30,
    )
    return r.stdout.strip()


def gpu_free_mb() -> Optional[int]:
    """VRAM livre (MiB) ou None se nvidia-smi indisponivel."""
    try:
        out = _smi("--query-gpu=memory.free")
        return int(out.splitlines()[0].split(",")[0].strip())
    except Exception:
        return None


def llm_present(names: Tuple[str, ...] = ("ollama", "llama")) -> bool:
    """True se um processo de LLM estiver usando a GPU."""
    try:
        r = subprocess.run(
            ["nvidia-smi", "--query-compute-apps=name", "--format=csv,noheader"],
            capture_output=True, text=True, timeout=30,
        )
        apps = r.stdout.lower()
        return any(n in apps for n in names)
    except Exception:
        return False


def acquire(
    required_mb: int = 6000,
    llm_names: Tuple[str, ...] = ("ollama", "llama"),
    poll_seconds: int = 30,
    max_wait_minutes: int = 120,
    status_path: Optional[Path] = None,
) -> bool:
    """Tenta obter o lease da GPU. Retorna True se liberada, False se esgotar espera.

    Nao levanta: o chamador usa o retorno para decidir GPU vs CPU.
    """
    t0 = time.time()
    deadline = t0 + max_wait_minutes * 60
    while True:
        free = gpu_free_mb()
        llm = llm_present(llm_names)
        ok = (free is not None) and (free >= required_mb) and (not llm)
        _write_status(status_path, holder="bonder" if ok else None,
                      free_mb=free, llm=llm, required_mb=required_mb, granted=ok)
        if ok:
            return True
        if free is None:
            return False  # sem nvidia-smi: deixa o chamador cair para CPU
        if time.time() >= deadline:
            return False
        time.sleep(poll_seconds)


def release(status_path: Optional[Path] = None) -> None:
    _write_status(status_path, holder=None, free_mb=gpu_free_mb(), llm=llm_present(),
                  required_mb=0, granted=False)


def apply_torch_limits(memory_fraction: float = 0.6) -> None:
    """Limita a fracao de VRAM do processo torch e reduz fragmentacao."""
    import os
    os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")
    try:
        import torch
        if torch.cuda.is_available():
            torch.cuda.set_per_process_memory_fraction(float(memory_fraction), 0)
    except Exception:
        pass


def _write_status(path, **kw) -> None:
    if not path:
        return
    try:
        payload = {"updated": time.strftime("%Y-%m-%dT%H:%M:%S"), **kw}
        Path(path).write_text(json.dumps(payload, indent=2, ensure_ascii=True), encoding="utf-8")
    except Exception:
        pass
