#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""detect.py - probe robusto de GPU com fallback CPU garantido.

Nao basta torch.cuda.is_available(): fazemos um teste ponta-a-ponta (versao CUDA,
capability sm_86 da RTX 3060 e um matmul descartavel em try/except) para pegar
mismatch de driver/runtime que so estoura no primeiro kernel. Qualquer falha ->
device='cpu' com mensagem ASCII-safe. O resultado e cacheado e registrado no
manifesto da corrida (reprodutibilidade).
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Any, Dict, Optional, Tuple


@dataclass
class DeviceDecision:
    device: str            # 'cuda' ou 'cpu'
    dtype: str             # 'float16' ou 'float32'
    reason: str
    name: Optional[str] = None
    capability: Optional[Tuple[int, int]] = None
    total_memory_mb: Optional[int] = None

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


_CACHE: Optional[DeviceDecision] = None


def probe(min_capability: Tuple[int, int] = (8, 6), force_cpu: bool = False) -> DeviceDecision:
    """Decide o device de forma segura. Nunca levanta: degrada para CPU."""
    global _CACHE
    if _CACHE is not None:
        return _CACHE

    if force_cpu:
        _CACHE = DeviceDecision("cpu", "float32", "forcado por configuracao (prefer=cpu)")
        return _CACHE

    try:
        import torch  # noqa: F401
    except Exception as e:  # torch ausente
        _CACHE = DeviceDecision("cpu", "float32", "torch indisponivel: %s" % e)
        return _CACHE

    try:
        import torch
        if not torch.cuda.is_available() or torch.version.cuda is None:
            _CACHE = DeviceDecision("cpu", "float32", "CUDA indisponivel (sem GPU ou build CPU do torch)")
            return _CACHE
        cap = torch.cuda.get_device_capability(0)
        name = torch.cuda.get_device_name(0)
        total = int(torch.cuda.get_device_properties(0).total_memory / (1024 * 1024))
        if cap < tuple(min_capability):
            _CACHE = DeviceDecision(
                "cpu", "float32",
                "capability %s < minima %s" % (cap, tuple(min_capability)),
                name=name, capability=cap, total_memory_mb=total,
            )
            return _CACHE
        # teste real de kernel (pega mismatch driver/runtime)
        x = torch.randn(8, 8, device="cuda")
        _ = (x @ x).sum().item()
        _CACHE = DeviceDecision(
            "cuda", "float16",
            "GPU pronta: %s (cap %s, %d MiB)" % (name, cap, total),
            name=name, capability=cap, total_memory_mb=total,
        )
        return _CACHE
    except Exception as e:
        _CACHE = DeviceDecision("cpu", "float32", "falha ao inicializar CUDA (fallback CPU): %s" % e)
        return _CACHE


def reset_cache() -> None:
    global _CACHE
    _CACHE = None
