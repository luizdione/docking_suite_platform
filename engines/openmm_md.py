# -*- coding: utf-8 -*-
"""
engines/openmm_md.py -- motor de dinamica molecular (OpenMM).

MD e cara (GPU, minutos a horas). Em vez de reimplementar o pipeline, este motor
INVOCA `md_simulation.py` -- que ja e um CLI completo e, desde a Fase 4 #1/#2/#3,
resolve o alvo pelo registry (resolve_target_paths) e falha alto em vez de cair
para TcTR. Assim a MD passa a ter uma entrada uniforme (Engine.run) sem duplicar
o codigo de build/simulacao.

dry_run monta o cmdline (proveniencia/plano) sem executar. Execucao real dispara
o script por subprocess.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from engines import toolchain
from engines.base import Engine, EngineRequest, EngineResult, register

ROOT = Path(__file__).resolve().parent.parent


def _production_ns(protocol, params):
    """ns de producao: params tem precedencia; senao, o maior 'ns' no schedule
    do protocolo de MD; senao None (o chamador decide o default)."""
    if params.get("ns") is not None:
        return params["ns"]
    ns = None
    for step in (protocol.get("schedule") or []):
        if isinstance(step, dict) and step.get("ns") is not None:
            ns = step["ns"] if ns is None else max(ns, step["ns"])
    return ns


@register
class OpenMMEngine(Engine):
    name = "openmm"

    def capabilities(self):
        return {"dynamics"}

    def check_toolchain(self):
        try:
            import openmm
            return {"available": True,
                    "tools": {"openmm": {"available": True,
                                         "version": getattr(openmm, "__version__", None)}},
                    "gpu": toolchain.gpu_entry()}
        except Exception as exc:  # openmm e opcional no env de dev
            return {"available": False,
                    "tools": {"openmm": {"available": False, "error": repr(exc)[:200]}},
                    "gpu": toolchain.gpu_entry()}

    def _cmdline(self, req):
        ns = _production_ns(req.protocol, req.params)
        script = ROOT / "md_simulation.py"
        cmd = [sys.executable, str(script), "--target", req.target]
        if ns is not None:
            cmd += ["--ns", str(ns)]
        tag = req.params.get("tag")
        if tag:
            cmd += ["--tag", str(tag)]
        mode = req.params.get("mode")
        if mode:
            cmd += ["--mode", str(mode)]
        if req.params.get("fresh"):
            cmd.append("--fresh")
        return cmd

    def run(self, req: EngineRequest) -> EngineResult:
        cmd = self._cmdline(req)
        if req.dry_run:
            return EngineResult(engine=self.name, status="ok", cmdline=cmd,
                                metrics={"mode": "dry_run",
                                         "protocol": req.protocol.get("id"),
                                         "ns": _production_ns(req.protocol, req.params)})
        tc = self.check_toolchain()
        if not tc.get("available"):
            return EngineResult(engine=self.name, status="error", cmdline=cmd,
                                error="OpenMM indisponivel neste ambiente "
                                      "(conda activate datascience/foldgpu com openmm+CUDA)")
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True,
                                  cwd=str(ROOT), timeout=req.params.get("timeout_s", 86400))
        except Exception as exc:
            return EngineResult(engine=self.name, status="error", cmdline=cmd,
                                error=f"md_simulation subprocess falhou: {exc!r}")
        if proc.returncode != 0:
            tail = (proc.stderr or proc.stdout or "")[-800:]
            return EngineResult(engine=self.name, status="error", cmdline=cmd,
                                error=f"md_simulation exit {proc.returncode}",
                                log=tail)
        return EngineResult(engine=self.name, status="ok", cmdline=cmd,
                            metrics={"returncode": 0},
                            log=(proc.stdout or "")[-2000:])
