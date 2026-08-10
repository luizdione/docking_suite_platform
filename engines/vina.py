# -*- coding: utf-8 -*-
"""
engines/vina.py -- motor AutoDock Vina (docking).

Envolve `docking_v2.run_docking()` SEM modifica-la (a assinatura ja esta certa) --
era exatamente o ponto de extensao previsto no plano. Toda a logica que vivia
inline em runner.vina_executor (Fase 6) migra para ca, agora atras do contrato
Engine. O runner passa a resolver o motor por nome via engines.get_engine.
"""
from __future__ import annotations

from pathlib import Path

from engines import toolchain
from engines.base import Engine, EngineRequest, EngineResult, register


@register
class VinaEngine(Engine):
    name = "vina"

    def capabilities(self):
        return {"pose", "affinity"}

    def check_toolchain(self):
        # docking_v2.find_vina() ja sabe resolver o vina.exe de forma absoluta
        # (nao depende do CWD); toolchain.resolve_known e o fallback.
        vina = None
        try:
            import docking_v2
            vina = docking_v2.find_vina()
        except Exception:
            vina = None
        if not vina:
            vina = toolchain.resolve_known("vina")
        return {
            "available": bool(vina),
            "tools": {"vina": toolchain.tool_entry(vina) if vina
                      else {"available": False}},
        }

    def prepare(self, req: EngineRequest) -> EngineRequest:
        # A pose vai para o workdir do run (proveniencia). Se o runner ja setou o
        # caminho, respeita; senao deriva do workdir.
        if not req.outputs.get("pose_pdbqt") and req.workdir:
            req.outputs["pose_pdbqt"] = str(
                Path(req.workdir) / f"{req.target}__{req.ligand}.pose.pdbqt")
        return req

    def run(self, req: EngineRequest) -> EngineResult:
        p = req.params or {}
        receptor = req.inputs.get("receptor_pdbqt")
        ligand = req.inputs.get("ligand_pdbqt")
        out = req.outputs.get("pose_pdbqt")
        center = p.get("grid_center")
        cmd_preview = ["vina", "--receptor", str(receptor), "--ligand", str(ligand),
                       "--exhaustiveness", str(p.get("exhaustiveness", 32)),
                       "--out", str(out)]

        if req.dry_run:
            return EngineResult(engine=self.name, status="ok", cmdline=cmd_preview,
                                metrics={"mode": "dry_run"})
        if not receptor or not Path(receptor).exists():
            return EngineResult(engine=self.name, status="error", cmdline=cmd_preview,
                                error=f"receptor_pdbqt ausente: {receptor} "
                                      "(prepare o alvo antes; ver run_full_pipeline)")
        if not ligand or not Path(ligand).exists():
            return EngineResult(engine=self.name, status="error", cmdline=cmd_preview,
                                error=f"ligand_pdbqt ausente: {ligand}")
        if not center or not p.get("grid_size"):
            return EngineResult(engine=self.name, status="error", cmdline=cmd_preview,
                                error="grid_center/grid_size nao resolvidos para o alvo")
        try:
            import docking_v2
            affinities, log = docking_v2.run_docking(
                receptor, ligand, out, center, p["grid_size"],
                exhaustiveness=p.get("exhaustiveness", 32))
        except Exception as exc:  # falha de infraestrutura -> erro, nao excecao
            return EngineResult(engine=self.name, status="error", cmdline=cmd_preview,
                                error=f"run_docking levantou: {exc!r}")
        if affinities is None:
            return EngineResult(engine=self.name, status="error", cmdline=cmd_preview,
                                error=f"vina falhou: {str(log)[:300]}")
        return EngineResult(
            engine=self.name, status="ok", cmdline=cmd_preview,
            outputs={"pose_pdbqt": out},
            metrics={"affinities": affinities,
                     "best_affinity_kcal_mol": min(affinities),
                     "n_modes": len(affinities)},
            log=str(log)[-2000:])
