# -*- coding: utf-8 -*-
"""
engines/haddock3.py -- motor HADDOCK3 (docking guiado por restricoes/AIRs).

O que a plataforma tinha era um `run-TcTR-2SCN-NQ.cfg` e um `ambig.tbl` escritos
a mao para UM par. Aqui esses arquivos passam a ser GERADOS: o `.cfg` a partir de
protocols/haddock3_default.json (modulos, sampling) e o `ambig.tbl` a partir dos
residuos do sitio do alvo no registry (resolve_site). Isso e o que torna o
HADDOCK3 utilizavel para qualquer par, nao so para o SUB01.

A geracao (prepare) e logica pura e testavel, sem o binario. run() real dispara
`haddock3 <cfg>`; se o haddock3 nao estiver instalado, erro claro.

Ressalva honesta: a numeracao de residuo do registry pode diferir da do PDB
(ex.: no 1BZL as Cys catalíticas aparecem com offset +1). O header do .tbl avisa
disso -- ajuste ao PDB fica a cargo de quem roda, nao e adivinhado aqui.
"""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

from engines import toolchain
from engines.base import Engine, EngineRequest, EngineResult, register


def _residue_numbers(residues):
    """['CYS52','HIS461'] -> [52, 461]."""
    out = []
    for r in residues or []:
        m = re.search(r"(\d+)", str(r))
        if m:
            out.append(int(m.group(1)))
    return out


def _site_residues(req):
    """Residuos do sitio do alvo, do registry (fallback_residues do site)."""
    explicit = (req.params or {}).get("active_residues")
    if explicit:
        return _residue_numbers(explicit)
    try:
        import registry
        site = registry.resolve_site(req.target, req.site)
        return _residue_numbers(site.get("fallback_residues"))
    except Exception:
        return []


def build_ambig_tbl(residues, *, receptor_segid="A", ligand_resid=1,
                    ligand_segid="C"):
    """ambig.tbl (restricoes AIR no formato CNS) a partir dos residuos do sitio."""
    lines = [
        "! AIRs HADDOCK3 geradas automaticamente a partir do registry (sites.*).",
        "! Numeracao = a do registry; pode exigir ajuste a numeracao do PDB.",
        "",
    ]
    for rn in residues:
        lines += [
            f"assign ( resid {rn} and segid {receptor_segid} )",
            "       (",
            f"        ( resid {ligand_resid} and segid {ligand_segid} )",
            "       )  2.0 2.0 0.0",
            "",
        ]
    return "\n".join(lines)


def build_cfg(run_dir_name, receptor_pdb, ligand_pdb, protocol_params, *,
              ambig_fname="ambig.tbl"):
    """run.cfg do HADDOCK3 a partir do protocolo do registry."""
    pp = protocol_params or {}
    modules = pp.get("modules") or [
        "topoaa", "rigidbody", "caprieval", "flexref", "emref", "caprieval"]
    sampling = pp.get("sampling_rigidbody", 1000)
    select_rb = pp.get("select_rigidbody", 200)
    select_fr = pp.get("select_flexref", 50)

    out = [
        "# HADDOCK3 config GERADA de registry/protocols/haddock3_default.json",
        f'run_dir = "{run_dir_name}"',
        'mode = "local"',
        "ncores = 8",
        "postprocess = true",
        "clean = true",
        "",
        "molecules = [",
        f'    "{receptor_pdb}",',
        f'    "{ligand_pdb}",',
        "    ]",
        "",
    ]
    for m in modules:
        out.append(f"[{m}]")
        if m == "topoaa":
            out += ['ligand_param_fname = "lig.param"', 'ligand_top_fname = "lig.top"']
        elif m == "rigidbody":
            out += [f'ambig_fname = "{ambig_fname}"', f"sampling = {sampling}"]
        elif m in ("flexref", "emref"):
            out += [f'ambig_fname = "{ambig_fname}"']
        elif m == "seletop":
            out += [f"select = {select_rb}"]
        elif m == "seletopclusts":
            out += ["top_models = 4"]
        elif m == "flexref_select":
            out += [f"select = {select_fr}"]
        out.append("")
    return "\n".join(out)


@register
class Haddock3Engine(Engine):
    name = "haddock3"

    def capabilities(self):
        return {"pose"}

    def check_toolchain(self):
        hd = toolchain.resolve_known("haddock3")
        return {"available": bool(hd),
                "tools": {"haddock3": toolchain.tool_entry(hd) if hd
                          else {"available": False}}}

    def prepare(self, req: EngineRequest) -> EngineRequest:
        residues = _site_residues(req)
        pp = (req.protocol or {}).get("params") or {}
        receptor = req.inputs.get("receptor_pdb", f"{req.target}_receptor.pdb")
        ligand = req.inputs.get("ligand_pdb", f"ligand_{req.ligand}.pdb")
        run_dir_name = f"run-{req.target}-{req.ligand}"
        cfg = build_cfg(run_dir_name, receptor, ligand, pp)
        tbl = build_ambig_tbl(residues)
        req.outputs["_cfg_preview"] = cfg
        req.outputs["_tbl_preview"] = tbl
        req.outputs["_n_restraints"] = len(residues)
        if req.workdir and not req.dry_run:
            wd = Path(req.workdir)
            wd.mkdir(parents=True, exist_ok=True)
            cfg_path = wd / f"{run_dir_name}.cfg"
            tbl_path = wd / "ambig.tbl"
            cfg_path.write_text(cfg, encoding="utf-8")
            tbl_path.write_text(tbl, encoding="utf-8")
            req.outputs["cfg"] = str(cfg_path)
            req.outputs["ambig_tbl"] = str(tbl_path)
        return req

    def run(self, req: EngineRequest) -> EngineResult:
        req = self.prepare(req)
        if req.outputs.get("_n_restraints", 0) == 0:
            return EngineResult(engine=self.name, status="error",
                                error=f"nenhum residuo de sitio para o alvo "
                                      f"'{req.target}' (registry sites.* vazio); "
                                      "sem AIRs o HADDOCK3 nao tem o que guiar")
        cfg = req.outputs.get("cfg", f"run-{req.target}-{req.ligand}.cfg")
        cmd = ["haddock3", str(cfg)]
        if req.dry_run:
            return EngineResult(engine=self.name, status="ok", cmdline=cmd,
                                metrics={"mode": "dry_run",
                                         "n_restraints": req.outputs["_n_restraints"]})
        hd = toolchain.resolve_known("haddock3")
        if not hd:
            return EngineResult(engine=self.name, status="error", cmdline=cmd,
                                error="haddock3 nao encontrado (instale no env conda; "
                                      "a config .cfg/.tbl ja foi gerada no workdir)")
        cmd[0] = hd
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True,
                                  cwd=str(req.workdir) if req.workdir else None,
                                  timeout=req.params.get("timeout_s", 86400))
        except Exception as exc:
            return EngineResult(engine=self.name, status="error", cmdline=cmd,
                                error=f"haddock3 subprocess falhou: {exc!r}")
        if proc.returncode != 0:
            return EngineResult(engine=self.name, status="error", cmdline=cmd,
                                error=f"haddock3 exit {proc.returncode}",
                                log=(proc.stderr or proc.stdout or "")[-800:])
        return EngineResult(engine=self.name, status="ok", cmdline=cmd,
                            outputs={"cfg": cfg, "ambig_tbl": req.outputs.get("ambig_tbl")},
                            metrics={"returncode": 0}, log=(proc.stdout or "")[-2000:])
