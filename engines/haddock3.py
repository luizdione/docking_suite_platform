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


def build_ambig_tbl(residues, partner_residues=None, *, receptor_segid="A",
                    ligand_resid=1, ligand_segid="C"):
    """ambig.tbl (restricoes AIR no formato CNS) a partir dos residuos do sitio.

    Dois modos:
    - proteina-ligante pequeno (retrocompatibilidade, default):
      `partner_residues` None. Cada residuo ativo do receptor (`residues`,
      segid=receptor_segid) aponta para o UNICO resid do ligante
      (`ligand_resid`, segid=ligand_segid) -- o ligante pequeno e tratado como
      1 residuo so (comportamento original, inalterado).
    - proteina-proteina (2 cadeias, `partner_residues` fornecido -- lista de
      residuos ativos do PARCEIRO, ex.: TcFLA-1BP): AIR ativo-ativo
      BIDIRECIONAL -- cada residuo ativo do receptor aponta para o OR logico
      de TODOS os residuos de `partner_residues` (segid=ligand_segid), e cada
      residuo de `partner_residues` aponta de volta para o OR logico de todos
      os residuos ativos do receptor (segid=receptor_segid). E o formato AIR
      classico do HADDOCK (ativo-ativo, sem expansao por raio de passivos).
      Para restricoes mais ricas (active+passive por raio de contato, como no
      run real do fla_fla/run1-ecto), prefira o CLI oficial
      `haddock3-restraints active_passive_to_ambig` -- esta funcao e o minimo
      necessario para a PLATAFORMA gerar algo sozinha, sem depender do CLI
      externo.
    """
    lines = [
        "! AIRs HADDOCK3 geradas automaticamente a partir do registry (sites.*).",
        "! Numeracao = a do registry; pode exigir ajuste a numeracao do PDB.",
        "",
    ]
    if partner_residues:
        def _or_block(resids, segid):
            parts = [f"( resid {r} and segid {segid} )" for r in resids]
            return "\n     or ".join(parts)

        for rn in residues:
            lines += [
                f"assign ( resid {rn} and segid {receptor_segid} )",
                "       (",
                "        " + _or_block(partner_residues, ligand_segid),
                "       )  2.0 2.0 0.0",
                "",
            ]
        for rn in partner_residues:
            lines += [
                f"assign ( resid {rn} and segid {ligand_segid} )",
                "       (",
                "        " + _or_block(residues, receptor_segid),
                "       )  2.0 2.0 0.0",
                "",
            ]
    else:
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
              ambig_fname="ambig.tbl", protein_protein=False,
              caprieval_reference=None):
    """run.cfg do HADDOCK3 a partir do protocolo do registry.

    protein_protein: True pula `ligand_param_fname`/`ligand_top_fname` no
    modulo `topoaa` -- essas linhas so existem porque um ligante pequeno e um
    residuo NAO-padrao que precisa de topologia/parametro CNS explicitos; duas
    proteinas (cadeias de aminoacido padrao) tem a topologia resolvida sozinha
    pelo HADDOCK3, sem arquivo extra. Pode ser passado explicitamente ou ligado
    implicitamente por `protocol_params["protein_protein"] = true` no registry.
    Default False preserva o comportamento antigo (proteina-ligante).

    caprieval_reference: se fornecido (caminho pro PDB/CIF de referencia, ex.:
    o modelo do Boltz), adiciona um bloco [caprieval] FINAL (depois do ultimo
    modulo do protocolo) com `reference_fname` apontando pra ele. E o que faz
    o irmsd/lrmsd/dockq do capri_clt.tsv final virarem RMSD/DockQ VS. essa
    referencia, em vez de auto-referenciados (comparados contra o proprio
    melhor modelo do HADDOCK -- ver ressalva em haddock_analysis.parse_capri_clt
    e haddock_analysis.decide). Sem isso (default None), nenhum bloco extra e
    adicionado -- comportamento antigo preservado.
    """
    pp = protocol_params or {}
    protein_protein = bool(protein_protein) or bool(pp.get("protein_protein"))
    caprieval_reference = caprieval_reference or pp.get("caprieval_reference")
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
            if not protein_protein:
                out += ['ligand_param_fname = "lig.param"',
                        'ligand_top_fname = "lig.top"']
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

    if caprieval_reference:
        out += [
            "[caprieval]",
            f'reference_fname = "{caprieval_reference}"',
            "",
        ]
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
        params = req.params or {}
        protein_protein = bool(params.get("protein_protein"))
        caprieval_reference = params.get("caprieval_reference")
        partner_residues = _residue_numbers(params.get("partner_active_residues"))
        cfg = build_cfg(run_dir_name, receptor, ligand, pp,
                        protein_protein=protein_protein,
                        caprieval_reference=caprieval_reference)
        tbl = build_ambig_tbl(residues, partner_residues or None)
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
