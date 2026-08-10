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
        """Resolve receptor_pdbqt/ligand_pdbqt/grid_center antes de run().

        Ate T-0021 isto so ajustava o caminho de saida da pose -- nada aqui
        preparava entradas de verdade, e _engine_executor() nem chamava
        prepare() (ver runner.py). Resultado: nenhuma campanha conseguia
        rodar docking real (dry_run=False) via runner.execute(). Corrigido
        nos dois lados; este metodo fecha a metade que cabe ao Vina.
        """
        # A pose vai para o workdir do run (proveniencia). Se o runner ja setou o
        # caminho, respeita; senao deriva do workdir.
        if not req.outputs.get("pose_pdbqt") and req.workdir:
            req.outputs["pose_pdbqt"] = str(
                Path(req.workdir) / f"{req.target}__{req.ligand}.pose.pdbqt")

        if req.dry_run:
            return req  # run() ja retorna cedo em dry_run; nao ha o que preparar

        import registry
        import docking_v2

        # 1) target_info por PAR, nao por campanha inteira: campos gerais da
        # ponte legada + o SITIO DESTE PAR resolvido direto via
        # registry.resolve_site(target, site) -- NUNCA
        # as_legacy_target_db(campaign, site_id=req.site), que aplicaria esse
        # site_id a TODOS os alvos da campanha na mesma chamada e quebraria
        # qualquer alvo sem aquele sitio (T-0021, achado #6).
        try:
            legacy = registry.as_legacy_target_db(req.campaign).get(req.target, {})
        except Exception:
            legacy = {}
        target_info = dict(legacy)
        site = {}
        if req.site:
            try:
                site = registry.resolve_site(req.target, req.site)
            except Exception as exc:
                raise RuntimeError(
                    f"nao consegui resolver o sitio '{req.site}' de "
                    f"'{req.target}': {exc!r}") from exc
            if site.get("definition"):
                target_info["grid_strategy"] = site["definition"]
            if site.get("ligand_resnames"):
                target_info["ligand_resnames"] = site["ligand_resnames"]
            if site.get("fallback_residues"):
                target_info["grid_fallback_residues"] = site["fallback_residues"]
            if site.get("box_size_A"):
                target_info["grid_size"] = tuple(site["box_size_A"])
            if site.get("box_center_A"):
                target_info["box_center_A"] = tuple(site["box_center_A"])

        workdir = Path(req.workdir) if req.workdir else Path(".")
        workdir.mkdir(parents=True, exist_ok=True)

        # 2) receptor
        if not req.inputs.get("receptor_pdbqt"):
            prepared = site.get("prepared_receptor_local_file")
            if prepared and Path(prepared).exists():
                # Escape hatch DECLARADO no registry (dado, nao codigo): o
                # sitio aponta pra um receptor JA preparado externamente --
                # ex.: superposicao estrutural local que a plataforma ainda
                # nao automatiza (ver T-0021 Fase D). Usa direto; grid_center
                # vem de box_center_A, que so faz sentido fisico JUNTO com
                # este arquivo especifico (mesmo referencial).
                req.inputs["receptor_pdbqt"] = str(prepared)
                if not req.params.get("grid_center") and target_info.get("box_center_A"):
                    req.params["grid_center"] = tuple(target_info["box_center_A"])
            else:
                target = registry.get_target(req.target)
                struct = target.get("structure") or {}
                source = struct.get("source", "pdb")
                if source == "pdb":
                    pdb_path, msg = docking_v2.download_pdb(
                        target_info.get("pdb"), workdir)
                    if pdb_path is None:
                        raise RuntimeError(
                            f"download do PDB {target_info.get('pdb')} falhou: {msg}")
                elif source in ("alphafold", "local"):
                    local_file = struct.get("local_file")
                    if not local_file or not Path(local_file).exists():
                        raise RuntimeError(
                            f"structure.local_file ausente/inexistente p/ "
                            f"'{req.target}': {local_file!r}")
                    # Converte .cif->.pdb se preciso (Bio.PDB.PDBParser, usado
                    # por clean_protein, le .cif SEM erro mas com 0 atomos --
                    # falha silenciosa confirmada em T-0021 Fase 0.3).
                    pdb_path = docking_v2.ensure_pdb_format(
                        local_file, workdir, req.target)
                else:
                    raise RuntimeError(
                        f"structure.source desconhecido p/ '{req.target}': {source!r}")

                clean_path, pdbqt_path, grid_center, pdbqt_ok, pdbqt_msg = \
                    docking_v2.clean_protein(pdb_path, target_info, workdir)
                if not pdbqt_ok:
                    raise RuntimeError(
                        f"preparo do receptor falhou p/ '{req.target}': {pdbqt_msg}")
                req.inputs["receptor_pdbqt"] = str(pdbqt_path)
                if not req.params.get("grid_center"):
                    req.params["grid_center"] = grid_center

        # 3) ligante -- reusa a mesma cascata/cache do SUB01 (Fase B),
        # parametrizada pelo SMILES/nome do ligante deste par, nao mais
        # hardcoded no corpo de run_full_pipeline(). O cache fica na pasta de
        # ligantes DA CAMPANHA (paths.campaign_dirs, compartilhada entre runs
        # -- e onde docking_project/ligand/2-SCN-NQ/ ja vive p/ o SUB01), NAO
        # no workdir do run especifico (que e novo/vazio a cada execucao).
        if not req.inputs.get("ligand_pdbqt"):
            import paths
            ligand = registry.get_ligand(req.ligand)
            smiles = ligand.get("smiles")
            if not smiles:
                raise RuntimeError(f"ligante '{req.ligand}' sem SMILES no registry")
            ligand_dir = paths.campaign_dirs(req.campaign)["ligand"]
            pdbqt = docking_v2.prepare_ligand_pdbqt(
                smiles, ligand.get("name") or req.ligand, ligand_dir)
            if pdbqt is None:
                raise RuntimeError(
                    f"preparo do ligante '{req.ligand}' falhou (ver log 'pipeline')")
            req.inputs["ligand_pdbqt"] = str(pdbqt)

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
            # docking_v2.find_vina() (default de run_docking quando vina_path=None)
            # so procura em caminhos relativos ao CWD -- NAO consulta
            # registry/toolchain.json. check_toolchain() ja resolve certo
            # (find_vina() -> toolchain.resolve_known(), que consulta o
            # override); reusar aqui, senao o override cadastrado no
            # toolchain.json fica sem efeito na execucao real, so no relatorio
            # de status.
            vina_bin = docking_v2.find_vina() or toolchain.resolve_known("vina")
            affinities, log = docking_v2.run_docking(
                receptor, ligand, out, center, p["grid_size"],
                exhaustiveness=p.get("exhaustiveness", 32), vina_path=vina_bin)
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
