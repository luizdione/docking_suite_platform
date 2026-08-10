# -*- coding: utf-8 -*-
"""
web_api.py -- camada de servico das rotas /api/* do dashboard (Fase 5.2 + T-0019).

Separa a LOGICA das rotas Flask: montar o bundle do registry para o front
(/api/config), resolver arquivos de estrutura para o visualizador 3Dmol.js
(/api/structure) e disparar/acompanhar o runner em background (/api/run*).

Fica fora do dashboard_v2.py de proposito -- assim e testavel sem subir o Flask
e o diff no dashboard (89 KB) se mantem minimo. As rotas em dashboard_v2.py sao
finas: parseiam a query, chamam uma funcao daqui e devolvem jsonify(...).

Nada aqui e especifico do SUB01: tudo vem do registry por campanha. As unicas
constantes de caminho sao as convencoes historicas de disco (docking_project/),
lidas em cascata junto com a raiz de dados da campanha (paths.py).
"""
from __future__ import annotations

import json
import threading
from dataclasses import asdict
from pathlib import Path

import paths
import registry
import runner

ROOT = Path(__file__).resolve().parent
DOCKING_PROJECT = ROOT / "docking_project"

# Extensao -> formato que o 3Dmol.js entende em addModel(data, format).
_FMT_BY_EXT = {
    ".pdb": "pdb",
    ".ent": "pdb",
    ".pdbqt": "pdbqt",
    ".cif": "cif",
    ".mmcif": "cif",
    ".sdf": "sdf",
    ".mol": "sdf",
    ".mol2": "mol2",
    ".xyz": "xyz",
}


# ===================================================================== helpers

def default_campaign():
    try:
        return registry.load_program().get("default_campaign", "SUB01")
    except Exception:
        return "SUB01"


def _first_existing(candidates):
    """Primeiro caminho existente da lista; se nenhum, tenta um glob
    case-insensitive no diretorio de cada candidato (portabilidade Linux, onde
    '1BZL.pdb' != '1bzl.pdb'). Retorna Path ou None."""
    cands = [Path(c) for c in candidates if c]
    for c in cands:
        if c.exists():
            return c
    for c in cands:
        parent = c.parent
        if not parent.exists():
            continue
        want = c.name.lower()
        for found in parent.iterdir():
            if found.name.lower() == want:
                return found
    return None


# =============================================================== /api/config

def _target_view(target, campaign_id):
    """Alvo achatado para o front. Superset do formato de dashboard._targets_json
    (mantido para a migracao incremental do JS nao quebrar), com extras usados
    pelo visualizador 3D e pelo runner."""
    struct = target.get("structure") or {}
    chem = target.get("chemistry") or {}
    tid = target["id"]
    pdb = struct.get("pdb_id")
    return {
        "id": tid,
        "name": target.get("name", ""),
        "target": target.get("name", ""),        # compat com o formato legado
        "organism": target.get("organism", ""),
        "category": target.get("category", "parasite"),
        "pdb": pdb or "N/A",
        "keep_cofactors": chem.get("keep_cofactors", []),
        "notes": target.get("notes", ""),
        "structure_source": struct.get("source", "pdb"),
        "ortholog_of": target.get("ortholog_of"),
        "default_site": target.get("default_site"),
        "sites": target.get("sites", {}),
        "has_receptor": _find_receptor(campaign_id, tid) is not None,
        "has_pose": _find_legacy_pose(campaign_id, tid) is not None,
    }


def _ligand_view(ligand):
    topo = ligand.get("topology") or {}
    return {
        "id": ligand["id"],
        "name": ligand.get("name", ""),
        "display_name": ligand.get("display_name", ""),
        "smiles": ligand.get("smiles", ""),
        "formula": ligand.get("formula", ""),
        "mw_g_mol": ligand.get("mw_g_mol"),
        "class": ligand.get("class", ""),
        "chemotype_tags": ligand.get("chemotype_tags", []),
        "is_macrocycle": topo.get("is_macrocycle", False),
    }


def _available_engines():
    try:
        import engines
        names = list(engines.list_engines())
    except Exception:
        # fallback: os quatro engines do pacote (Fase 4)
        names = ["vina", "boltz2", "openmm", "haddock3"]
    return sorted(names)


def _campaign_summaries():
    out = []
    for cid in registry.list_campaigns():
        try:
            c = registry.get_campaign(cid)
        except Exception:
            continue
        out.append({
            "id": cid,
            "sigla": c.get("sigla", ""),
            "titulo": c.get("titulo_paper_atual") or c.get("nome_especifico", ""),
            "status": c.get("status", ""),
        })
    return out


def build_config_bundle(campaign_id=None):
    """Bundle data-driven de uma campanha para o front (Fase 5.2).

    Um unico GET substitui os placeholders hard-coded do app.html
    (__TARGETS__, CURRENT_COMPOUND, TOX_DATA...). Levanta registry.RegistryError
    se a campanha nao existe -- a rota traduz em 404.
    """
    camp_id = campaign_id or default_campaign()
    campaign = registry.get_campaign(camp_id)   # valida a campanha (erra alto)

    targets = [_target_view(t, camp_id) for t in registry.list_targets(camp_id)]
    ligands = [_ligand_view(l) for l in registry.list_ligands(camp_id)]
    pairs = registry.list_pairs(camp_id)         # objetos crus: pequenos e uteis
    protocols = registry.list_protocols()

    return {
        "default_campaign": default_campaign(),
        "campaign_id": camp_id,
        "campaign": campaign,
        "campaigns": _campaign_summaries(),
        "ui": campaign.get("ui", {}),
        "default_ligand": campaign.get("default_ligand"),
        "targets": targets,
        "ligands": ligands,
        "pairs": pairs,
        "protocols": protocols,
        "engines": _available_engines(),
        "registry_sha256": _safe_registry_sha(),
    }


def _safe_registry_sha():
    try:
        return registry.snapshot_sha256()
    except Exception:
        return None


def bootstrap_for(campaign_id=None):
    """Bootstrap LEVE injetado no app.html (placeholder __BOOTSTRAP__), Fase 5.2/5.4.

    So o essencial para o cabecalho e o seletor de campanha serem data-driven -- o
    bundle completo (alvos/pares/protocolos) fica no /api/config, buscado sob demanda
    (aba 3D). Robusto: campanha invalida cai na default em vez de estourar o index()."""
    camp_id = campaign_id or default_campaign()
    try:
        campaign = registry.get_campaign(camp_id)
    except Exception:
        camp_id = default_campaign()
        try:
            campaign = registry.get_campaign(camp_id)
        except Exception:
            campaign = {}
    ui = campaign.get("ui") or {}
    # Resolve o ligante default (id -> nome/smiles) p/ o front nao precisar de fetch.
    default_lig_id = campaign.get("default_ligand")
    lig_name, lig_smiles = None, None
    if default_lig_id:
        try:
            lg = registry.get_ligand(default_lig_id)
            lig_name = lg.get("name") or lg.get("display_name")
            lig_smiles = lg.get("smiles")
        except Exception:
            pass
    return {
        "campaign_id": camp_id,
        "default_campaign": default_campaign(),
        "campaigns": _campaign_summaries(),
        "ui": ui,
        "recommended_targets": ui.get("recommended_targets", []),
        "primary_targets": campaign.get("primary_targets", []),
        "counter_targets": campaign.get("counter_targets", []),
        "default_ligand": default_lig_id,
        "default_ligand_name": lig_name,
        "default_ligand_smiles": lig_smiles,
        "sigla": campaign.get("sigla", ""),
        "titulo": campaign.get("titulo_paper_atual") or campaign.get("nome_especifico", ""),
    }


# ============================================================ /api/structure

def _find_receptor(campaign_id, target_id):
    """Arquivo de receptor (.pdb) de um alvo, preferindo a versao _clean.

    Procura na raiz de dados da campanha e, em cascata, na arvore historica
    docking_project/proteins/. O nome vem do PDB do alvo no registry.
    """
    try:
        target = registry.get_target(target_id)
    except Exception:
        return None
    pdb = ((target.get("structure") or {}).get("pdb_id") or "").strip()
    if not pdb:
        return None
    roots = []
    try:
        roots.append(paths.data_root(campaign_id) / "proteins")
    except Exception:
        pass
    roots.append(DOCKING_PROJECT / "proteins")
    cands = []
    for r in roots:
        cands += [r / f"{pdb}_clean.pdb", r / f"{pdb}.pdb"]
    return _first_existing(cands)


def _find_legacy_pose(campaign_id, target_id):
    """Pose de docking historica (docking_v2): docking/result_<ALVO>.pdbqt.

    O sufixo do arquivo casa com o id do alvo no registry (TcTR, hGR, ...), que
    preserva o case original -- por isso nao ha mapa de nomes aqui."""
    roots = []
    try:
        roots.append(paths.data_root(campaign_id) / "docking")
    except Exception:
        pass
    roots.append(DOCKING_PROJECT / "docking")
    cands = [r / f"result_{target_id}.pdbqt" for r in roots]
    return _first_existing(cands)


def _find_run_pose(campaign_id, run_id, target_id, ligand_id):
    """Pose gerada por um run do runner: <run_dir>/<alvo>__<ligante>.pose.pdbqt."""
    try:
        rd = paths.run_dir(campaign_id, run_id)
    except Exception:
        return None
    return _first_existing([rd / f"{target_id}__{ligand_id}.pose.pdbqt"])


def resolve_structure(campaign_id, what="receptor", target=None, run_id=None,
                      ligand=None):
    """Localiza um arquivo de estrutura e devolve seu conteudo + formato.

    what: 'receptor' | 'pose' | 'run_pose'
    Retorna dict {ok, format, data, filename, source, kind}. Levanta
    FileNotFoundError com mensagem clara se nao houver arquivo -- a rota
    traduz em 404 (o front mostra "estrutura ausente" em vez de tela preta).
    """
    camp = campaign_id or default_campaign()
    if not target:
        raise ValueError("parametro 'target' e obrigatorio")

    if what == "receptor":
        path = _find_receptor(camp, target)
        kind = "receptor"
    elif what == "pose":
        path = _find_legacy_pose(camp, target)
        kind = "pose_legacy"
    elif what == "run_pose":
        if not (run_id and ligand):
            raise ValueError("what='run_pose' exige 'run_id' e 'ligand'")
        path = _find_run_pose(camp, run_id, target, ligand)
        kind = "pose_run"
    else:
        raise ValueError(f"what invalido: {what!r} (use receptor|pose|run_pose)")

    if path is None:
        raise FileNotFoundError(
            f"estrutura ausente: what={what} campanha={camp} alvo={target}"
            + (f" run={run_id} ligante={ligand}" if what == "run_pose" else "")
        )

    fmt = _FMT_BY_EXT.get(path.suffix.lower(), "pdb")
    data = path.read_text(encoding="utf-8", errors="replace")
    return {
        "ok": True,
        "format": fmt,
        "data": data,
        "filename": path.name,
        "source": str(path),
        "kind": kind,
        "campaign": camp,
        "target": target,
    }


# =============================================================== /api/run*

# Registro em memoria dos jobs disparados neste processo. O estado autoritativo
# de progresso e o run.json (escrito incrementalmente por runner.execute); este
# dict guarda so o que o run.json ainda nao sabe: se a THREAD morreu antes mesmo
# de gravar o manifesto, e se e dry_run.
_JOBS: dict = {}
_JOBS_LOCK = threading.Lock()


def _write_starting_manifest(campaign, run_id, planned):
    """Grava um run.json 'starting' com as tasks planejadas ANTES de disparar a
    thread. Sem isso, o front veria 'starting' opaco por ~10 s -- o tempo que
    execute() gasta no primeiro toolchain snapshot antes de escrever o manifesto.
    execute() sobrescreve este placeholder com o manifesto real (toolchain,
    registry_snapshot, progresso). Nunca sobrescreve um run.json ja existente."""
    rd = paths.run_dir(campaign, run_id, create=True)
    mf = rd / runner.MANIFEST_NAME
    if mf.exists():
        return
    tasks = [asdict(t) for t in planned.tasks]
    manifest = {
        "schema": runner.RUN_SCHEMA,
        "run_id": run_id,
        "created_at": runner._now_iso(),
        "finished_at": None,
        "status": "starting",
        "spec": asdict(planned.spec),
        "tasks": tasks,
        "tasks_summary": {"planned": len(tasks)},
        "_placeholder": True,
    }
    mf.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
                  encoding="utf-8")


def start_run(campaign, kind="docking", pairs=None, exhaustiveness=None,
              dry_run=True, seed=None):
    """Planeja e dispara um run em background; devolve o run_id de imediato.

    dry_run=True (default seguro): os engines planejam o cmdline sem invocar
    binarios (vina.exe fora do PATH nesta maquina) -- roda em ms e nao escreve
    poses. Desmarque para execucao real.

    plan() roda SINCRONO aqui de proposito: se a campanha/par for invalida, o
    erro volta na hora (HTTP 400) em vez de sumir numa thread.
    """
    spec = runner.RunSpec(
        campaign_id=campaign,
        pairs=list(pairs or []),
        kind=kind,
        exhaustiveness=exhaustiveness,
        seed=seed,
    )
    planned = runner.plan(spec)                 # valida contra o registry (erra alto)
    run_id = runner.new_run_id(spec)

    with _JOBS_LOCK:
        _JOBS[run_id] = {
            "campaign": campaign,
            "kind": kind,
            "state": "starting",
            "error": None,
            "n_tasks": len(planned.tasks),
            "dry_run": bool(dry_run),
        }

    # Feedback imediato: grava o plano como run.json 'starting' ANTES de disparar
    # a thread. Sem isso o front veria ~10 s de 'starting' opaco -- o tempo do 1o
    # toolchain snapshot dentro de execute() antes do primeiro manifesto. Falha
    # aqui nao aborta o run: execute() regrava o manifesto real de qualquer forma.
    try:
        _write_starting_manifest(campaign, run_id, planned)
    except Exception:
        pass

    def _worker():
        try:
            runner.execute(planned, run_id=run_id, dry_run=dry_run)
            with _JOBS_LOCK:
                _JOBS[run_id]["state"] = "finished"
        except Exception as exc:  # infraestrutura: a task ja captura erro de engine
            with _JOBS_LOCK:
                _JOBS[run_id]["state"] = "crashed"
                _JOBS[run_id]["error"] = repr(exc)

    threading.Thread(target=_worker, name=f"run-{run_id}", daemon=True).start()

    return {
        "run_id": run_id,
        "campaign": campaign,
        "kind": kind,
        "n_tasks": len(planned.tasks),
        "dry_run": bool(dry_run),
        "state": "starting",
        "tasks": [t.task_id for t in planned.tasks],
    }


def _read_manifest(campaign, run_id):
    mf = paths.run_dir(campaign, run_id) / runner.MANIFEST_NAME
    return json.loads(mf.read_text(encoding="utf-8"))


def _task_detail(manifest):
    """Tasks resumidas para a UI de progresso: uma linha por par, com gates."""
    out = []
    for t in manifest.get("tasks", []):
        gates = [{"name": g.get("name"), "status": g.get("status"),
                  "value": g.get("value"), "note": g.get("note")}
                 for g in t.get("gate_results", [])]
        out.append({
            "task_id": t.get("task_id"),
            "target": t.get("target"),
            "ligand": t.get("ligand"),
            "engine": t.get("engine"),
            "status": t.get("status"),
            "wall_s": t.get("wall_s"),
            "error": t.get("error"),
            "metrics": t.get("metrics", {}),
            "gates": gates,
        })
    return out


def run_status(campaign, run_id):
    """Estado agregado + detalhe por task de um run, para o polling do front."""
    with _JOBS_LOCK:
        job = dict(_JOBS.get(run_id, {}))

    try:
        manifest = _read_manifest(campaign, run_id)
    except FileNotFoundError:
        # thread ainda nao gravou o primeiro manifesto (janela de ms), ou crashou antes
        return {
            "run_id": run_id,
            "status": job.get("state", "starting") if job else "unknown",
            "tasks": {},
            "gates": {},
            "job": job,
            "tasks_detail": [],
        }
    except Exception as exc:
        return {"run_id": run_id, "status": "error", "error": repr(exc), "job": job}

    agg = runner.status(campaign, run_id)   # resumo canonico, lido do mesmo run.json
    agg["job"] = job
    agg["dry_run"] = job.get("dry_run")
    agg["tasks_detail"] = _task_detail(manifest)
    agg["created_at"] = manifest.get("created_at")
    agg["finished_at"] = manifest.get("finished_at")
    agg["n_tasks"] = len(manifest.get("tasks", []))
    agg["n_done"] = sum(1 for t in manifest.get("tasks", [])
                        if t.get("status") in ("done", "gate_failed", "error"))
    return agg


def list_recent_runs(campaign, limit=25):
    """Runs recentes de uma campanha, do mais novo ao mais velho."""
    try:
        runs_root = paths.data_root(campaign) / "runs"
    except Exception:
        return []
    if not runs_root.exists():
        return []
    out = []
    for d in sorted(runs_root.iterdir(), key=lambda p: p.name, reverse=True):
        mf = d / runner.MANIFEST_NAME
        if not mf.exists():
            continue
        try:
            m = json.loads(mf.read_text(encoding="utf-8"))
        except Exception:
            continue
        out.append({
            "run_id": m.get("run_id", d.name),
            "status": m.get("status"),
            "kind": (m.get("spec") or {}).get("kind"),
            "created_at": m.get("created_at"),
            "finished_at": m.get("finished_at"),
            "tasks_summary": m.get("tasks_summary", {}),
        })
        if len(out) >= limit:
            break
    return out
