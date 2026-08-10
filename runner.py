# -*- coding: utf-8 -*-
"""
runner.py -- Orquestracao de execucoes e proveniencia (Fase 6).

Antes disto, a proveniencia de um resultado era "o Dione rodou o script certo":
nada no disco dizia qual versao do registry, qual binario do Vina ou quais gates
de qualidade produziram um numero. Este modulo torna a execucao um objeto de
primeira classe.

Conceitos
---------
RunSpec  O QUE rodar: campanha, pares, etapa (kind), overrides.
plan()   RunSpec resolvido contra o registry -> lista de Task (o "DAG"; hoje as
         tasks sao independentes, mas a estrutura ja permite dependencias).
Task     Uma unidade atomica: um par (alvo x ligante) numa etapa, com seus
         parametros, gates, entradas/saidas e -- apos rodar -- hashes, tempos e
         veredito dos gates.
execute()Materializa um Run em <data_root>/runs/<run_id>/run.json, contendo:
           - registry_snapshot embutido + sha256 (reprodutibilidade retroativa);
           - toolchain (versao + hash de cada executavel + GPU);
           - tasks[] com cmdline, hashes de entrada/saida, wall time e gates.
resume() Retoma um run: reexecuta so as tasks que ainda nao concluiram.
status() Le o run.json e resume o estado agregado.

Duas decisoes de projeto que valem o comentario
------------------------------------------------
1. Gate reprovado NAO e excecao. A task recebe status "gate_failed" e o run
   segue. So falha de infraestrutura (o executor explodiu) vira "error". Isto
   espelha docking_v2.run_docking(), que ja devolve (None, msg) em vez de levantar.
2. O executor e injetavel. O default para Vina chama docking_v2.run_docking() sem
   toca-la -- exatamente a fronteira com a Fase 4 (pacote engines/), que
   substituira o executor sem reescrever este arquivo. E o que permite testar a
   proveniencia sem depender de um binario de docking instalado.
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
from dataclasses import asdict, dataclass, field, fields
from datetime import datetime
from pathlib import Path

import paths
import registry

ROOT = Path(__file__).resolve().parent
SUBPROJETO_ROOT = ROOT / "subprojeto"

MANIFEST_NAME = "run.json"
RUN_SCHEMA = "docking-run/1"

# --------------------------------------------------------------------- gates
# Limiares de qualidade. Sao regras da plataforma, nao de uma campanha; um par
# pode sobrepor o de RMSD via validation.rmsd_pass_A no registry.
GATE_MIN_TORSIONS = 15       # PDBQT do ligante rotacionavel o bastante (plano F6/F7.2)
GATE_REDOCK_RMSD_A = 2.0     # redock <= 2 A e "reproduz a pose cristalografica"
GATE_SITE_PLDDT = 70.0       # confianca minima no sitio (engines estruturais)

# Etapa do pipeline_svg que cada kind alimenta. Default do SUB01; uma campanha
# pode declarar o proprio mapa em ui.pipeline_steps. NAO inventamos etapa nova:
# se o subprojeto nao a tiver, o status.json fica intocado (ver _update_subproject_status).
STEP_OF_KIND = {"docking": "docking_painel_9_alvos"}


# ===================================================================== modelos

@dataclass
class RunSpec:
    """O que se quer rodar. Serializavel; vai inteiro para o manifesto."""
    campaign_id: str
    pairs: list = field(default_factory=list)   # pair_ids; vazio = todos da campanha
    kind: str = "docking"
    engine: str | None = None                   # override; senao vem do protocolo/par
    exhaustiveness: int | None = None           # override do protocolo (Vina)
    seed: int | None = None
    notes: str = ""


@dataclass
class Task:
    """Uma unidade de execucao, com o veredito preenchido em runtime."""
    task_id: str
    kind: str
    engine: str
    campaign: str
    target: str
    ligand: str
    site: str | None = None
    protocol_id: str | None = None
    params: dict = field(default_factory=dict)
    gates: list = field(default_factory=list)
    inputs: dict = field(default_factory=dict)
    outputs: dict = field(default_factory=dict)
    # --- preenchido ao executar
    status: str = "planned"   # planned|running|done|gate_failed|error
    cmdline: list | None = None
    input_hashes: dict = field(default_factory=dict)
    output_hashes: dict = field(default_factory=dict)
    gate_results: list = field(default_factory=list)
    metrics: dict = field(default_factory=dict)
    wall_s: float | None = None
    started_at: str | None = None
    finished_at: str | None = None
    error: str | None = None


@dataclass
class Plan:
    spec: RunSpec
    tasks: list


@dataclass
class Run:
    run_id: str
    run_dir: Path
    manifest: dict
    plan: Plan


# ===================================================================== helpers

def _now_iso():
    return datetime.now().isoformat(timespec="seconds")


def _sha256_file(path):
    """sha256 hex de um arquivo, em blocos. None se ilegivel (nunca levanta)."""
    try:
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()
    except OSError:
        return None


def _cmp(value, op, threshold):
    return {
        ">=": value >= threshold,
        "<=": value <= threshold,
        ">": value > threshold,
        "<": value < threshold,
        "==": value == threshold,
    }[op]


def _rehydrate(cls, data):
    """Reconstroi uma dataclass a partir de um dict, ignorando chaves extras."""
    keep = {f.name for f in fields(cls)}
    return cls(**{k: v for k, v in data.items() if k in keep})


# ===================================================================== plano

def _resolve_pairs(spec):
    """Pares da campanha filtrados pelo spec (todos, se spec.pairs vazio)."""
    all_pairs = registry.list_pairs(spec.campaign_id)
    if not spec.pairs:
        return all_pairs
    want = set(spec.pairs)
    chosen = [p for p in all_pairs
              if registry.pair_id(p["campaign"], p["target"], p["ligand"]) in want]
    missing = want - {registry.pair_id(p["campaign"], p["target"], p["ligand"])
                      for p in chosen}
    if missing:
        raise registry.RegistryError(
            f"pares nao encontrados na campanha {spec.campaign_id}: {sorted(missing)}"
        )
    return chosen


def _build_gates(pair, protocol, engine):
    """Gates aplicaveis, derivados do par + protocolo. Declarativos: cada um
    nomeia a metrica que um avaliador sabe medir sobre os artefatos reais."""
    caps = protocol.get("capabilities") or []
    gates = []
    if engine == "vina" or "pose" in caps:
        gates.append({"name": "torsions_min", "metric": "ligand_torsions",
                      "op": ">=", "threshold": GATE_MIN_TORSIONS})
        gates.append({"name": "protonation", "metric": "pose_has_hydrogens",
                      "op": ">=", "threshold": 1})
    val = pair.get("validation") or {}
    if val.get("redock_reference_pdb"):
        gates.append({"name": "redock_rmsd", "metric": "redock_rmsd_A", "op": "<=",
                      "threshold": val.get("rmsd_pass_A", GATE_REDOCK_RMSD_A),
                      "reference_pdb": val["redock_reference_pdb"],
                      "rmsd_function": val.get("rmsd_function", "CalcRMS")})
    if "confidence" in caps or engine in ("boltz", "boltz2"):
        gates.append({"name": "site_plddt_min", "metric": "site_plddt",
                      "op": ">=", "threshold": GATE_SITE_PLDDT})
    return gates


def _task_params(spec, protocol):
    p = dict(protocol.get("params") or {})
    if spec.exhaustiveness is not None:
        p["exhaustiveness"] = spec.exhaustiveness
    if spec.seed is not None:
        p["seed"] = spec.seed
    return p


def plan(spec):
    """Resolve um RunSpec contra o registry num Plan de Tasks. Puro: nao executa,
    nao toca disco -- por isso e testavel sem Vina."""
    if isinstance(spec, dict):
        spec = _rehydrate(RunSpec, spec)
    legacy = registry.as_legacy_target_db(spec.campaign_id)
    tasks = []
    for pair in _resolve_pairs(spec):
        target, ligand = pair["target"], pair["ligand"]
        protocol = registry.resolve_protocol(spec.campaign_id, spec.kind, target, ligand)
        engine = (spec.engine
                  or protocol.get("engine")
                  or (pair.get("engines") or {}).get("primary_pose"))
        tgt = legacy.get(target, {})
        params = _task_params(spec, protocol)
        params["grid_size"] = list(tgt.get("grid_size") or [])
        params["grid_center"] = None            # resolvido pelo engine no runtime
        params["pdb"] = tgt.get("pdb")
        tasks.append(Task(
            task_id=registry.pair_id(spec.campaign_id, target, ligand),
            kind=spec.kind,
            engine=engine,
            campaign=spec.campaign_id,
            target=target,
            ligand=ligand,
            site=pair.get("site"),
            protocol_id=protocol.get("id"),
            params=params,
            gates=_build_gates(pair, protocol, engine),
        ))
    return Plan(spec=spec, tasks=tasks)


# ============================================================ medicao de gates

def _read_text(path):
    try:
        return Path(path).read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None


def _pdbqt_for(task):
    """Melhor PDBQT do ligante disponivel: a pose de saida, senao a entrada."""
    return (task.outputs.get("pose_pdbqt")
            or task.inputs.get("ligand_pdbqt"))


def _measure(metric, task, gspec):
    """(valor, nota). valor None => gate 'skipped' (nao da para medir agora)."""
    if metric == "ligand_torsions":
        text = _read_text(_pdbqt_for(task))
        if text is None:
            return None, "sem PDBQT do ligante para contar torsoes"
        n = sum(1 for ln in text.splitlines() if ln.startswith("BRANCH"))
        return n, f"{n} torsoes (BRANCH) no PDBQT"
    if metric == "pose_has_hydrogens":
        text = _read_text(_pdbqt_for(task))
        if text is None:
            return None, "sem pose para checar protonacao"
        has_h = 0
        for ln in text.splitlines():
            if ln.startswith(("ATOM", "HETATM")):
                toks = ln.split()
                if toks and toks[-1] in ("HD", "H"):
                    has_h = 1
                    break
        return has_h, ("pose com hidrogenios polares" if has_h
                       else "pose sem hidrogenios (HD/H)")
    if metric == "redock_rmsd_A":
        # A medicao real exige um engine de redock + RDKit CalcRMS; isso chega
        # com a Fase 4. Aqui o gate fica declarado e explicitamente adiado, em
        # vez de fingir um numero.
        return None, "redock RMSD requer engine de redock (Fase 4); gate declarado"
    if metric == "site_plddt":
        return None, "pLDDT so se aplica a engine estrutural (Boltz/AF)"
    return None, f"metrica desconhecida: {metric}"


def _run_gates(task):
    results = []
    for g in task.gates:
        try:
            value, note = _measure(g["metric"], task, g)
        except Exception as exc:  # avaliador nunca derruba o run
            results.append({**g, "status": "error", "value": None,
                            "note": f"avaliador falhou: {exc!r}"})
            continue
        if value is None:
            results.append({**g, "status": "skipped", "value": None, "note": note})
            continue
        ok = _cmp(value, g["op"], g["threshold"])
        results.append({**g, "status": "passed" if ok else "failed",
                        "value": value, "note": note})
    return results


# ================================================================== executores

def _engine_request(task, dry_run):
    """Monta um EngineRequest a partir de uma Task, resolvendo protocolo e
    target_info do registry. A pose ja tem caminho (o execute() setou), de onde
    sai o workdir do engine."""
    from engines import EngineRequest
    protocol = {"id": task.protocol_id}
    try:
        if task.protocol_id:
            protocol = registry.get_protocol(task.protocol_id)
    except Exception:
        pass
    target_info = {}
    try:
        target_info = registry.as_legacy_target_db(task.campaign).get(task.target, {})
    except Exception:
        pass
    pose = task.outputs.get("pose_pdbqt")
    workdir = Path(pose).parent if pose else None
    return EngineRequest(
        campaign=task.campaign, target=task.target, ligand=task.ligand,
        kind=task.kind, protocol=protocol, site=task.site,
        params=dict(task.params or {}), inputs=dict(task.inputs or {}),
        outputs=dict(task.outputs or {}), target_info=target_info,
        workdir=workdir, dry_run=dry_run)


def _engine_executor(task, dry_run=False):
    """Executor default: resolve o engine por nome (pacote engines/) e roda.

    Chama o contrato completo -- prepare() ANTES de run() (engines/base.py
    documenta esse ciclo desde a Fase 4, mas nada aqui o respeitava: ia direto
    para run(), que so consome entradas ja resolvidas e nunca as prepara). Sem
    isso, nenhuma campanha conseguia rodar docking real (dry_run=False): Vina
    recebia sempre inputs={} / grid_center=None e falhava alto (por design,
    nunca em silencio) em "receptor_pdbqt ausente".

    Converte EngineResult -> dict {cmdline, outputs, inputs, metrics, error?},
    superset da forma que _run_task ja consumia do antigo vina_executor -- por
    isso os testes de proveniencia, que injetam um executor fake, continuam
    validos sem mudanca (so ganham a chave 'inputs' extra, que ignoram).
    """
    try:
        from engines import get_engine
        engine = get_engine(task.engine)
    except Exception as exc:  # engine desconhecido/quebrado nao derruba o run
        return {"error": f"engine '{task.engine}' indisponivel: {exc!r}"}
    req = _engine_request(task, dry_run)
    try:
        req = engine.prepare(req) or req  # resolve entradas/saidas; nunca executa
    except Exception as exc:  # arquivo ausente, sitio invalido etc. -- erro claro
        return {"error": f"engine {task.engine} prepare() levantou: {exc!r}"}
    result = engine.run(req)
    out = {"cmdline": result.cmdline,
           "metrics": dict(result.metrics or {}),
           "outputs": dict(result.outputs or {}),
           "inputs": dict(req.inputs or {}),
           "params": dict(req.params or {})}
    if result.log:
        out["log_tail"] = result.log
    if result.status == "error" or result.error:
        out["error"] = result.error or f"engine {task.engine} status={result.status}"
    return out


def _executor_for(task, override):
    if override is not None:
        return override
    return _engine_executor


# =================================================================== toolchain

def toolchain_snapshot(planned=None):
    """Versao + hash de cada executavel + estado da GPU, via engines.toolchain
    (cascata PATH -> registry/toolchain.json -> heuristica conda)."""
    from engines import toolchain
    return toolchain.snapshot()


# ==================================================================== execucao

def _new_run_id(spec):
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    digest = hashlib.sha256(
        registry._canonical_json(asdict(spec)).encode("utf-8")).hexdigest()[:6]
    return f"{stamp}-{spec.kind}-{digest}"


def new_run_id(spec):
    """Id de um run, publico. Necessario para quem precisa do id ANTES de
    execute() terminar -- por exemplo o dashboard, que dispara execute() numa
    thread e ja devolve o id ao navegador para poder consultar o status. Chamar
    execute(run_id=new_run_id(spec)) reaproveita exatamente este id."""
    if isinstance(spec, dict):
        spec = _rehydrate(RunSpec, spec)
    return _new_run_id(spec)


def _manifest(spec, run_id, tasks, snap):
    return {
        "schema": RUN_SCHEMA,
        "run_id": run_id,
        "created_at": _now_iso(),
        "finished_at": None,
        "status": "running",
        "spec": asdict(spec),
        "platform": {"python": sys.version.split()[0], "os": sys.platform},
        "registry_sha256": registry.snapshot_sha256(snap),
        "toolchain": toolchain_snapshot(tasks),
        "registry_snapshot": snap,
        "tasks": [asdict(t) for t in tasks],
    }


def _write_manifest(run_dir, manifest):
    manifest["tasks_summary"] = _summarize([t["status"] for t in manifest["tasks"]])
    (run_dir / MANIFEST_NAME).write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def _summarize(statuses):
    out = {}
    for s in statuses:
        out[s] = out.get(s, 0) + 1
    return out


def _aggregate_status(tasks):
    st = {t.status for t in tasks}
    if not st or st <= {"done"}:
        return "completed"
    if st & {"error"}:
        return "error"
    if st & {"gate_failed"}:
        return "gate_failed"
    if st & {"planned", "running"}:
        return "incomplete"
    return "completed"


def _run_task(task, executor, dry_run):
    task.started_at = _now_iso()
    task.status = "running"
    task.input_hashes = {role: _sha256_file(pth)
                         for role, pth in task.inputs.items()
                         if pth and Path(pth).exists()}
    t0 = time.perf_counter()
    try:
        result = executor(task, dry_run=dry_run) or {}
    except Exception as exc:  # so infraestrutura chega aqui
        task.wall_s = round(time.perf_counter() - t0, 3)
        task.status = "error"
        task.error = f"executor levantou: {exc!r}"
        task.finished_at = _now_iso()
        return
    task.wall_s = round(time.perf_counter() - t0, 3)
    task.cmdline = result.get("cmdline")
    task.metrics = result.get("metrics", {}) or {}
    # params tambem e atualizado (nao so inputs/outputs): prepare() resolve
    # grid_center a partir de None -- sem isto, run.json fica com o valor de
    # planejamento (None), nao o que foi de fato usado no comando do Vina.
    for key, val in (result.get("params") or {}).items():
        task.params[key] = val
    for role, pth in (result.get("inputs") or {}).items():
        task.inputs[role] = str(pth)
    for role, pth in (result.get("outputs") or {}).items():
        task.outputs[role] = str(pth)
    # Recalcula input_hashes (nao so no topo desta funcao, onde task.inputs
    # normalmente ainda estava vazio): o executor pode ter preenchido entradas
    # via engine.prepare() nesse meio tempo -- sem isto, o run.json alegava
    # proveniencia sobre hash de arquivo que nunca foi de fato hasheado.
    task.input_hashes = {role: _sha256_file(pth)
                         for role, pth in task.inputs.items()
                         if pth and Path(pth).exists()}
    task.output_hashes = {role: _sha256_file(pth)
                          for role, pth in task.outputs.items()
                          if pth and Path(pth).exists()}
    if result.get("error"):
        task.status = "error"
        task.error = result["error"]
        task.finished_at = _now_iso()
        return
    task.gate_results = _run_gates(task)
    failed = [g for g in task.gate_results if g["status"] in ("failed", "error")]
    task.status = "gate_failed" if failed else "done"
    task.finished_at = _now_iso()


def execute(planned, executor=None, run_id=None, dry_run=False,
            write_subproject=True):
    """Executa um Plan, grava a proveniencia e (opcional) alimenta o SVG.

    executor: se dado, roda TODAS as tasks (injecao usada em teste e pela Fase 4);
              senao, escolhe por task.engine via ENGINE_EXECUTORS.
    run_id:   reaproveitado por resume(); senao gerado.
    """
    spec = planned.spec
    run_id = run_id or _new_run_id(spec)
    run_dir = paths.run_dir(spec.campaign_id, run_id, create=True)

    snap = registry.snapshot()
    manifest = _manifest(spec, run_id, planned.tasks, snap)
    _write_manifest(run_dir, manifest)  # grava antes de rodar: resume sobrevive a crash

    for task in planned.tasks:
        if task.status == "done":          # resume pula o que ja concluiu
            continue
        if not task.outputs.get("pose_pdbqt"):
            task.outputs["pose_pdbqt"] = str(
                run_dir / f"{task.target}__{task.ligand}.pose.pdbqt")
        _run_task(task, _executor_for(task, executor), dry_run)
        manifest["tasks"] = [asdict(t) for t in planned.tasks]
        _write_manifest(run_dir, manifest)  # incremental

    manifest["status"] = _aggregate_status(planned.tasks)
    manifest["finished_at"] = _now_iso()
    _write_manifest(run_dir, manifest)

    if write_subproject:
        info = _update_subproject_status(spec, planned.tasks)
        if info:
            manifest["subproject_status"] = info
            _write_manifest(run_dir, manifest)

    return Run(run_id=run_id, run_dir=run_dir, manifest=manifest, plan=planned)


def resume(campaign_id, run_id, executor=None, dry_run=False):
    """Retoma um run: reexecuta as tasks que nao concluiram, preservando o resto."""
    run_dir = paths.run_dir(campaign_id, run_id)
    manifest = json.loads((run_dir / MANIFEST_NAME).read_text(encoding="utf-8"))
    spec = _rehydrate(RunSpec, manifest["spec"])
    tasks = [_rehydrate(Task, d) for d in manifest["tasks"]]
    return execute(Plan(spec=spec, tasks=tasks), executor=executor,
                   run_id=run_id, dry_run=dry_run)


def status(campaign_id, run_id):
    """Estado agregado de um run, lido do run.json."""
    run_dir = paths.run_dir(campaign_id, run_id)
    manifest = json.loads((run_dir / MANIFEST_NAME).read_text(encoding="utf-8"))
    tasks = manifest["tasks"]
    gates = []
    for t in tasks:
        for g in t.get("gate_results", []):
            gates.append(g["status"])
    return {
        "run_id": run_id,
        "status": manifest["status"],
        "created_at": manifest.get("created_at"),
        "finished_at": manifest.get("finished_at"),
        "registry_sha256": manifest.get("registry_sha256"),
        "tasks": _summarize([t["status"] for t in tasks]),
        "gates": _summarize(gates),
        "path": str(run_dir / MANIFEST_NAME),
    }


# ====================================================== reuso: SVG do subprojeto

def _subproject_step(campaign_id, kind):
    try:
        ui = registry.get_campaign(campaign_id).get("ui") or {}
        steps = ui.get("pipeline_steps") or {}
        if kind in steps:
            return steps[kind]
    except Exception:
        pass
    return STEP_OF_KIND.get(kind)


def _etapa_status(tasks):
    st = {t.status for t in tasks}
    if st <= {"done"}:
        return "done"
    if st & {"gate_failed", "error"}:
        return "blocked"
    if "done" in st:
        return "data"          # parcial: ha resultado, mas nem tudo passou
    return "pending"


def _update_subproject_status(spec, tasks, root=None):
    """Atualiza SO a etapa executada em subprojeto/<ID>/status.json, preservando
    todas as demais (o SUB01 curou 11 etapas a mao, incluindo manuscrito). Nunca
    cria subprojeto nem etapa que nao existam -- o runner reflete, nao inventa."""
    root = Path(root) if root else SUBPROJETO_ROOT
    sp_id = spec.campaign_id
    step = _subproject_step(sp_id, spec.kind)
    if not step:
        return None
    st_path = root / sp_id / "status.json"
    if not st_path.exists():
        return None
    data = json.loads(st_path.read_text(encoding="utf-8"))
    etapas = data.get("etapas") or []
    novo = _etapa_status(tasks)
    for etapa in etapas:
        if etapa.get("etapa") == step:
            etapa["status"] = novo
            break
    else:
        return None   # etapa inexistente neste subprojeto: nao poluir
    data["etapas"] = etapas
    data.setdefault("_last_run", {})[step] = {
        "status": novo, "at": _now_iso(),
        "n_tasks": len(tasks),
        "by": {k: v for k, v in _summarize([t.status for t in tasks]).items()},
    }
    st_path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n",
                       encoding="utf-8")
    return {"path": str(st_path), "etapa": step, "status": novo}


# ========================================================================= CLI

def _cli(argv=None):
    import argparse

    ap = argparse.ArgumentParser(description="Runner de execucoes com proveniencia.")
    sub = ap.add_subparsers(dest="cmd", required=True)

    sp = sub.add_parser("plan", help="resolve e imprime o plano (nao executa)")
    sp.add_argument("--campaign", required=True)
    sp.add_argument("--kind", default="docking")
    sp.add_argument("--pairs", nargs="*", default=[])

    sr = sub.add_parser("run", help="executa um plano e grava a proveniencia")
    sr.add_argument("--campaign", required=True)
    sr.add_argument("--kind", default="docking")
    sr.add_argument("--pairs", nargs="*", default=[])
    sr.add_argument("--exhaustiveness", type=int, default=None)
    sr.add_argument("--dry-run", action="store_true")
    sr.add_argument("--no-subproject", action="store_true")

    st = sub.add_parser("status", help="resume um run existente")
    st.add_argument("--campaign", required=True)
    st.add_argument("--run", required=True)

    rs = sub.add_parser("resume", help="retoma um run existente")
    rs.add_argument("--campaign", required=True)
    rs.add_argument("--run", required=True)
    rs.add_argument("--dry-run", action="store_true")

    args = ap.parse_args(argv)

    if args.cmd == "plan":
        pl = plan(RunSpec(campaign_id=args.campaign, kind=args.kind, pairs=args.pairs))
        print(f"plano: {len(pl.tasks)} task(s) | campanha {args.campaign} | {args.kind}")
        for t in pl.tasks:
            gs = ",".join(g["name"] for g in t.gates) or "-"
            print(f"  {t.task_id:32s} engine={t.engine or '?':7s} "
                  f"proto={t.protocol_id or '?':14s} gates=[{gs}]")
        return 0

    if args.cmd == "run":
        pl = plan(RunSpec(campaign_id=args.campaign, kind=args.kind, pairs=args.pairs,
                          exhaustiveness=args.exhaustiveness))
        run = execute(pl, dry_run=args.dry_run,
                      write_subproject=not args.no_subproject)
        print(f"run {run.run_id} -> {run.manifest['status']}")
        print(f"  manifesto: {run.run_dir / MANIFEST_NAME}")
        print(f"  registry_sha256: {run.manifest['registry_sha256']}")
        for t in run.plan.tasks:
            print(f"  {t.task_id:32s} {t.status}")
        return 0 if run.manifest["status"] in ("completed",) else 1

    if args.cmd == "status":
        info = status(args.campaign, args.run)
        print(json.dumps(info, ensure_ascii=False, indent=2))
        return 0

    if args.cmd == "resume":
        run = resume(args.campaign, args.run, dry_run=args.dry_run)
        print(f"resume {run.run_id} -> {run.manifest['status']}")
        return 0 if run.manifest["status"] == "completed" else 1

    return 0


if __name__ == "__main__":
    raise SystemExit(_cli())
