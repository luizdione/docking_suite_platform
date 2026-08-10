# -*- coding: utf-8 -*-
"""
Gate da Fase 6 -- runner e proveniencia.

Prova o criterio de aceite ("um run de docking do SUB01 produz runs/<id>/run.json
com snapshot do registry e hashes; o SVG de /subprojects reflete esse run") sem
depender do Vina: injeta um executor falso deterministico que escreve um PDBQT
plausivel. O que se verifica aqui e a CAMADA DE PROVENIENCIA -- o motor de
docking real e a Fase 4.

Cobre:
  - run.json com registry_snapshot embutido + sha256 recomputavel;
  - toolchain presente;
  - tasks[] com cmdline, hashes de entrada/saida, wall time e gates;
  - gate reprovado -> status 'gate_failed' (NAO excecao);
  - erro de executor -> status 'error' capturado;
  - subprojeto/<ID>/status.json atualizado por MERGE (as 11 etapas curadas do
    SUB01 sobrevivem) e consumivel por pipeline_svg.status_of;
  - resume() pula tasks ja concluidas.

Uso:
    pytest tests/test_runner_provenance.py -v
    python tests/test_runner_provenance.py        # sem pytest instalado
"""
import importlib.util
import json
import shutil
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import paths           # noqa: E402
import registry        # noqa: E402
import runner          # noqa: E402

_ORIG_DATA_ROOT = paths.data_root
_ORIG_SUBPROJETO_ROOT = runner.SUBPROJETO_ROOT
_REAL_STATUS = ROOT / "subprojeto" / "SUB01" / "status.json"


# --------------------------------------------------------------- infraestrutura

def _fake_pdbqt(n_branch, with_h=True):
    """PDBQT minimo com n_branch torsoes (linhas BRANCH) e, opcional, um H polar."""
    out = ["REMARK  pose sintetica de teste", "ROOT",
           "ATOM      1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 C"]
    if with_h:
        out.append(
            "ATOM      2  H1  LIG A   1       1.000   0.000   0.000  1.00  0.00     0.000 HD")
    out.append("ENDROOT")
    for i in range(n_branch):
        a = i + 3
        out.append(f"BRANCH   1 {a:3d}")
        out.append(
            f"ATOM  {a:5d}  C{a} LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 C")
        out.append(f"ENDBRANCH   1 {a:3d}")
    out.append(f"TORSDOF {n_branch}")
    return "\n".join(out) + "\n"


def _make_executor(n_branch=20, with_h=True, fail=False):
    def _ex(task, dry_run=False):
        if fail:
            return {"error": "simulated engine failure"}
        dest = Path(task.outputs["pose_pdbqt"])
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(_fake_pdbqt(n_branch, with_h), encoding="utf-8")
        return {
            "cmdline": ["fake-vina", "--out", str(dest)],
            "outputs": {"pose_pdbqt": str(dest)},
            "metrics": {"affinities": [-9.5, -9.1, -8.7],
                        "best_affinity_kcal_mol": -9.5, "n_modes": 3},
        }
    return _ex


def _workspace():
    """Redireciona data_root e o subprojeto para um tmp isolado; devolve o tmp."""
    tmp = Path(tempfile.mkdtemp(prefix="runner_test_"))
    paths.data_root = lambda campaign_id=None, **kw: tmp / "data"
    sp = tmp / "subprojeto" / "SUB01"
    sp.mkdir(parents=True, exist_ok=True)
    shutil.copy(_REAL_STATUS, sp / "status.json")   # copia real: prova o merge
    runner.SUBPROJETO_ROOT = tmp / "subprojeto"
    return tmp


def _teardown(tmp):
    paths.data_root = _ORIG_DATA_ROOT
    runner.SUBPROJETO_ROOT = _ORIG_SUBPROJETO_ROOT
    try:
        shutil.rmtree(tmp, ignore_errors=True)
    except OSError:
        pass


def _load_manifest(run):
    return json.loads((run.run_dir / runner.MANIFEST_NAME).read_text(encoding="utf-8"))


# ---------------------------------------------------------------------- testes

def test_plan_resolves_sub01_pairs():
    """plan() e puro: resolve os pares do SUB01 em tasks de Vina com gates."""
    pl = runner.plan(runner.RunSpec(campaign_id="SUB01", kind="docking"))
    assert len(pl.tasks) >= 3, f"esperava >=3 pares, veio {len(pl.tasks)}"
    for t in pl.tasks:
        assert t.engine == "vina", f"{t.task_id}: engine {t.engine}"
        assert t.protocol_id == "vina_exh32", f"{t.task_id}: proto {t.protocol_id}"
        names = {g["name"] for g in t.gates}
        assert "torsions_min" in names and "protonation" in names, names


def test_run_produces_manifest_with_snapshot_and_hashes():
    tmp = _workspace()
    try:
        pl = runner.plan(runner.RunSpec(campaign_id="SUB01", kind="docking"))
        run = runner.execute(pl, executor=_make_executor(20), run_id="TEST-OK")
        m = _load_manifest(run)

        assert m["schema"] == runner.RUN_SCHEMA
        assert m["status"] == "completed", m["status"]
        # snapshot embutido + sha256 recomputavel e consistente
        assert m["registry_snapshot"], "registry_snapshot ausente"
        assert m["registry_sha256"] == registry.snapshot_sha256(), "sha256 nao bate"
        assert m["registry_snapshot"]["campaigns"].get("SUB01"), "snapshot sem SUB01"
        # toolchain
        assert m["toolchain"]["python"]["available"] is True
        assert "vina" in m["toolchain"] and "gpu" in m["toolchain"]
        # tasks
        for t in m["tasks"]:
            assert t["status"] == "done", (t["task_id"], t["status"], t.get("error"))
            assert t["wall_s"] is not None and t["wall_s"] >= 0
            assert t["cmdline"], "cmdline ausente"
            assert t["output_hashes"].get("pose_pdbqt"), "hash da pose ausente"
            g = {x["name"]: x["status"] for x in t["gate_results"]}
            assert g["torsions_min"] == "passed", g
            assert g["protonation"] == "passed", g
            if "redock_rmsd" in g:
                assert g["redock_rmsd"] == "skipped", g
    finally:
        _teardown(tmp)


def test_gate_failure_is_not_exception():
    """PDBQT com poucas torsoes: gate reprova, run segue, status gate_failed."""
    tmp = _workspace()
    try:
        pl = runner.plan(runner.RunSpec(campaign_id="SUB01", kind="docking"))
        run = runner.execute(pl, executor=_make_executor(3), run_id="TEST-GATE")
        m = _load_manifest(run)
        assert m["status"] == "gate_failed", m["status"]
        for t in m["tasks"]:
            assert t["status"] == "gate_failed", (t["task_id"], t["status"])
            g = {x["name"]: x["status"] for x in t["gate_results"]}
            assert g["torsions_min"] == "failed", g
            # a saida ainda foi hasheada (o artefato existe, so nao passou no gate)
            assert t["output_hashes"].get("pose_pdbqt")
    finally:
        _teardown(tmp)


def test_executor_error_is_captured():
    tmp = _workspace()
    try:
        pl = runner.plan(runner.RunSpec(campaign_id="SUB01", kind="docking"))
        run = runner.execute(pl, executor=_make_executor(fail=True), run_id="TEST-ERR")
        m = _load_manifest(run)
        assert m["status"] == "error", m["status"]
        for t in m["tasks"]:
            assert t["status"] == "error", (t["task_id"], t["status"])
            assert "simulated" in (t["error"] or ""), t["error"]
    finally:
        _teardown(tmp)


def test_subproject_status_merge_preserves_curated_steps():
    """O run atualiza SO a etapa de docking; as outras 10 curadas do SUB01 ficam."""
    tmp = _workspace()
    try:
        before = json.loads(_REAL_STATUS.read_text(encoding="utf-8"))
        n_before = len(before["etapas"])
        pl = runner.plan(runner.RunSpec(campaign_id="SUB01", kind="docking"))
        runner.execute(pl, executor=_make_executor(20), run_id="TEST-SP")

        after_path = tmp / "subprojeto" / "SUB01" / "status.json"
        data = json.loads(after_path.read_text(encoding="utf-8"))
        etapas = {e["etapa"]: e["status"] for e in data["etapas"]}

        assert len(data["etapas"]) == n_before, "numero de etapas mudou (merge destrutivo)"
        assert etapas["docking_painel_9_alvos"] == "done", etapas
        assert etapas["manuscrito_build_paper"] == "pending", "etapa curada perdida"
        assert data.get("metricas_chave"), "bloco curado metricas_chave sumiu"
        assert data.get("_last_run", {}).get("docking_painel_9_alvos"), "sem rastro do run"

        # e continua consumivel pelo pipeline_svg (a razao de existir do status.json)
        try:
            spec_ = importlib.util.spec_from_file_location(
                "pipeline_svg", ROOT / "subprojeto" / "pipeline_svg.py")
            mod = importlib.util.module_from_spec(spec_)
            spec_.loader.exec_module(mod)
            assert mod.status_of(data["etapas"], "docking_painel_9_alvos") == "done"
        except ImportError as exc:
            print(f"       (pipeline_svg nao importavel neste ambiente: {exc})")
    finally:
        _teardown(tmp)


def test_resume_skips_done_tasks():
    tmp = _workspace()
    try:
        pl = runner.plan(runner.RunSpec(campaign_id="SUB01", kind="docking"))
        run1 = runner.execute(pl, executor=_make_executor(20), run_id="TEST-RES")
        assert all(t.status == "done" for t in run1.plan.tasks)

        def _boom(task, dry_run=False):
            raise AssertionError(f"executor chamado para task ja concluida: {task.task_id}")

        run2 = runner.resume("SUB01", "TEST-RES", executor=_boom)
        m = _load_manifest(run2)
        assert m["status"] == "completed", m["status"]
        assert all(t["status"] == "done" for t in m["tasks"])
    finally:
        _teardown(tmp)


if __name__ == "__main__":
    failures = 0
    for _name, _fn in sorted(globals().items()):
        if not _name.startswith("test_") or not callable(_fn):
            continue
        try:
            _fn()
            print(f"[PASS] {_name}")
        except AssertionError as exc:
            failures += 1
            print(f"[FAIL] {_name}\n       {exc}")
        except Exception as exc:  # noqa: BLE001 -- bug no proprio teste/runner
            failures += 1
            print(f"[ERRO] {_name}\n       {type(exc).__name__}: {exc}")
    print(f"\n{failures} falha(s)")
    sys.exit(1 if failures else 0)
