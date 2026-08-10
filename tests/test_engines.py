# -*- coding: utf-8 -*-
"""
Gate da Fase 4 -- pacote engines/.

Verifica o contrato dos motores sem depender de Vina/Boltz/HADDOCK3 instalados:
  - os 4 engines se registram;
  - capabilities e dry_run (cmdline planejado) de cada um;
  - erro CLARO quando falta entrada (nunca resultado plausivel e errado);
  - geracao PURA de config do HADDOCK3 (.cfg + ambig.tbl) e do YAML do Boltz-2;
  - toolchain.find_tool (cascata) e o snapshot;
  - INTEGRACAO runner -> engines: execute(dry_run) usa o pacote de verdade
    (nao o executor fake dos testes de proveniencia).

Uso:
    pytest tests/test_engines.py -v
    python tests/test_engines.py
"""
import shutil
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import engines  # noqa: E402
from engines import EngineRequest, get_engine, list_engines, toolchain  # noqa: E402
from engines.haddock3 import build_ambig_tbl, build_cfg  # noqa: E402


def _req(**kw):
    base = dict(campaign="SUB01", target="TcTR", ligand="2-scn-nq", kind="docking")
    base.update(kw)
    return EngineRequest(**base)


# ------------------------------------------------------------------ registro

def test_all_engines_register():
    names = set(list_engines())
    assert {"vina", "boltz2", "openmm", "haddock3"} <= names, names


# ---------------------------------------------------------------------- vina

def test_vina_capabilities_and_dry_run():
    e = get_engine("vina")
    assert e.capabilities() == {"pose", "affinity"}
    r = e.run(_req(kind="docking",
                   params={"exhaustiveness": 32, "grid_center": [1, 2, 3],
                           "grid_size": [25, 25, 25]},
                   inputs={"receptor_pdbqt": "r.pdbqt", "ligand_pdbqt": "l.pdbqt"},
                   outputs={"pose_pdbqt": "out.pdbqt"}, dry_run=True))
    assert r.ok and r.metrics.get("mode") == "dry_run", r
    assert "--exhaustiveness" in r.cmdline


def test_vina_missing_input_is_clear_error():
    e = get_engine("vina")
    r = e.run(_req(kind="docking",
                   params={"grid_center": [1, 2, 3], "grid_size": [25, 25, 25]},
                   inputs={}, outputs={"pose_pdbqt": "out.pdbqt"}, dry_run=False))
    assert r.status == "error" and "receptor" in (r.error or ""), r


# ----------------------------------------------------------------- toolchain

def test_toolchain_find_tool_and_snapshot():
    assert toolchain.find_tool("tool_inexistente_xyz_123", conda_envs=()) is None
    snap = toolchain.snapshot()
    assert snap["python"]["available"] is True
    assert "vina" in snap and "gpu" in snap


# ----------------------------------------------------------------- haddock3

def test_haddock3_config_generation():
    cfg = build_cfg("run-X-Y", "rec.pdb", "lig.pdb",
                    {"modules": ["topoaa", "rigidbody", "caprieval"],
                     "sampling_rigidbody": 500})
    assert 'run_dir = "run-X-Y"' in cfg
    assert '"rec.pdb"' in cfg and '"lig.pdb"' in cfg
    assert "[rigidbody]" in cfg and "sampling = 500" in cfg
    assert 'ambig_fname = "ambig.tbl"' in cfg
    tbl = build_ambig_tbl([53, 58, 461])
    assert tbl.count("assign (") == 3
    assert "segid A" in tbl and "segid C" in tbl


def test_haddock3_dry_run_uses_registry_site_residues():
    e = get_engine("haddock3")
    r = e.run(_req(kind="pose", target="TcTR",
                   protocol={"id": "haddock3_default",
                             "params": {"modules": ["topoaa", "rigidbody", "caprieval"]}},
                   dry_run=True))
    assert r.ok, r
    assert r.metrics.get("n_restraints", 0) > 0, "TcTR tem fallback_residues no registry"


# ------------------------------------------------------------------- boltz2

def test_boltz2_missing_sequence_is_clear_error():
    e = get_engine("boltz2")
    r = e.run(_req(kind="structure",
                   protocol={"id": "boltz2_5samples", "params": {"diffusion_samples": 5}},
                   params={"smiles": "O=C1C=CC(=O)c2ccccc21"}, dry_run=True))
    assert r.status == "error" and "sequencia" in (r.error or ""), r


def test_boltz2_yaml_with_sequence():
    e = get_engine("boltz2")
    r = e.run(_req(kind="structure",
                   protocol={"id": "boltz2_5samples", "params": {"diffusion_samples": 5}},
                   inputs={"sequence": "MADEEKLRQ"},
                   params={"smiles": "O=C1C=CC(=O)c2ccccc21"}, dry_run=True))
    assert r.ok, r
    yaml = r.metrics.get("yaml", "")
    assert "MADEEKLRQ" in yaml and "O=C1C=CC(=O)c2ccccc21" in yaml, yaml


# ------------------------------------------------------------------- openmm

def test_openmm_dry_run_cmdline():
    e = get_engine("openmm")
    assert e.capabilities() == {"dynamics"}
    r = e.run(_req(kind="md", target="TcTR",
                   protocol={"id": "md_ff14sb_hmr_300k", "schedule": [{"ns": 100}]},
                   params={"tag": "TcTR_test"}, dry_run=True))
    assert r.ok, r
    assert any("md_simulation.py" in str(c) for c in r.cmdline)
    assert "--target" in r.cmdline and "TcTR" in r.cmdline
    assert r.metrics.get("ns") == 100


# --------------------------------------------------- integracao runner->engines

def test_runner_default_executor_uses_engines():
    """execute(dry_run) SEM executor injetado resolve o VinaEngine de verdade."""
    import paths
    import runner
    orig = paths.data_root
    tmp = Path(tempfile.mkdtemp(prefix="engines_test_"))
    paths.data_root = lambda campaign_id=None, **kw: tmp / "data"
    try:
        pl = runner.plan(runner.RunSpec(campaign_id="SUB01", kind="docking"))
        run = runner.execute(pl, dry_run=True, run_id="ENG-DRY",
                             write_subproject=False)
        assert run.manifest["status"] == "completed", run.manifest["status"]
        for t in run.plan.tasks:
            assert t.status == "done", (t.task_id, t.status, t.error)
            assert t.cmdline and t.cmdline[0] == "vina", t.cmdline
        # toolchain do manifesto agora vem de engines.toolchain
        assert run.manifest["toolchain"]["python"]["available"] is True
    finally:
        paths.data_root = orig
        shutil.rmtree(tmp, ignore_errors=True)


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
        except Exception as exc:  # noqa: BLE001
            failures += 1
            print(f"[ERRO] {_name}\n       {type(exc).__name__}: {exc}")
    print(f"\n{failures} falha(s)")
    sys.exit(1 if failures else 0)
