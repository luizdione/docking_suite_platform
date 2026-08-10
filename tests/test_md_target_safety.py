# -*- coding: utf-8 -*-
"""
Gate da Fase 4 -- seguranca de alvo / force field em md_simulation.py.

Prova que as 3 correcoes valem:
  #1 alvo desconhecido ABORTA (nao cai mais para o alvo default -- rodar outro
     sistema produziria trajetoria plausivel e cientificamente errada);
  #2 o alvo e seus metadados (pdb, n_fads, organism) vem do REGISTRY, nao mais de
     um TARGET_MAP fixo; so tem MD preparada quem declara structure.cofactor_copies
     -- e `n_fads` e o VALOR curado (hGR=1, nao derivado de n_protomers);
  #3 XML de FAD/ligante ausente ABORTA (nao reusa mais o XML de outro alvo sem
     declaracao explicita via --fad-xml/--lig-xml).

Nao requer OpenMM (opcional no ambiente de dev): injeta um modulo 'openmm' falso
antes de importar md_simulation. resolve_target_paths nao toca OpenMM; consulta o
registry real do repo.

Uso:
    pytest tests/test_md_target_safety.py -v
    python tests/test_md_target_safety.py
"""
import shutil
import sys
import tempfile
import types
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

# md_simulation faz `from openmm import ...` no topo e sai com codigo 2 se faltar.
# Um modulo falso deixa o import passar sem tocar nada de OpenMM (as funcoes que o
# usam nao sao exercitadas aqui).
try:
    import openmm  # noqa: F401
except Exception:
    _fake = types.ModuleType("openmm")
    for _n in ("app", "unit", "LangevinMiddleIntegrator", "VerletIntegrator",
               "AndersenThermostat", "MonteCarloBarostat", "Platform",
               "CustomExternalForce"):
        setattr(_fake, _n, object())
    sys.modules["openmm"] = _fake

import md_simulation  # noqa: E402


def _args(target, fad_xml=None, lig_xml=None, complex=None):
    return types.SimpleNamespace(target=target, complex=complex,
                                 fad_xml=fad_xml, lig_xml=lig_xml)


def _exit_code(fn, *a, **k):
    try:
        fn(*a, **k)
    except SystemExit as exc:
        return exc.code
    return None


# ------------------------------------------------------------------ #1 fallback

def test_unknown_target_aborts():
    """Alvo fora do registry aborta com exit != 0."""
    code = _exit_code(md_simulation.resolve_target_paths,
                      _args("ALVO_INEXISTENTE"), Path(tempfile.gettempdir()))
    assert code not in (None, 0), f"esperava abortar, exit={code}"


def test_unknown_target_not_silently_switched():
    """O alvo NAO pode ser trocado por outro (fallback silencioso morto)."""
    a = _args("ALVO_INEXISTENTE")
    _exit_code(md_simulation.resolve_target_paths, a, Path(tempfile.gettempdir()))
    assert a.target == "ALVO_INEXISTENTE", "alvo trocado em silencio (fallback vivo!)"


# ------------------------------------------------------------- #2 registry-driven

def test_registry_target_without_md_aborts():
    """Alvo que EXISTE no registry mas nao declara cofactor_copies aborta.

    TcCruzain esta no painel do SUB01 mas nao tem MD preparada -- nao pode rodar
    por default como se fosse um alvo com FAD."""
    code = _exit_code(md_simulation.resolve_target_paths, _args("TcCruzain"),
                      Path(tempfile.gettempdir()))
    assert code not in (None, 0), f"alvo sem MD preparada deveria abortar, exit={code}"


def test_md_ready_targets_are_the_three_flavoenzymes():
    """So TcTR/LiTR/hGR declaram cofactor_copies (MD preparada) no registry."""
    ready = set(md_simulation._md_ready_targets())
    assert ready == {"TcTR", "LiTR", "hGR"}, ready


def test_md_info_comes_from_registry():
    """pdb/n_fads/organism vem do registry (nao do TARGET_MAP, que foi removido)."""
    hgr = md_simulation._target_md_info("hGR")
    assert hgr["pdb"] == "3DJJ_clean.pdb", hgr
    assert hgr["n_fads"] == 1, "hGR e curado como 1 FAD (nao 2); ver cofactor_copies"
    assert "sapiens" in (hgr["organism"] or "").lower(), hgr
    tctr = md_simulation._target_md_info("TcTR")
    assert tctr["n_fads"] == 2 and tctr["pdb"] == "1BZL_clean.pdb", tctr


def test_target_map_is_gone():
    """O TARGET_MAP fixo foi removido do modulo (fonte unica = registry)."""
    assert not hasattr(md_simulation, "TARGET_MAP"), "TARGET_MAP ainda existe"


# ---------------------------------------------------------------- #3 reuso de XML

def test_missing_fad_xml_aborts_no_silent_reuse():
    """XML de FAD ausente aborta, sem reusar o de outro alvo (alvo com MD valida)."""
    tmp = Path(tempfile.mkdtemp())
    try:
        # LiTR tem MD preparada no registry, mas a base tmp nao tem o XML dele
        code = _exit_code(md_simulation.resolve_target_paths, _args("LiTR"), tmp)
        assert code not in (None, 0), f"XML ausente deveria abortar, exit={code}"
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_explicit_xml_is_respected():
    """Com --fad-xml/--lig-xml explicitos, resolve sem erro e sem reuso implicito."""
    tmp = Path(tempfile.mkdtemp())
    try:
        fad = tmp / "meu_fad.xml"; fad.write_text("<ForceField/>", encoding="utf-8")
        lig = tmp / "meu_lig.xml"; lig.write_text("<ForceField/>", encoding="utf-8")
        a = _args("hGR", fad_xml=fad, lig_xml=lig)
        target_info, sub, tagged = md_simulation.resolve_target_paths(a, tmp)
        assert a.fad_xml == fad and a.lig_xml == lig, "path explicito ignorado"
        assert tagged == "hGR" and target_info["n_fads"] == 1, target_info
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_known_target_resolves_to_own_files():
    """Alvo com MD e XMLs no lugar resolve para os arquivos do PROPRIO alvo."""
    tmp = Path(tempfile.mkdtemp())
    try:
        sub = tmp / "smoke_test_litr"
        sub.mkdir(parents=True)
        (sub / "cofactor_FAD_LiTR.xml").write_text("<ForceField/>", encoding="utf-8")
        (sub / "lig_LiTR_smirnoff.xml").write_text("<ForceField/>", encoding="utf-8")
        a = _args("LiTR")
        target_info, _sub, tagged = md_simulation.resolve_target_paths(a, tmp)
        assert a.fad_xml == sub / "cofactor_FAD_LiTR.xml", a.fad_xml
        assert a.lig_xml == sub / "lig_LiTR_smirnoff.xml", a.lig_xml
        assert "TcTR" not in str(a.fad_xml), "resolveu para arquivo de outro alvo"
    finally:
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
