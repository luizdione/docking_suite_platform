# -*- coding: utf-8 -*-
"""Gate de seguranca -- alvos OLIGOMERICOS (pentamero de F-actina) e o contra-alvo
tubulina declaram estado oligomerico e tem centro de caixa EXPLICITO (nunca caem no
centro de massa em silencio).

O sitio inter-protomerico da jasplaquinolida no pentamero e o sitio da colchicina na
tubulina usam a POSE NATIVA do ligante embutido (definition=explicit_center,
box_center_A gravado por site_derivation.validate_embedded_ligand). Este teste garante
que (a) os campos de assembly estao vivos e validados, (b) todo sitio oligomerico tem
box_center_A -> docking_v2 nao tem como cair no ultimo recurso (centro de massa).

Uso:
    pytest tests/test_oligomer_site_safety.py -v
    python tests/test_oligomer_site_safety.py
"""
import sys
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import registry  # noqa: E402

TARGETS = ROOT / "registry" / "targets"
PENTAMERS = ["TcAct1_factin_5mer", "TcAct2_factin_5mer", "TcAct3_factin_5mer",
             "TcAct4_factin_5mer", "Hs_ACTB_factin_5mer"]


def _t(tid):
    return json.loads((TARGETS / f"{tid}.json").read_text(encoding="utf-8"))


def test_pentamer_targets_declare_oligomer():
    """Os 5 pentameros declaram n_protomers_required=5 e assembly (campos antes mortos)."""
    for tid in PENTAMERS:
        s = _t(tid)["structure"]
        assert s.get("n_protomers_required") == 5, (tid, s.get("n_protomers_required"))
        assert isinstance(s.get("assembly"), str) and s["assembly"], f"{tid}: assembly vazio"


def test_pentamer_jasp_site_has_explicit_center():
    """Sitio da jasplaquinolida no pentamero: explicit_center com box_center_A [x,y,z]
    -- nunca ausente, logo docking nunca cai no centro de massa."""
    for tid in PENTAMERS:
        site = _t(tid)["sites"]["jasplaquinolida"]
        assert site["definition"] == "explicit_center", (tid, site["definition"])
        c = site.get("box_center_A")
        assert isinstance(c, list) and len(c) == 3, (tid, c)
        assert site.get("box_center_frame") == "receptor_native", tid


def test_tubulin_counter_target_colchicine_site():
    """Contra-alvo tubulina: categoria counter_target e sitio da colchicina explicito."""
    t = _t("tubulin_ab")
    assert t["category"] == "counter_target", t.get("category")
    s = t["structure"]
    assert (s.get("n_protomers_required") or 1) > 1 and s.get("assembly"), "assembly da tubulina"
    site = t["sites"]["colchicina"]
    assert site["definition"] == "explicit_center", site["definition"]
    assert isinstance(site.get("box_center_A"), list) and len(site["box_center_A"]) == 3


def test_every_multiprotomer_target_declares_assembly():
    """Invariante da validacao WS1: n_protomers_required>1 EXIGE assembly descrito."""
    for tid in registry._ids("targets"):
        s = (registry._read(registry._dir("targets") / f"{tid}.json").get("structure") or {})
        if (s.get("n_protomers_required") or 1) > 1:
            assert s.get("assembly"), f"{tid}: n_protomers>1 sem assembly (validacao WS1 falhou)"


def test_registry_valid_with_oligomers():
    """O registry inteiro (com os oligomeros + tubulina) valida sem erro."""
    errs = registry.validate_all()
    assert errs == [], "registry invalido:\n" + "\n".join(errs)


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
