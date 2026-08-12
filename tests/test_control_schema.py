# -*- coding: utf-8 -*-
"""Gate do vocabulario de CONTROLE (matriz de cross-docking de especificidade).

Um par de controle negativo (droga sem sitio no alvo, ex.: colchicina x actina) tem de
declarar isso EXPLICITAMENTE via bloco `control`, senao "afinidade fraca" fica
indistinguivel de "controle que funcionou". Este teste garante coerencia
(type<->expected_binding), a presenca do controle POSITIVO (colchicina x tubulina) e dos
NEGATIVOS de especificidade, e que o registry inteiro valida.

Uso:
    pytest tests/test_control_schema.py -v
    python tests/test_control_schema.py
"""
import sys
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import registry  # noqa: E402

PAIRS = ROOT / "registry" / "pairs"


def _pairs_with_control():
    out = {}
    for f in PAIRS.glob("*.json"):
        d = json.loads(f.read_text(encoding="utf-8"))
        if d.get("control"):
            out[f.stem] = d
    return out


def test_control_blocks_wellformed():
    """type valido, expected_binding booleano e COERENTE com o tipo, rationale presente."""
    pcs = _pairs_with_control()
    assert pcs, "nenhum par com bloco control encontrado"
    for name, d in pcs.items():
        c = d["control"]
        assert c.get("type") in ("positive", "negative"), (name, c.get("type"))
        assert isinstance(c.get("expected_binding"), bool), f"{name}: expected_binding nao booleano"
        assert (c["type"] == "positive") == c["expected_binding"], f"{name}: incoerente"
        assert c.get("rationale"), f"{name}: sem rationale"


def test_positive_control_present():
    """Braco POSITIVO da especificidade: colchicina x tubulina (deve ligar)."""
    pcs = _pairs_with_control()
    name = "SUB04__tubulin_ab__colchicina"
    assert name in pcs, f"{name} ausente"
    assert pcs[name]["control"]["type"] == "positive"
    assert pcs[name]["control"]["expected_binding"] is True


def test_negative_controls_present():
    """Negativos: colchicina x actina (mono e pentamero) e drogas de actina x tubulina."""
    negs = {n for n, d in _pairs_with_control().items() if d["control"]["type"] == "negative"}
    for expected in ("SUB04__TcAct1__colchicina",
                     "SUB04__TcAct1_factin_5mer__colchicina",
                     "SUB04__tubulin_ab__citocalasina-d",
                     "SUB04__tubulin_ab__jasplaquinolida"):
        assert expected in negs, f"{expected} ausente entre negativos {sorted(negs)}"


def test_registry_valid_with_controls():
    """O registry inteiro (com os pares de controle) valida sem erro."""
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
