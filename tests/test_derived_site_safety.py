# -*- coding: utf-8 -*-
"""
Gate de seguranca -- sitios 'derived_from_reference_complex' (SUB04) nunca caem
em centro de massa da proteina em SILENCIO.

Antes desta correcao, docking_v2._calculate_grid_center() nao tinha ramo para a
estrategia 'derived_from_reference_complex' e, sem box_center_A, caia no ULTIMO
RECURSO (centro de massa da proteina inteira) -- um valor plausivel e
cientificamente errado, sem nenhum aviso. O grid ficaria no meio da proteina em
vez do bolso de ligacao derivado do complexo de referencia.

A correcao trata 'derived_from_reference_complex' com a MESMA regra fail-loud do
'explicit_center': le box_center_A do registry (derivado offline por
tools/derive_sites.py) ou LEVANTA erro claro -- NUNCA centro de massa.

Uso:
    pytest tests/test_derived_site_safety.py -v
    python tests/test_derived_site_safety.py
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import docking_v2  # noqa: E402


# --- estrutura fake: COM = (1,1,1), detectavel se o ramo cair no ultimo recurso -
class _Vec:
    def __init__(self, a): self.a = np.array(a, dtype=float)
    def get_array(self): return self.a


class _Atom:
    def __init__(self, xyz, name="CA", element="C"):
        self._xyz = xyz; self.name = name; self.element = element
    def get_vector(self): return _Vec(self._xyz)
    def get_coord(self): return np.array(self._xyz, dtype=float)


class _Res:
    def __init__(self, atoms, resname="ALA", resnum=1):
        self._atoms = atoms; self._resname = resname
        self.id = (" ", resnum, " ")
    def get_resname(self): return self._resname
    def get_atoms(self): return list(self._atoms)
    def __iter__(self): return iter(self._atoms)


class _Chain:
    def __init__(self, residues, cid="A"):
        self._res = residues; self.id = cid
    def __iter__(self): return iter(self._res)


class _Model:
    def __init__(self, chains): self._chains = chains
    def __iter__(self): return iter(self._chains)


class _Struct:
    """Todos os atomos em (1,1,1) -> centro de massa = (1,1,1)."""
    def __init__(self):
        atoms = [_Atom([1.0, 1.0, 1.0]) for _ in range(5)]
        self._models = [_Model([_Chain([_Res(atoms)])])]
    def __iter__(self): return iter(self._models)
    def get_atoms(self):
        for m in self._models:
            for c in m:
                for r in c:
                    for a in r:
                        yield a


def _raises(fn, *a, **k):
    try:
        out = fn(*a, **k)
    except Exception as exc:  # noqa: BLE001
        return exc
    return ("RETORNOU", out)


# ---------------------------------------------------------------- testes

def test_derived_without_center_raises_not_com():
    """derived_from_reference_complex SEM box_center_A: erro claro, NAO (1,1,1)."""
    r = _raises(docking_v2._calculate_grid_center, _Struct(),
                {"grid_strategy": "derived_from_reference_complex"})
    assert not (isinstance(r, tuple) and r and r[0] == "RETORNOU"), (
        f"caiu no centro de massa em silencio: {r}")
    assert isinstance(r, Exception), f"esperava excecao, veio {r!r}"


def test_derived_with_center_returns_it():
    """Com box_center_A gravado (derivado offline), retorna esse centro."""
    out = docking_v2._calculate_grid_center(
        _Struct(),
        {"grid_strategy": "derived_from_reference_complex",
         "box_center_A": [5.0, 6.0, 7.0]})
    assert tuple(round(float(x), 2) for x in out) == (5.0, 6.0, 7.0), out


def test_explicit_center_still_fail_loud():
    """Regressao: explicit_center sem box_center_A continua erro claro (nao COM)."""
    r = _raises(docking_v2._calculate_grid_center, _Struct(),
                {"grid_strategy": "explicit_center"})
    assert isinstance(r, Exception), f"explicit_center sem centro deveria erro, veio {r!r}"


def test_com_path_still_reachable_for_legacy_strategies():
    """Sanidade: uma estrategia legada sem match ainda pode usar o ultimo recurso
    (o guard fail-loud e SO para explicit/derived, nao quebrou o resto)."""
    out = docking_v2._calculate_grid_center(_Struct(), {"grid_strategy": "ligand_center"})
    assert tuple(round(float(x), 2) for x in out) == (1.0, 1.0, 1.0), out


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
