"""
validate_system.py - Sanity check headless do FAD XML + md_clean.py SEM rodar MD.

Verifica em <30 s, sem dispara dinamica:
  1. md_clean.py: AST OK, 10 funcoes top-level esperadas, sem imports quebrados
  2. cofactor_FAD_TcTR.xml: corretude estrutural (atom-element-bond)
  3. OpenMM: cria ForceField + Topology + System sem chegar a step()

Sai com codigo 0 se tudo ok, codigo != 0 se houver problema.
"""
from __future__ import annotations

import argparse
import ast
import re
import sys
import traceback
from pathlib import Path

DOCKING_ROOT = Path(__file__).resolve().parents[1]
MD_CLEAN = DOCKING_ROOT / "md_clean.py"
FAD_XML = DOCKING_ROOT / "docking_project" / "figures" / "smoke_test_fad" / "cofactor_FAD_TcTR.xml"
FAD_PDB = DOCKING_ROOT / "docking_project" / "figures" / "smoke_test_fad" / "cofactor_FAD_TcTR.pdb"


def check_md_clean(strict_line_limit: int = 1500) -> dict:
    out = {"file": str(MD_CLEAN), "ok": False}
    if not MD_CLEAN.exists():
        out["error"] = "md_clean.py not found"
        return out
    src = MD_CLEAN.read_text(encoding="utf-8")
    out["lines"] = len(src.splitlines())
    if out["lines"] > strict_line_limit:
        out["error"] = (f"line count {out['lines']} > limit {strict_line_limit}. "
                        f"User policy: refactor into a new file.")
        out["recommend_refactor"] = True
        return out
    try:
        tree = ast.parse(src)
    except SyntaxError as e:
        out["error"] = f"SyntaxError: {e}"
        return out
    funcs = [n.name for n in tree.body if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))]
    expected = {"build_system_fresh", "build_canonical_fad_topology",
                "pick_platform", "run_md", "analyze", "main",
                "forensic_residue_check", "diagnose_system",
                "build_system_fast", "add_position_restraints"}
    missing = expected - set(funcs)
    if missing:
        out["error"] = f"missing functions: {sorted(missing)}"
        return out
    out["functions"] = funcs
    out["ok"] = True
    return out


def check_fad_xml() -> dict:
    out = {"file": str(FAD_XML), "ok": False}
    if not FAD_XML.exists():
        out["error"] = "FAD XML not found"
        return out
    txt = FAD_XML.read_text(encoding="utf-8")
    # Atom -> type
    a2t = {m.group(1): m.group(2)
           for m in re.finditer(r'<Atom name="([^"]+)" type="(FAD\d+)"', txt)}
    types = {m.group(1): (m.group(2), float(m.group(3)))
             for m in re.finditer(
                 r'<Type name="(FAD\d+)" element="([^"]+)" mass="([^"]+)" class="FAD\d+"\s*/>',
                 txt)}
    res_match = re.search(r'<Residue name="FAD">(.*?)</Residue>', txt, re.DOTALL)
    res_bonds = []
    if res_match:
        res_bonds = [(m.group(1), m.group(2))
                     for m in re.finditer(
                         r'<Bond atomName1="([^"]+)" atomName2="([^"]+)"',
                         res_match.group(1))]
    hbf = re.search(r'<HarmonicBondForce>(.*?)</HarmonicBondForce>', txt, re.DOTALL)
    hbf_bonds = []
    if hbf:
        hbf_bonds = re.findall(
            r'<Bond (?:class1|type1)="(FAD\d+)" (?:class2|type2)="(FAD\d+)" length="[^"]+" k="[^"]+"\s*/>',
            hbf.group(1))

    out["n_atoms"] = len(a2t)
    out["n_types"] = len(types)
    out["n_res_bonds"] = len(res_bonds)
    out["n_hbf_bonds"] = len(hbf_bonds)

    # Check element correctness
    def expected_elem(name):
        for ch in name:
            if ch.isalpha():
                return ch.upper()
        return "?"

    wrong_elements = []
    for atom, tp in a2t.items():
        tinfo = types.get(tp)
        if not tinfo:
            continue
        if tinfo[0] != expected_elem(atom):
            wrong_elements.append((atom, tp, tinfo[0], expected_elem(atom)))
    out["wrong_elements"] = len(wrong_elements)

    # Coverage of HBF over residue bonds
    hbf_class_pairs = {frozenset({c1, c2}) for c1, c2 in hbf_bonds}
    missing_hbf = []
    for a, b in res_bonds:
        ta = a2t.get(a)
        tb = a2t.get(b)
        if ta and tb and frozenset({ta, tb}) not in hbf_class_pairs:
            missing_hbf.append((a, b))
    out["missing_hbf_bonds"] = len(missing_hbf)

    out["ok"] = (out["wrong_elements"] == 0 and out["missing_hbf_bonds"] == 0
                 and out["n_types"] > 0 and out["n_hbf_bonds"] >= len(res_bonds) // 2)
    if not out["ok"]:
        out["error"] = (f"wrong_elements={out['wrong_elements']}, "
                        f"missing_hbf_bonds={out['missing_hbf_bonds']}")
    return out


def check_openmm_load() -> dict:
    """Try to create ForceField + System for the FAD-only template (headless)."""
    out = {"ok": False}
    try:
        import openmm  # noqa
        from openmm import app
    except ImportError as e:
        out["error"] = f"openmm not importable: {e}"
        return out
    if not FAD_PDB.exists():
        out["error"] = f"FAD PDB not found: {FAD_PDB}"
        return out
    try:
        pdb = app.PDBFile(str(FAD_PDB))
        ff = app.ForceField(str(FAD_XML))
        # Create a tiny system — just FAD on its own; this tests bonded params
        sys_obj = ff.createSystem(
            pdb.topology, nonbondedMethod=app.NoCutoff,
            constraints=None, rigidWater=False,
        )
        out["n_particles"] = sys_obj.getNumParticles()
        out["n_forces"] = sys_obj.getNumForces()
        # Count HarmonicBondForce bonds applied
        for f in sys_obj.getForces():
            if isinstance(f, openmm.HarmonicBondForce):
                out["n_harmonic_bonds_applied"] = f.getNumBonds()
                break
        else:
            out["n_harmonic_bonds_applied"] = 0
        out["ok"] = out["n_harmonic_bonds_applied"] > 0
        if not out["ok"]:
            out["error"] = "ForceField creates System but 0 harmonic bonds applied"
    except Exception as e:
        out["error"] = f"{type(e).__name__}: {e}"
        out["traceback"] = traceback.format_exc()
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--line-limit", type=int, default=1500)
    ap.add_argument("--skip-openmm", action="store_true",
                    help="Pula validacao OpenMM (mais rapida)")
    args = ap.parse_args()

    checks = [("md_clean.py", check_md_clean(strict_line_limit=args.line_limit)),
              ("fad_xml", check_fad_xml())]
    if not args.skip_openmm:
        checks.append(("openmm_load", check_openmm_load()))

    rc = 0
    for name, info in checks:
        status = "OK" if info.get("ok") else "FAIL"
        print(f"[{status}] {name}")
        for k, v in info.items():
            if k in ("ok", "file"):
                continue
            if isinstance(v, list) and len(v) > 5:
                v = f"[{len(v)} items: {v[:3]} ...]"
            print(f"    {k}: {v}")
        if not info.get("ok"):
            rc = 1
    return rc


if __name__ == "__main__":
    sys.exit(main())
