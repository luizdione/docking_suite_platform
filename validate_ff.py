#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
validate_ff.py - PORTAO DE VALIDACAO de force-field (roda em segundos).

Pega os bugs que custaram semanas ANTES de rodar a MD (que so falhava como NaN
cript­ico depois de minutos). Checa qualquer XML de FF OpenMM de residuo custom
(FAD, LIG, cofatores) e, opcionalmente, o PDB template:

  1. PeriodicTorsionForce NAO pode estar vazio        (bug: 0 torsions)
  2. HarmonicBondForce: bonds presentes, todas com k>0 (bug: HBF vazio/k=0)
  3. HarmonicAngleForce: angles presentes
  4. Cargas: soma ~ inteiro, sem |q| absurdo (>2 e)    (bug: cargas ruins)
  5. LJ sigma coerente com o ELEMENTO de cada Type     (bug: elemento/LJ trocado)
  6. (se nomes quimicos) elemento <-> 1a letra do nome  (bug: scramble de nomes)
  7. (se --pdb) coluna element do PDB == Type; e smoke-test de createSystem
     (com a FF completa) se openmm estiver disponivel  (bug: bonds nao casam)

Uso:
    python validate_ff.py docking_project/figures/smoke_test_fad/cofactor_FAD_TcTR_v3.xml \
        --pdb docking_project/figures/smoke_test_fad/cofactor_FAD_TcTR_v3.pdb
    python validate_ff.py lig_TcTR_smirnoff.xml          # so o XML

Sai com codigo 0 (PASS) ou 1 (FAIL). Use no pipeline ANTES da MD.
"""
from __future__ import annotations

import argparse
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

# Faixas PLAUSIVEIS de sigma LJ (nm) por elemento (GAFF/OpenFF, cobrindo
# alifatico+aromatico). Usamos FAIXA, nao "vizinho mais proximo", para nao dar
# falso-positivo em C aromatico (~0.348) que fica entre C e Cl/S.
# H: piso 0.0 — hidrogenios polares/ligados a N,O (e tipo agua) tem sigma LJ
# muito pequeno ou zero (ex.: 0.05 nm). So o TETO importa para pegar scramble
# (um H com sigma de C/N/O >=0.30 continua sinalizado).
SIGMA_RANGE = {"H": (0.0, 0.28), "C": (0.325, 0.367), "N": (0.30, 0.345),
               "O": (0.27, 0.317), "P": (0.355, 0.40), "S": (0.335, 0.375),
               "F": (0.285, 0.335), "Cl": (0.34, 0.36), "Br": (0.38, 0.41)}


def _implied_elem(name):
    for c in name:
        if c.isalpha():
            return c.upper()
    return None


def sigma_ok_for(elem, sigma):
    """True se sigma cabe na faixa plausivel do elemento (ou elemento desconhecido)."""
    rng = SIGMA_RANGE.get(elem)
    if rng is None:
        return True
    return rng[0] <= sigma <= rng[1]


def validate(xml_path: Path, pdb_path: Path | None):
    issues, warns, oks = [], [], []
    root = ET.parse(xml_path).getroot()

    # --- AtomTypes ---
    at = root.find("AtomTypes")
    types = {}
    if at is None or not len(at):
        issues.append("AtomTypes ausente/vazio")
    else:
        for t in at.findall("Type"):
            types[t.get("name")] = (t.get("element"), float(t.get("mass") or 0))
        oks.append(f"AtomTypes: {len(types)}")

    # --- Residue: atoms, names, bonds, charges ---
    resroot = root.find("Residues")
    res = resroot.find("Residue") if resroot is not None else None
    atom_type, atom_charge, res_bonds = {}, {}, []
    if res is None:
        issues.append("Residue ausente")
    else:
        for a in res.findall("Atom"):
            atom_type[a.get("name")] = a.get("type")
            if a.get("charge") is not None:
                atom_charge[a.get("name")] = float(a.get("charge"))
        res_bonds = [(b.get("atomName1"), b.get("atomName2")) for b in res.findall("Bond")]
        oks.append(f"Residue '{res.get('name')}': {len(atom_type)} atomos, {len(res_bonds)} bonds")

    # --- (1) Torsions ---
    tf = root.find("PeriodicTorsionForce")
    ntor = 0 if tf is None else len(list(tf))
    if ntor == 0:
        issues.append("PeriodicTorsionForce VAZIO (0 torsions) -> rings sem rigidez, NaN provavel")
    else:
        oks.append(f"PeriodicTorsionForce: {ntor} termos")

    # --- (2) Bonds (HBF) ---
    hbf = root.find("HarmonicBondForce")
    if hbf is None or not len(hbf):
        issues.append("HarmonicBondForce ausente/vazio (nenhuma mola harmonica)")
    else:
        ks = [float(b.get("k")) for b in hbf.findall("Bond") if b.get("k")]
        nzero = sum(1 for k in ks if k <= 0)
        if nzero:
            issues.append(f"HarmonicBondForce: {nzero} bond(s) com k<=0")
        oks.append(f"HarmonicBondForce: {len(ks)} bonds, k {min(ks):.0f}..{max(ks):.0f}")

    # --- (3) Angles ---
    af = root.find("HarmonicAngleForce")
    nang = 0 if af is None else len(list(af))
    if nang == 0:
        warns.append("HarmonicAngleForce vazio (sem termos de angulo)")
    else:
        oks.append(f"HarmonicAngleForce: {nang} angles")

    # --- (4) Cargas ---
    if atom_charge:
        s = sum(atom_charge.values())
        if abs(s - round(s)) > 0.02:
            warns.append(f"soma de cargas = {s:.3f} (longe de inteiro)")
        extremes = {n: q for n, q in atom_charge.items() if abs(q) > 2.0}
        if extremes:
            issues.append(f"cargas extremas |q|>2: {list(extremes.items())[:5]}")
        oks.append(f"cargas: soma={s:.3f}, |q|max={max(abs(q) for q in atom_charge.values()):.2f}")

    # --- (5) LJ sigma vs elemento do Type ---
    nb = root.find("NonbondedForce")
    if nb is not None and types:
        bad_lj = []
        for a in nb.findall("Atom"):
            tp = a.get("type") or a.get("class")
            sig = a.get("sigma")
            if tp in types and sig is not None:
                decl = types[tp][0]
                if decl and not sigma_ok_for(decl, float(sig)):
                    bad_lj.append((tp, decl, round(float(sig), 4)))
        if bad_lj:
            issues.append(f"LJ sigma fora da faixa do elemento em {len(bad_lj)} tipos "
                          f"(ex.: {bad_lj[0][0]} element={bad_lj[0][1]} sigma={bad_lj[0][2]} nm)")
        else:
            oks.append("LJ sigma coerente com os elementos")

    # --- (6) elemento <-> nome (so se nomes parecem quimicos) ---
    if atom_type and types:
        chem_names = [n for n in atom_type if _implied_elem(n) in SIGMA_RANGE]
        if len(chem_names) > 0.6 * len(atom_type):  # maioria tem nome quimico
            wrong = []
            for n in chem_names:
                decl = types.get(atom_type[n], (None,))[0]
                imp = _implied_elem(n)
                if decl and imp and decl != imp:
                    wrong.append((n, decl, imp))
            if wrong:
                issues.append(f"elemento != nome em {len(wrong)} atomos (scramble) "
                              f"(ex.: {wrong[0][0]} element={wrong[0][1]} mas nome sugere {wrong[0][2]})")
            else:
                oks.append("elemento <-> nome consistente")

    # --- (7) PDB ---
    if pdb_path:
        if not pdb_path.exists():
            warns.append(f"PDB {pdb_path.name} nao encontrado")
        else:
            mm = pdb_elem_mismatch(pdb_path, atom_type, types)
            if mm:
                issues.append(f"coluna element do PDB diverge do Type em {len(mm)} atomos "
                              f"(ex.: {mm[0]})")
            else:
                oks.append("coluna element do PDB consistente com os Types")
            smoke = createsystem_smoke(xml_path, pdb_path)
            if smoke is True:
                oks.append("createSystem (smoke-test): OK")
            elif smoke is not None:
                issues.append(f"createSystem falhou: {smoke}")

    return issues, warns, oks


def pdb_elem_mismatch(pdb_path, atom_type, types):
    out = []
    for ln in pdb_path.read_text(errors="replace").splitlines():
        if ln.startswith(("ATOM", "HETATM")):
            name = ln[12:16].strip()
            el = ln[76:78].strip()
            if name in atom_type and atom_type[name] in types:
                decl = types[atom_type[name]][0]
                if el and decl and el.capitalize() != decl.capitalize():
                    out.append(f"{name}: PDB={el} Type={decl}")
    return out


def createsystem_smoke(xml_path, pdb_path):
    """Tenta createSystem com a FF COMPLETA (amber14+tip3p+este XML). Retorna
    True (OK), None (openmm indisponivel) ou string de erro."""
    try:
        from openmm import app, unit
    except Exception:
        return None
    try:
        ff = app.ForceField("amber14-all.xml", "amber14/tip3p.xml", str(xml_path))
        pdb = app.PDBFile(str(pdb_path))
        ff.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff,
                        constraints=None, rigidWater=False)
        return True
    except Exception as e:
        return str(e)[:160]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("xml", type=Path, help="XML de ForceField OpenMM a validar")
    ap.add_argument("--pdb", type=Path, default=None, help="PDB template (opcional, recomendado)")
    args = ap.parse_args()
    if not args.xml.exists():
        print(f"FAIL: XML nao encontrado: {args.xml}")
        sys.exit(2)

    print(f"== validate_ff.py: {args.xml.name} ==")
    issues, warns, oks = validate(args.xml, args.pdb)
    for o in oks:
        print(f"  [ok]   {o}")
    for w in warns:
        print(f"  [warn] {w}")
    for i in issues:
        print(f"  [FAIL] {i}")
    print("-" * 60)
    if issues:
        print(f"RESULTADO: FAIL ({len(issues)} problema(s) critico(s)). NAO rode a MD ainda.")
        sys.exit(1)
    print("RESULTADO: PASS" + (f" ({len(warns)} aviso(s))" if warns else "") + ". Force-field OK p/ MD.")
    sys.exit(0)


if __name__ == "__main__":
    main()
