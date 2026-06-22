#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
check_ligand_binding.py - PRE-VOO: o ligante esta ANCORADO na proteina?

Antes de gastar horas de GPU numa MD longa, confirma que o LIG (composto) esta
no BOLSO da proteina (em contato), e nao solto no solvente / fora do sitio.
Mede: distancia minima LIG<->proteina, atomos/residuos de contato (<4.5 A),
distancia ao FAD (cofator do sitio ativo da TR), e da um veredito.

Uso:
    python check_ligand_binding.py                       # complex_fresh.pdb
    python check_ligand_binding.py docking_project/figures/md_clean/minimized_prod_v3.pdb
    python check_ligand_binding.py <complexo.pdb> --lig LIG --cut 4.5
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path
import numpy as np

WORK = Path("docking_project/figures/md_clean")
AA = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
      "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
      "HID", "HIE", "HIP", "CYX"}
SOLV = {"HOH", "WAT", "NA", "CL", "K"}


def parse(pdb):
    lig, prot, fad = [], [], []   # (xyz)  ; prot/fad guardam (resname,resid,chain,name,xyz)
    for ln in pdb.read_text(errors="replace").splitlines():
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        resn = ln[17:20].strip()
        name = ln[12:16].strip()
        el = (ln[76:78].strip() or name[:1]).upper()
        if el == "H" or name.startswith("H"):
            continue
        try:
            xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
        except ValueError:
            continue
        chain = ln[21]
        resid = ln[22:26].strip()
        if resn == args_lig:
            lig.append(xyz)
        elif resn in AA:
            prot.append((resn, resid, chain, name, xyz))
        elif resn == "FAD":
            fad.append(xyz)
    return lig, prot, fad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pdb", nargs="?", default=str(WORK / "complex_fresh.pdb"))
    ap.add_argument("--lig", default="LIG")
    ap.add_argument("--cut", type=float, default=4.5, help="corte de contato (A)")
    a = ap.parse_args()
    global args_lig
    args_lig = a.lig

    pdb = Path(a.pdb)
    if not pdb.exists():
        sys.exit(f"PDB nao encontrado: {pdb}")
    print(f"== check_ligand_binding: {pdb.name} (ligante={a.lig}, corte={a.cut} A) ==")

    lig, prot, fad = parse(pdb)
    if not lig:
        sys.exit(f"FAIL: nenhum atomo do ligante '{a.lig}' encontrado no PDB.")
    if not prot:
        sys.exit("FAIL: nenhum atomo de proteina encontrado.")
    L = np.array(lig)
    P = np.array([p[4] for p in prot])
    print(f"  ligante: {len(lig)} heavy atoms | proteina: {len(prot)} heavy atoms")

    # distancias LIG x proteina
    D = np.sqrt(((L[:, None, :] - P[None, :, :]) ** 2).sum(-1))   # (nLig, nProt)
    dmin = float(D.min())
    contact_atoms = int((D < a.cut).any(axis=0).sum())
    # residuos em contato (qualquer atomo do residuo < corte de qualquer atomo do LIG)
    res_in_contact = {}
    pmin_per_atom = D.min(axis=0)   # menor dist de cada atomo de proteina ao LIG
    for k, (resn, resid, chain, name, _) in enumerate(prot):
        if pmin_per_atom[k] < a.cut:
            key = f"{resn}{resid}/{chain}"
            res_in_contact[key] = min(res_in_contact.get(key, 9e9), float(pmin_per_atom[k]))
    res_sorted = sorted(res_in_contact.items(), key=lambda x: x[1])

    print(f"  dist minima LIG<->proteina : {dmin:.2f} A")
    print(f"  atomos de proteina < {a.cut} A: {contact_atoms}")
    print(f"  residuos em contato        : {len(res_sorted)}")
    if res_sorted:
        print("    " + ", ".join(f"{r}({d:.1f})" for r, d in res_sorted[:18]))

    # dist ao FAD (sitio ativo da TR fica perto do flavin)
    if fad:
        F = np.array(fad)
        dfad = float(np.sqrt(((L[:, None, :] - F[None, :, :]) ** 2).sum(-1)).min())
        print(f"  dist minima LIG<->FAD      : {dfad:.2f} A")

    # veredito (contato vdW de atomos pesados nao-ligados e ~2.8-3.5 A; <2.0 = clash)
    n_res = len(res_sorted)
    print("-" * 56)
    if dmin < 2.0:
        print(f"VEREDITO: possivel CLASH (dist min {dmin:.2f} A < 2.0). "
              "Minimize melhor antes de rodar.")
        sys.exit(2)
    if contact_atoms >= 12 and n_res >= 5 and dmin <= 3.8:
        print(f"VEREDITO: LIGANTE ANCORADO NO BOLSO ({n_res} residuos em contato, "
              f"{contact_atoms} atomos; dist min {dmin:.2f} A = contato vdW normal). "
              "OK p/ rodar a MD longa.")
        sys.exit(0)
    if n_res >= 2 and dmin <= 4.5:
        print(f"VEREDITO: contato fraco/superficial ({n_res} residuos, dist min "
              f"{dmin:.2f} A). Revise a pose antes de 50 ns (pode dessorver).")
        sys.exit(1)
    print(f"VEREDITO: LIGANTE DESTACADO (dist min {dmin:.1f} A, {n_res} residuos) "
          "— NAO rode 50 ns; a pose esta fora do bolso.")
    sys.exit(2)


if __name__ == "__main__":
    main()
