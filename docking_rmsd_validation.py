#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
docking_rmsd_validation.py - Validacao do protocolo de docking por REDOCKING.

Torna REAL a afirmacao do Metodos ("re-docking de ligante co-cristalizado,
RMSD < 2.0 A") que hoje nao tem codigo por tras. Protocolo classico de
validacao: extrai o ligante co-cristalizado, re-docka na MESMA estrutura com a
caixa centrada nele, e mede o RMSD heavy-atom da melhor pose vs. o cristal.

Caso TcTR: 1BZL tem GCG (glutationa) co-cristalizada -> referencia.

Roda no env datascience (precisa: vina.exe na pasta do projeto; openbabel/obabel
para preparar pdbqt). Para outros alvos, ajuste PDB/ligand code.

Uso:
    python docking_rmsd_validation.py                 # 1BZL / GCG
    python docking_rmsd_validation.py --pdb 1BZL --lig GCG --box 22
"""
from __future__ import annotations
import argparse, subprocess, shutil, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
PROT = HERE / "docking_project" / "proteins"
VINA = HERE / "vina.exe"


def _first_chain(pdb: Path, resname: str):
    """Primeira cadeia (e resSeq) que contem o ligante -> evita pegar as 2 copias
    do dimero (que viram 2 ROOT no pdbqt e quebram o Vina)."""
    for ln in pdb.read_text(errors="replace").splitlines():
        if ln.startswith("HETATM") and ln[17:20].strip() == resname:
            return ln[21], ln[22:26]
    return None, None


def read_hetatm(pdb: Path, resname: str):
    """Coords + nomes do ligante co-cristalizado (heavy atoms), 1a copia apenas."""
    ch, rs = _first_chain(pdb, resname)
    names, xyz = [], []
    for ln in pdb.read_text(errors="replace").splitlines():
        if ln.startswith("HETATM") and ln[17:20].strip() == resname \
                and ln[21] == ch and ln[22:26] == rs:
            el = ln[76:78].strip()
            nm = ln[12:16].strip()
            if el.upper() == "H" or nm.startswith("H"):
                continue
            names.append(nm)
            xyz.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
    return names, np.array(xyz)


def write_lig_pdb(pdb: Path, resname: str, out: Path):
    ch, rs = _first_chain(pdb, resname)
    lines = [ln for ln in pdb.read_text(errors="replace").splitlines()
             if ln.startswith("HETATM") and ln[17:20].strip() == resname
             and ln[21] == ch and ln[22:26] == rs]
    out.write_text("\n".join(lines) + "\nEND\n")


def run(cmd):
    print("   $", " ".join(str(c) for c in cmd))
    return subprocess.run([str(c) for c in cmd], capture_output=True, text=True)


def rmsd_named(ref_names, ref_xyz, pose_names, pose_xyz):
    """RMSD heavy-atom emparelhando por NOME (sem re-superposicao: mesmo
    referencial do receptor)."""
    idx = {n: i for i, n in enumerate(pose_names)}
    common = [n for n in ref_names if n in idx]
    if len(common) < 0.5 * len(ref_names):
        # fallback: assume mesma ordem
        m = min(len(ref_xyz), len(pose_xyz))
        d = ref_xyz[:m] - pose_xyz[:m]
        return float(np.sqrt((d * d).sum(1).mean())), m
    a = np.array([ref_xyz[ref_names.index(n)] for n in common])
    b = np.array([pose_xyz[idx[n]] for n in common])
    d = a - b
    return float(np.sqrt((d * d).sum(1).mean())), len(common)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", default="1BZL")
    ap.add_argument("--lig", default="GCG", help="codigo do ligante co-cristalizado")
    ap.add_argument("--box", type=float, default=22.0, help="aresta da caixa (A)")
    ap.add_argument("--exhaustiveness", type=int, default=32)
    args = ap.parse_args()

    raw = PROT / f"{args.pdb}.pdb"
    receptor_pdbqt = PROT / f"{args.pdb}_receptor.pdbqt"
    if not raw.exists() or not receptor_pdbqt.exists():
        sys.exit(f"FALTA {raw.name} ou {receptor_pdbqt.name} em {PROT}")
    if not VINA.exists():
        sys.exit(f"vina.exe nao encontrado em {VINA}")
    obabel = shutil.which("obabel") or shutil.which("obabel.exe")
    if not obabel:
        sys.exit("openbabel (obabel) nao encontrado no PATH (conda install -c conda-forge openbabel)")

    work = HERE / "docking_project" / "validation"
    work.mkdir(parents=True, exist_ok=True)

    # 1) referencia cristalografica
    ref_names, ref_xyz = read_hetatm(raw, args.lig)
    if len(ref_xyz) == 0:
        sys.exit(f"ligante {args.lig} nao encontrado em {raw.name}")
    center = ref_xyz.mean(0)
    print(f"Ref {args.lig}: {len(ref_names)} heavy atoms; centro caixa = "
          f"{center[0]:.2f} {center[1]:.2f} {center[2]:.2f}")

    # 2) preparar ligante pdbqt (a partir do cristal, com Hs e cargas)
    lig_pdb = work / f"{args.lig}_xtal.pdb"
    write_lig_pdb(raw, args.lig, lig_pdb)
    lig_pdbqt = work / f"{args.lig}_lig.pdbqt"
    run([obabel, lig_pdb, "-O", lig_pdbqt, "-p", "7.4", "--partialcharge", "gasteiger"])

    # 3) vina redocking (caixa centrada no ligante cristalografico)
    out_pdbqt = work / f"redock_{args.pdb}_{args.lig}.pdbqt"
    r = run([VINA, "--receptor", receptor_pdbqt, "--ligand", lig_pdbqt,
             "--center_x", center[0], "--center_y", center[1], "--center_z", center[2],
             "--size_x", args.box, "--size_y", args.box, "--size_z", args.box,
             "--exhaustiveness", args.exhaustiveness, "--out", out_pdbqt])
    if r.returncode != 0:
        print(r.stdout[-1000:]); sys.exit("vina falhou")

    # 4) melhor pose -> pdb -> RMSD vs cristal
    best_pdb = work / f"redock_{args.pdb}_{args.lig}_best.pdb"
    run([obabel, out_pdbqt, "-O", best_pdb, "-f", "1", "-l", "1"])
    pose_names, pose_xyz = [], []
    for ln in best_pdb.read_text(errors="replace").splitlines():
        if ln.startswith(("ATOM", "HETATM")):
            el = ln[76:78].strip(); nm = ln[12:16].strip()
            if el.upper() == "H" or nm.startswith("H"):
                continue
            pose_names.append(nm)
            pose_xyz.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
    pose_xyz = np.array(pose_xyz)

    rmsd, n = rmsd_named(ref_names, ref_xyz, pose_names, pose_xyz)
    print("\n" + "=" * 56)
    print(f"REDOCKING {args.pdb}/{args.lig}: RMSD melhor pose vs cristal = "
          f"{rmsd:.2f} A  ({n} atomos)")
    print("  -> protocolo VALIDADO" if rmsd < 2.0 else
          "  -> RMSD > 2 A: ligante flexivel/box; reportar e discutir")
    print("=" * 56)
    (work / f"validation_{args.pdb}_{args.lig}.txt").write_text(
        f"redocking {args.pdb}/{args.lig}\nRMSD={rmsd:.3f} A over {n} heavy atoms\n"
        f"box center={center.tolist()} edge={args.box}\n")


if __name__ == "__main__":
    main()
