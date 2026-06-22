# -*- coding: utf-8 -*-
"""
prep_4new_receptor.py -- Prepara 4NEW (TcTR + inibidor EP127/2JR) para o redocking
de validacao: limpa (proteina + FAD; remove HOH e o ligante 2JR), converte para
receptor PDBQT (OpenBabel + Gasteiger, sem tags ROOT/BRANCH). Coloca 4NEW.pdb e
4NEW_receptor.pdbqt em docking_project/proteins/ para docking_rmsd_validation.py.
"""
import os, shutil
from openbabel import pybel, openbabel as ob

DOCK = r"C:\Users\luizd\OneDrive\Data_Science\Projetos\docking"
PROT = os.path.join(DOCK, "docking_project", "proteins")
VAL = os.path.join(DOCK, "docking_project", "validation")
SRC = os.path.join(VAL, "4NEW.pdb")
RAWDST = os.path.join(PROT, "4NEW.pdb")
CLEAN = os.path.join(PROT, "4NEW_clean.pdb")
RECPDBQT = os.path.join(PROT, "4NEW_receptor.pdbqt")

shutil.copy(SRC, RAWDST)

# 1) limpar: proteina (ATOM) + FAD; remover HOH e ligante 2JR
keep = []
for ln in open(SRC, errors="ignore"):
    if ln.startswith("ATOM"):
        keep.append(ln)
    elif ln.startswith("HETATM") and ln[17:20].strip() == "FAD":
        keep.append(ln)
    elif ln.startswith("TER"):
        keep.append(ln)
with open(CLEAN, "w") as fh:
    fh.writelines(keep); fh.write("END\n")
print("clean: %d linhas (ATOM+FAD)" % len(keep))

# 2) PDB -> PDBQT receptor (OpenBabel + Gasteiger)
mol = next(pybel.readfile("pdb", CLEAN))
mol.addh()
cm = ob.OBChargeModel.FindType("gasteiger")
if cm:
    cm.ComputeCharges(mol.OBMol)
tmp = RECPDBQT + ".tmp"
mol.write("pdbqt", tmp, overwrite=True)

bad = ("ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF")
kept = 0
with open(tmp) as fin, open(RECPDBQT, "w") as fout:
    for line in fin:
        if line.strip().startswith(bad):
            continue
        fout.write(line); kept += 1
os.remove(tmp)
print("receptor PDBQT: %d linhas -> %s" % (kept, RECPDBQT))
