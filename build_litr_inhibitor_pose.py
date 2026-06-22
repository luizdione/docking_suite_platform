# -*- coding: utf-8 -*-
"""
build_litr_inhibitor_pose.py -- Semente do MODO INIBIDOR da LiTR a partir da pose
de docking Vina (result_LiTR.pdbqt). Mesma estrategia rigida do hGR: casamento por
impressao digital de distancias (Hungarian, por elemento) + Kabsch do template
lig_TcTR.pdb. Saida: smoke_test_litr/lig_LiTR.pdb.

Uso (env datascience):  python build_litr_inhibitor_pose.py
"""
import os, numpy as np
from scipy.optimize import linear_sum_assignment

DOCK = r"C:\Users\luizd\OneDrive\Data_Science\Projetos\docking"
FIG = os.path.join(DOCK, "docking_project", "figures")
SMLI = os.path.join(FIG, "smoke_test_litr"); os.makedirs(SMLI, exist_ok=True)
LIG_TMPL = os.path.join(FIG, "smoke_test_fad", "lig_TcTR.pdb")
PDBQT = os.path.join(DOCK, "docking_project", "docking", "result_LiTR.pdbqt")
PROT = os.path.join(DOCK, "docking_project", "proteins", "2JK6_clean.pdb")
AD2EL = {"A": "C", "C": "C", "N": "N", "NA": "N", "O": "O", "OA": "O", "S": "S", "SA": "S"}
DITHIOL_RESSEQ = (52, 57); POSE_CHAIN = "A"


def read_heavy(path):
    els, xyz = [], []
    for ln in open(path, errors="ignore"):
        if ln.startswith(("ATOM", "HETATM")):
            nm = ln[12:16].strip(); el = (ln[76:78].strip() or nm[:1])
            if el.upper() == "H" or nm[:1] == "H": continue
            els.append(el.capitalize()); xyz.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
    return els, np.array(xyz)


def read_vina_heavy():
    els, xyz = [], []
    for ln in open(PDBQT, errors="ignore"):
        if ln.startswith("ENDMDL"): break
        if ln.startswith(("ATOM", "HETATM")):
            adt = ln[77:79].strip(); el = AD2EL.get(adt, (adt[:1] or ln[12:14].strip()[:1]))
            if el == "H": continue
            els.append(el.capitalize()); xyz.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
    return els, np.array(xyz)


def fp(xyz):
    return np.sort(np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=-1), axis=1)


def correspondence(elsA, A, elsB, B):
    fA, fB = fp(A), fp(B); corr = [-1] * len(elsA)
    for el in set(elsA):
        ia = [i for i in range(len(elsA)) if elsA[i] == el]
        ib = [j for j in range(len(elsB)) if elsB[j] == el]
        C = np.array([[np.linalg.norm(fA[i] - fB[j]) for j in ib] for i in ia])
        r, c = linear_sum_assignment(C)
        for k in range(len(r)): corr[ia[r[k]]] = ib[c[k]]
    return corr


def kabsch(P, Q):
    cP, cQ = P.mean(0), Q.mean(0); H = (P - cP).T @ (Q - cQ)
    U, S, Vt = np.linalg.svd(H); d = np.sign(np.linalg.det(Vt.T @ U.T))
    return Vt.T @ np.diag([1, 1, d]) @ U.T, cP, cQ


def main():
    elsT, T = read_heavy(LIG_TMPL); elsV, V = read_vina_heavy()
    print("template heavy=%d  vina heavy=%d" % (len(elsT), len(elsV)))
    corr = correspondence(elsT, T, elsV, V)
    if any(c < 0 for c in corr): raise SystemExit("correspondencia incompleta")
    Vm = V[corr]; R, cP, cQ = kabsch(T, Vm)
    rmsd = np.sqrt(np.mean(np.sum(((T - cP) @ R + cQ - Vm) ** 2, axis=1)))
    print("RMSD rigido template->Vina = %.2f A" % rmsd)
    out_lines, conect = [], []
    for ln in open(LIG_TMPL, errors="ignore"):
        if ln.startswith(("ATOM", "HETATM")):
            p = np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
            q = (p - cP) @ R + cQ
            out_lines.append("%s%8.3f%8.3f%8.3f%s" % (ln[:30], q[0], q[1], q[2], ln[54:]))
        elif ln.startswith("CONECT"):
            conect.append(ln)
    with open(os.path.join(SMLI, "lig_LiTR.pdb"), "w") as fh:
        fh.writelines(out_lines); fh.writelines(conect); fh.write("END\n")
    # diagnostico FAD/ditiol (cadeia A)
    iso_names = {"N1", "C2", "N3", "C4", "C4A", "N5", "C5A", "C6", "C7", "C8", "C9", "C9A", "N10", "C10"}
    iso, sg = [], []
    for ln in open(PROT, errors="ignore"):
        if ln.startswith("HETATM") and ln[17:20].strip() == "FAD" and ln[21] == POSE_CHAIN and ln[12:16].strip() in iso_names:
            iso.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
        if ln.startswith("ATOM") and ln[12:16].strip() == "SG" and ln[17:20].strip() == "CYS" and ln[21] == POSE_CHAIN:
            try: rs = int(ln[22:26])
            except: continue
            if rs in DITHIOL_RESSEQ: sg.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
    Lh = (T - cP) @ R + cQ; iso = np.array(iso); dit = np.array(sg).mean(0)
    print("OK -> lig_LiTR.pdb  | pose: LIG->FAD(min)=%.1f A  LIG->ditiol(min)=%.1f A"
          % (np.min(np.linalg.norm(Lh[:, None] - iso[None], axis=-1)), np.min(np.linalg.norm(Lh - dit, axis=-1))))


if __name__ == "__main__":
    main()
