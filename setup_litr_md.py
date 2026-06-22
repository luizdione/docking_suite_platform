# -*- coding: utf-8 -*-
"""
setup_litr_md.py -- Setup das DMs da LiTR (2JK6, L. infantum) para md_simulation.py
--target LiTR (modos subversive E inhibitor). Dimero com 2 FADs; ditiol Cys52/57.

Gera em docking_project/figures/smoke_test_litr/:
  1. cofactor_FAD_TcTR.pdb         (template canonico 86 atomos, copiado)
  2. cofactor_FAD_0_xtal.pdb       (FAD cadeia A do 2JK6, 53 heavy, ordem do template)
  3. cofactor_FAD_1_xtal.pdb       (FAD cadeia B do 2JK6)
  4. lig_LiTR_fadsite.pdb          (subversivo: 2-SCN-NQ na si-face do FAD da cadeia A,
                                    voltada p/ ditiol Cys52/57 da cadeia A)

(O modo inibidor usa lig_LiTR.pdb, gerado por build_litr_inhibitor_pose.py.)
FAD/LIG XML reaproveitados do TcTR pelo md_simulation.py. Nomes de FAD == template.

Uso (env datascience):  python setup_litr_md.py
"""
import os, shutil, numpy as np

DOCK = r"C:\Users\luizd\OneDrive\Data_Science\Projetos\docking"
FIG = os.path.join(DOCK, "docking_project", "figures")
PROT = os.path.join(DOCK, "docking_project", "proteins")
SMTC = os.path.join(FIG, "smoke_test_fad")
SMLI = os.path.join(FIG, "smoke_test_litr"); os.makedirs(SMLI, exist_ok=True)
TEMPLATE = os.path.join(SMTC, "cofactor_FAD_TcTR.pdb")
LIG_SRC = os.path.join(SMTC, "lig_TcTR.pdb")
RECEPTOR = os.path.join(PROT, "2JK6_clean.pdb")

DITHIOL_RESSEQ = (52, 57)
POSE_CHAIN = "A"               # protomero usado p/ semear a pose subversiva
STACK = 3.6
ISO_NAMES = {"N1", "C2", "N3", "C4", "C4A", "N5", "C5A", "C6", "C7", "C8", "C9", "C9A", "N10", "C10"}
NAPH_NAMES = {"LG%d" % i for i in range(12)}


def is_h(ln):
    nm = ln[12:16].strip(); el = ln[76:78].strip()
    return el == "H" or nm[:1] == "H"


def fad_heavy_by_chain(chain):
    d = {}
    for ln in open(RECEPTOR, errors="ignore"):
        if ln.startswith("HETATM") and ln[17:20].strip() == "FAD" and ln[21] == chain and not is_h(ln):
            d[ln[12:16].strip()] = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
    return d


def write_xtal(chain, out_name):
    coords = fad_heavy_by_chain(chain)
    out = []
    for ln in open(TEMPLATE, errors="ignore"):
        if not ln.startswith(("ATOM", "HETATM")) or is_h(ln):
            continue
        nm = ln[12:16].strip()
        if nm not in coords:
            raise SystemExit("nome %s do template ausente na cadeia %s do 2JK6" % (nm, chain))
        x, y, z = coords[nm]
        out.append(ln[:30] + "%8.3f%8.3f%8.3f" % (x, y, z) + ln[54:])
    if len(out) != 53:
        raise SystemExit("cadeia %s: esperados 53 heavy, %d" % (chain, len(out)))
    with open(os.path.join(SMLI, out_name), "w") as fh:
        fh.writelines(out); fh.write("END\n")
    print("  FAD cadeia %s -> %s" % (chain, out_name))


def read_atoms(path, resn=None, chain=None):
    out = []
    for ln in open(path, errors="ignore"):
        if not ln.startswith(("ATOM", "HETATM")): continue
        if resn and ln[17:20].strip() != resn: continue
        if chain and ln[21] != chain: continue
        out.append(dict(name=ln[12:16].strip(),
                        xyz=np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])))
    return out


def plane(c):
    cen = c.mean(0); _, _, vh = np.linalg.svd(c - cen); return cen, vh[2]


def rot_matrix(a, b):
    a = a / np.linalg.norm(a); b = b / np.linalg.norm(b)
    v = np.cross(a, b); s = np.linalg.norm(v); cth = np.dot(a, b)
    if s < 1e-8: return np.eye(3) if cth > 0 else -np.eye(3)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * ((1 - cth) / (s * s))


def build_pose():
    # ditiol catalitico Cys52/57 da cadeia escolhida
    sg = []
    for ln in open(RECEPTOR, errors="ignore"):
        if ln.startswith("ATOM") and ln[12:16].strip() == "SG" and ln[17:20].strip() == "CYS" \
                and ln[21] == POSE_CHAIN:
            try: rs = int(ln[22:26])
            except: continue
            if rs in DITHIOL_RESSEQ:
                sg.append(np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]))
    dithiol = np.mean(sg, axis=0)
    print("  ditiol LiTR Cys52/57 (cadeia %s) = %s" % (POSE_CHAIN, np.round(dithiol, 2)))

    fad = read_atoms(RECEPTOR, resn="FAD", chain=POSE_CHAIN)
    iso = [a for a in fad if a["name"] in ISO_NAMES and np.linalg.norm(a["xyz"] - dithiol) < 12.0]
    if len(iso) < 6:
        raise SystemExit("poucos atomos de isoaloxazina perto do ditiol (%d)" % len(iso))
    c_iso, n_iso = plane(np.array([a["xyz"] for a in iso]))
    if np.dot(dithiol - c_iso, n_iso) < 0: n_iso = -n_iso

    lig = read_atoms(LIG_SRC)
    naph = [a for a in lig if a["name"] in NAPH_NAMES]
    c_lig, n_lig = plane(np.array([a["xyz"] for a in naph]))
    R = rot_matrix(n_lig, n_iso); target = c_iso + STACK * n_iso
    # reescreve preservando linhas originais do template do ligante
    out_lines = []
    for ln in open(LIG_SRC, errors="ignore"):
        if ln.startswith(("ATOM", "HETATM")):
            p = np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
            q = R @ (p - c_lig) + target
            out_lines.append("%s%8.3f%8.3f%8.3f%s" % (ln[:30], q[0], q[1], q[2], ln[54:]))
    conect = [ln for ln in open(LIG_SRC, errors="ignore") if ln.startswith("CONECT")]
    p_out = os.path.join(SMLI, "lig_LiTR_fadsite.pdb")
    with open(p_out, "w") as fh:
        fh.writelines(out_lines); fh.writelines(conect); fh.write("END\n")
    nq_new = np.array([R @ (a["xyz"] - c_lig) + target for a in naph]).mean(0)
    print("  pose si-face -> lig_LiTR_fadsite.pdb | dist->isoaloxazina=%.2f A | CONECT=%d"
          % (np.linalg.norm(nq_new - c_iso), len(conect)))


if __name__ == "__main__":
    shutil.copy(TEMPLATE, os.path.join(SMLI, "cofactor_FAD_TcTR.pdb"))
    print("[1] template -> smoke_test_litr/cofactor_FAD_TcTR.pdb")
    print("[2] FAD xtais (2 cadeias):")
    write_xtal("A", "cofactor_FAD_0_xtal.pdb")
    write_xtal("B", "cofactor_FAD_1_xtal.pdb")
    print("[3] pose subversiva:")
    build_pose()
    print("\nSETUP LiTR (subversivo) OK.")
