# -*- coding: utf-8 -*-
"""
setup_hgr_md.py -- Prepara o setup da DM subversiva da hGR (controle de
seletividade) para md_simulation.py --target hGR --mode subversive --fresh.

Gera em docking_project/figures/smoke_test_hgr/:
  1. cofactor_FAD_TcTR.pdb       (template canonico 86 atomos, copiado)
  2. cofactor_FAD_0_xtal.pdb     (53 heavy do FAD do 3DJJ, na ORDEM/NOMES do
                                  template -> casa com build_canonical_fad_topology)
  3. lig_hGR_fadsite.pdb         (2-SCN-NQ empilhada na si-face da isoaloxazina do
                                  FAD do 3DJJ, voltada para o ditiol Cys58/Cys63)

FAD XML e LIG SMIRNOFF sao reaproveitados do TcTR pelo proprio md_simulation.py
(mesmo cofator/molecula). Nomes de atomo do FAD 3DJJ == template (conferido).

Uso (env datascience):  python setup_hgr_md.py
"""
import os, shutil, numpy as np

DOCK = r"C:\Users\luizd\OneDrive\Data_Science\Projetos\docking"
FIG  = os.path.join(DOCK, "docking_project", "figures")
PROT = os.path.join(DOCK, "docking_project", "proteins")
SMTC = os.path.join(FIG, "smoke_test_fad")      # templates TcTR
SMHG = os.path.join(FIG, "smoke_test_hgr")      # destino hGR
os.makedirs(SMHG, exist_ok=True)

TEMPLATE = os.path.join(SMTC, "cofactor_FAD_TcTR.pdb")   # 86 atomos canonicos
LIG_SRC  = os.path.join(SMTC, "lig_TcTR.pdb")            # nomes/CONECT do ligante
RECEPTOR = os.path.join(PROT, "3DJJ_clean.pdb")          # hGR + FAD + Cys58/63

FAD_RESID = " 479"                                        # FAD do 3DJJ (chain A)
DITHIOL_RESSEQ = (58, 63)                                 # Cys catalitica hGR
STACK = 3.6                                               # A, empilhamento pi
ISO_NAMES = {"N1","C2","N3","C4","C4A","N5","C5A","C6","C7","C8","C9","C9A","N10","C10"}
NAPH_NAMES = {"LG%d" % i for i in range(12)}


def is_h(ln):
    nm = ln[12:16].strip(); el = ln[76:78].strip()
    return el == "H" or (nm[:1] == "H")


def fad_heavy_3djj():
    """{name: (x,y,z)} dos 53 heavy do FAD 479 do 3DJJ."""
    d = {}
    for ln in open(RECEPTOR, errors="ignore"):
        if ln.startswith("HETATM") and ln[17:20].strip() == "FAD" and ln[22:26] == FAD_RESID:
            if is_h(ln): continue
            d[ln[12:16].strip()] = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
    return d


def step_template():
    shutil.copy(TEMPLATE, os.path.join(SMHG, "cofactor_FAD_TcTR.pdb"))
    print("[1] template canonico copiado -> smoke_test_hgr/cofactor_FAD_TcTR.pdb")


def step_xtal():
    """Reescreve as 53 linhas heavy do TEMPLATE com as coordenadas do FAD do 3DJJ
    (mesmos nomes/colunas do template -> match garantido)."""
    coords = fad_heavy_3djj()
    out, miss = [], []
    for ln in open(TEMPLATE, errors="ignore"):
        if not ln.startswith(("ATOM", "HETATM")) or is_h(ln):
            continue
        nm = ln[12:16].strip()
        if nm not in coords:
            miss.append(nm); continue
        x, y, z = coords[nm]
        out.append(ln[:30] + "%8.3f%8.3f%8.3f" % (x, y, z) + ln[54:])
    if miss:
        raise SystemExit("ERRO: nomes do template ausentes no 3DJJ: %s" % miss)
    if len(out) != 53:
        raise SystemExit("ERRO: esperados 53 heavy, obtidos %d" % len(out))
    p = os.path.join(SMHG, "cofactor_FAD_0_xtal.pdb")
    with open(p, "w") as fh:
        fh.writelines(out); fh.write("END\n")
    print("[2] xtal FAD hGR (53 heavy, ordem do template) -> %s" % os.path.basename(p))


def read_atoms(path, resn=None):
    out = []
    for ln in open(path, errors="ignore"):
        if not ln.startswith(("ATOM", "HETATM")): continue
        if resn and ln[17:20].strip() != resn: continue
        out.append(dict(line=ln, name=ln[12:16].strip(),
                        xyz=np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])))
    return out


def plane(coords):
    c = coords.mean(0)
    _, _, vh = np.linalg.svd(coords - c)
    return c, vh[2]


def rot_matrix(a, b):
    a = a / np.linalg.norm(a); b = b / np.linalg.norm(b)
    v = np.cross(a, b); s = np.linalg.norm(v); cth = np.dot(a, b)
    if s < 1e-8:
        return np.eye(3) if cth > 0 else -np.eye(3)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * ((1 - cth) / (s * s))


def step_pose():
    # ditiol catalitico hGR: ponto medio dos SG de Cys58/Cys63 (3DJJ)
    sg = []
    for ln in open(RECEPTOR, errors="ignore"):
        if ln.startswith("ATOM") and ln[12:16].strip() == "SG" and ln[17:20].strip() == "CYS":
            try: rs = int(ln[22:26])
            except: continue
            if rs in DITHIOL_RESSEQ:
                sg.append(np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]))
    if len(sg) < 2:
        raise SystemExit("ERRO: SG de Cys58/Cys63 nao encontrados no 3DJJ (achei %d)" % len(sg))
    dithiol = np.mean(sg, axis=0)
    print("    ditiol hGR (Cys58/63) centro = %s" % np.round(dithiol, 2))

    # isoaloxazina do FAD do 3DJJ proxima ao ditiol
    fad = read_atoms(RECEPTOR, resn="FAD")
    iso = [a for a in fad if a["name"] in ISO_NAMES and np.linalg.norm(a["xyz"] - dithiol) < 12.0]
    if len(iso) < 6:
        raise SystemExit("Poucos atomos de isoaloxazina perto do ditiol (%d)" % len(iso))
    iso_xyz = np.array([a["xyz"] for a in iso])
    c_iso, n_iso = plane(iso_xyz)
    if np.dot(dithiol - c_iso, n_iso) < 0:
        n_iso = -n_iso
    print("    isoaloxazina: %d atomos, centroide=%s" % (len(iso), np.round(c_iso, 2)))

    # ligante (mesmos nomes/ordem do lig_TcTR -> casa com SMIRNOFF reaproveitado)
    lig = read_atoms(LIG_SRC)
    naph = [a for a in lig if a["name"] in NAPH_NAMES]
    naph_xyz = np.array([a["xyz"] for a in naph])
    c_lig, n_lig = plane(naph_xyz)
    R = rot_matrix(n_lig, n_iso)
    target = c_iso + STACK * n_iso

    new_lines = []
    for a in lig:
        new = R @ (a["xyz"] - c_lig) + target
        ln = a["line"]
        new_lines.append("%s%8.3f%8.3f%8.3f%s" % (ln[:30], new[0], new[1], new[2], ln[54:]))
    conect = [ln for ln in open(LIG_SRC, errors="ignore") if ln.startswith("CONECT")]

    p = os.path.join(SMHG, "lig_hGR_fadsite.pdb")
    with open(p, "w") as fh:
        fh.writelines(new_lines); fh.writelines(conect); fh.write("END\n")
    nq_new = np.array([R @ (a["xyz"] - c_lig) + target for a in naph]).mean(0)
    print("[3] pose si-face -> %s | dist naftoquinona->isoaloxazina=%.2f A | CONECT=%d"
          % (os.path.basename(p), np.linalg.norm(nq_new - c_iso), len(conect)))


if __name__ == "__main__":
    step_template()
    step_xtal()
    step_pose()
    print("\nSETUP hGR OK. Proximo: smoke test curto (md_simulation.py --target hGR "
          "--mode subversive --fresh --ns 0.1) antes do run de 50 ns.")
