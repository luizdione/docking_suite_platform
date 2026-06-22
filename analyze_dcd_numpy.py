#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
analyze_dcd_numpy.py - analise de trajetoria MD SEM mdtraj (so numpy+matplotlib).

Fallback leve para reproduzir md_clean.analyze() quando mdtraj nao esta instalado
(ex.: ambiente minimo, sandbox sem espaco para scipy). Le um DCD CHARMM/X-PLOR
('CORD', formato que o OpenMM DCDReporter grava) + o PDB de topologia e calcula:

  - RMSD da proteina (backbone CA), superposicao de Kabsch frame-a-frame vs frame 0;
  - RMSD do LIGANTE NO BOLSO: apos superpor a proteina, deslocamento do LIG vs
    frame 0 SEM re-superpor no proprio LIG -> capta deriva/saida do sitio;
  - RMSF por residuo (CA);
  - convergencia (media cumulativa) + estatistica de equilibrio (2a metade).

Robustez: se o PDB tiver MENOS atomos que o DCD (ex.: builds com numero de aguas
diferente, ou sync do OneDrive atrasado), usa o PREFIXO comum -- proteina+FAD+LIG
+ions sao um prefixo deterministico identico entre builds; so as aguas finais
variam, e CA/LIG vivem nesse prefixo. Salva PNG (3 paineis) + JSON com os MESMOS
nomes de md_clean.analyze() (md_<tag>_rmsd_rmsf.png / md_<tag>_rmsd_summary.json).

Uso:
    python analyze_dcd_numpy.py md_prod_v3.dcd minimized_prod_v3.pdb --tag prod_v3
    python analyze_dcd_numpy.py <dcd> <topo.pdb> --tag <tag> --out <dir> --lig LIG
"""
from __future__ import annotations
import argparse, json, struct, sys
from pathlib import Path
import numpy as np

AA = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
      "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
      "HID", "HIE", "HIP", "CYX"}


def read_dcd(path: Path) -> np.ndarray:
    """Le DCD CHARMM/X-PLOR. Retorna xyz[nframes, natom, 3] em Angstrom."""
    data = Path(path).read_bytes()
    p = 0

    def u(fmt):
        nonlocal p
        sz = struct.calcsize(fmt)
        v = struct.unpack(fmt, data[p:p + sz]); p += sz
        return v

    if u("<i")[0] != 84:
        raise ValueError("bloco inicial != 84 (nao parece DCD)")
    if data[p:p + 4] != b"CORD":
        raise ValueError("magic != CORD")
    p += 4
    icntrl = struct.unpack("<20i", data[p:p + 80]); p += 80
    nset = icntrl[0]
    has_cell = icntrl[10]                 # ICNTRL(11)=QCRYS: 1 -> unit cell por frame
    assert u("<i")[0] == 84
    tblk = u("<i")[0]; ntitle = u("<i")[0]; p += 80 * ntitle
    assert u("<i")[0] == tblk
    assert u("<i")[0] == 4
    natom = u("<i")[0]
    assert u("<i")[0] == 4

    frames = []
    rec = 4 * natom
    while p < len(data):
        if has_cell:
            cs = u("<i")[0]; p += cs; u("<i")        # 6 doubles (cell)
        coords = np.empty((natom, 3), dtype=np.float32)
        ok = True
        for d in range(3):
            try:
                s = u("<i")[0]
            except struct.error:
                ok = False; break
            if s != rec:
                ok = False; break
            coords[:, d] = np.frombuffer(data[p:p + rec], dtype="<f4"); p += rec
            if u("<i")[0] != rec:
                ok = False; break
        if not ok:
            break
        frames.append(coords)
    return np.asarray(frames)


def parse_pdb(path: Path):
    names, resn, resseq = [], [], []
    for ln in Path(path).read_text(errors="replace").splitlines():
        if ln.startswith(("ATOM", "HETATM")):
            names.append(ln[12:16].strip())
            resn.append(ln[17:20].strip())
            resseq.append(ln[22:26].strip())
    return np.array(names), np.array(resn), resseq


def superpose_on(frame, ref, sel):
    """Alinha `frame` em `ref` usando os atomos `sel` (Kabsch)."""
    pm = frame[sel].mean(0); qm = ref[sel].mean(0)
    H = (frame[sel] - pm).T @ (ref[sel] - qm)
    V, S, Wt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Wt.T @ V.T))
    R = Wt.T @ np.diag([1, 1, d]) @ V.T
    return (frame - pm) @ R.T + qm


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dcd", type=Path)
    ap.add_argument("topo", type=Path)
    ap.add_argument("--tag", default="prod")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--lig", default="LIG")
    a = ap.parse_args()
    out = a.out or a.dcd.parent

    print(f"== analyze_dcd_numpy: {a.dcd.name} | topo={a.topo.name} ==")
    xyz = read_dcd(a.dcd)
    n = len(xyz)
    names, resn, resseq = parse_pdb(a.topo)
    natom_pdb = len(names)
    print(f"  frames: {n} | DCD atoms: {xyz.shape[1] if n else 0} | PDB atoms: {natom_pdb}")
    if n == 0:
        sys.exit("FAIL: 0 frames lidos do DCD.")
    dcd_natom = xyz.shape[1]
    if dcd_natom != natom_pdb:
        if natom_pdb <= dcd_natom:
            print(f"  [aviso] DCD>{natom_pdb}: usando prefixo comum "
                  f"({dcd_natom - natom_pdb} aguas extras no build do DCD; "
                  f"proteina+FAD+LIG intactos).")
            xyz = xyz[:, :natom_pdb, :]
        else:
            sys.exit(f"FAIL: PDB {natom_pdb} > DCD {dcd_natom}; topologia incompativel.")
    if n < 2:
        sys.exit(f"AVISO: so {n} frame(s); analise temporal pulada.")

    ca = np.where((names == "CA") & np.isin(resn, list(AA)))[0]
    lig = np.where(resn == a.lig)[0]
    print(f"  CA: {len(ca)} | LIG '{a.lig}': {len(lig)}")
    if len(ca) < 3:
        sys.exit("FAIL: poucos CA para superposicao.")

    ref = xyz[0]
    aligned = np.empty_like(xyz)
    for i in range(n):
        aligned[i] = superpose_on(xyz[i], ref, ca)

    rmsd_prot = np.sqrt(((aligned[:, ca, :] - ref[ca]) ** 2).sum(2).mean(1))
    rmsd_lig = None
    if len(lig):
        dd = aligned[:, lig, :] - ref[lig]
        rmsd_lig = np.sqrt((dd ** 2).sum(2).mean(1))

    cummean = np.cumsum(rmsd_prot) / np.arange(1, n + 1)
    win = max(1, n // 20)
    runmean = np.convolve(rmsd_prot, np.ones(win) / win, mode="same")
    half = n // 2
    eq_mean = float(rmsd_prot[half:].mean()); eq_std = float(rmsd_prot[half:].std())
    mean_ca = aligned[:, ca, :].mean(0)
    rmsf = np.sqrt(((aligned[:, ca, :] - mean_ca) ** 2).sum(2).mean(0))
    res_ids = [int(resseq[i]) for i in ca]
    t = np.arange(n, dtype=float)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(18, 5))
    ax[0].plot(t, rmsd_prot, color="#1f77b4", lw=0.9, alpha=0.6, label="Proteina (CA)")
    ax[0].plot(t, runmean, color="#08306b", lw=1.8, label=f"media movel ({win}f)")
    if rmsd_lig is not None:
        ax[0].plot(t, rmsd_lig, color="#d62728", lw=1.2, alpha=0.9, label="Ligante no bolso")
    ax[0].axvspan(float(t[half]), float(t[-1]), color="green", alpha=0.06)
    ax[0].set(xlabel="Frame (DCD)", ylabel="RMSD (A)",
              title=f"RMSD - {a.tag} (equil. prot {eq_mean:.2f}+/-{eq_std:.2f} A)")
    ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)
    ax[1].plot(res_ids, rmsf, color="#2ca02c")
    ax[1].set(xlabel="Residuo (resSeq)", ylabel="RMSF (A)", title=f"RMSF - {a.tag}")
    ax[1].grid(alpha=0.3)
    ax[2].plot(t, cummean, color="#6a3d9a", lw=1.8)
    ax[2].axhline(eq_mean, ls="--", color="gray", lw=1, label=f"equil {eq_mean:.2f} A")
    ax[2].set(xlabel="Frame (DCD)", ylabel="RMSD cumul. (A)", title=f"Convergencia - {a.tag}")
    ax[2].legend(fontsize=8); ax[2].grid(alpha=0.3)
    fig.tight_layout()
    out.mkdir(parents=True, exist_ok=True)
    png = out / f"md_{a.tag}_rmsd_rmsf.png"
    fig.savefig(png, dpi=200, bbox_inches="tight"); plt.close(fig)

    summary = {
        "tag": a.tag, "n_frames": int(n),
        "rmsd_prot_final_A": float(rmsd_prot[-1]),
        "rmsd_prot_equil_mean_A": eq_mean, "rmsd_prot_equil_std_A": eq_std,
        "rmsd_lig_pocket_final_A": (None if rmsd_lig is None else float(rmsd_lig[-1])),
        "rmsd_lig_pocket_max_A": (None if rmsd_lig is None else float(rmsd_lig.max())),
        "rmsf_mean_A": float(rmsf.mean()), "rmsf_max_A": float(rmsf.max()),
        "note": "gerado por analyze_dcd_numpy.py (sem mdtraj); eixo X = indice de frame",
    }
    (out / f"md_{a.tag}_rmsd_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"  figura: {png}")
    print(f"  resumo: {json.dumps(summary)}")

    if rmsd_lig is not None:
        lp = float(rmsd_lig.max())
        print("-" * 50)
        if lp < 2.0:
            print(f"VEREDITO binding: LIGANTE ESTAVEL no bolso (desvio max {lp:.2f} A).")
        elif lp < 3.5:
            print(f"VEREDITO binding: leve deriva (max {lp:.2f} A) - aceitavel.")
        else:
            print(f"VEREDITO binding: ATENCAO, ligante derivou {lp:.2f} A - rode mais ns.")


if __name__ == "__main__":
    main()
