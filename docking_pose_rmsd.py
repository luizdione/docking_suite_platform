#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
docking_pose_rmsd.py - Analise de RMSD das poses do AutoDock Vina.

Responde a critica "nao usa grafico de RMSD" no LADO DOCKING (a parte de MD tem
RMSD/RMSF em md_clean.py::analyze()). Para cada result_*.pdbqt:
  - le todas as poses (MODEL) e o afinidade/rmsd do Vina (REMARK VINA RESULT);
  - calcula a matriz de RMSD pareado heavy-atom (no referencial do receptor,
    sem re-superposicao -- convencao AutoDock para clustering de poses);
  - agrupa por single-linkage num corte (default 2.0 A);
  - gera figura (afinidade x pose por cluster, RMSD vs melhor pose, heatmap) e
    um resumo CSV/JSON.

NAO precisa de RDKit/OpenFF -- so numpy + matplotlib.

Uso:
    python docking_pose_rmsd.py                 # todos os result_*.pdbqt
    python docking_pose_rmsd.py --target TcTR    # so um alvo
    python docking_pose_rmsd.py --cutoff 2.0
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
DOCK_DIR = HERE / "docking_project" / "docking"
FIG_DIR = HERE / "docking_project" / "figures"


def parse_pdbqt_poses(path: Path):
    """Retorna (list_of_coord_arrays[Nx3], list_of_atom_names, list_of_(affinity,
    rmsd_lb, rmsd_ub))."""
    poses, names, vina = [], None, []
    cur, cur_names = None, None
    for ln in path.read_text(errors="replace").splitlines():
        if ln.startswith("MODEL"):
            cur, cur_names = [], []
        elif ln.startswith("REMARK") and "VINA RESULT" in ln:
            nums = re.findall(r"[-+]?\d+\.\d+", ln)
            vina.append(tuple(float(x) for x in nums[:3]) if len(nums) >= 3
                        else (float(nums[0]) if nums else np.nan, np.nan, np.nan))
        elif ln.startswith(("ATOM", "HETATM")) and cur is not None:
            elem = ln[76:78].strip()
            name = ln[12:16].strip()
            # pdbqt: ignora H (apolar ja embutido); mantem heavy + polar
            if elem.upper() == "H" or name.startswith("H"):
                continue
            cur.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
            cur_names.append(name)
        elif ln.startswith("ENDMDL") and cur is not None:
            poses.append(np.array(cur))
            if names is None:
                names = cur_names
            cur, cur_names = None, None
    return poses, names, vina


def rmsd_matrix(poses):
    n = len(poses)
    M = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = poses[i] - poses[j]
            M[i, j] = M[j, i] = np.sqrt((d * d).sum(axis=1).mean())
    return M


def single_linkage(M, cutoff):
    """Clustering single-linkage simples (union-find) no corte dado."""
    n = M.shape[0]
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i in range(n):
        for j in range(i + 1, n):
            if M[i, j] <= cutoff:
                parent[find(i)] = find(j)
    roots = {}
    labels = []
    for i in range(n):
        r = find(i)
        roots.setdefault(r, len(roots))
        labels.append(roots[r])
    return np.array(labels)


def analyze_target(pdbqt: Path, cutoff: float):
    poses, names, vina = parse_pdbqt_poses(pdbqt)
    if not poses:
        print(f"  {pdbqt.name}: sem poses; pulando")
        return None
    M = rmsd_matrix(poses)
    rmsd_to_best = M[0]                 # vs pose 1 (melhor afinidade)
    labels = single_linkage(M, cutoff)
    aff = np.array([v[0] for v in vina]) if vina else np.full(len(poses), np.nan)
    n_clusters = len(set(labels))
    target = pdbqt.stem.replace("result_", "")

    # figura
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    cmap = plt.get_cmap("tab10")
    idx = np.arange(1, len(poses) + 1)
    ax[0].bar(idx, aff, color=[cmap(l % 10) for l in labels])
    ax[0].set(xlabel="Pose (ranqueada)", ylabel="Afinidade Vina (kcal/mol)",
              title=f"{target}: afinidade por pose\n(cor = cluster, corte {cutoff} A)")
    ax[0].invert_yaxis()
    ax[1].plot(idx, rmsd_to_best, "o-", color="#d62728")
    ax[1].axhline(2.0, ls="--", color="gray", lw=1)
    ax[1].set(xlabel="Pose", ylabel="RMSD vs melhor pose (A)",
              title=f"{target}: RMSD heavy-atom vs pose 1")
    im = ax[2].imshow(M, cmap="viridis", vmin=0)
    ax[2].set(title=f"{target}: matriz RMSD pareada (A)",
              xlabel="Pose", ylabel="Pose")
    ax[2].set_xticks(range(len(poses))); ax[2].set_xticklabels(idx, fontsize=7)
    ax[2].set_yticks(range(len(poses))); ax[2].set_yticklabels(idx, fontsize=7)
    fig.colorbar(im, ax=ax[2], fraction=0.046)
    fig.tight_layout()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    out_png = FIG_DIR / f"docking_pose_rmsd_{target}.png"
    fig.savefig(out_png, dpi=130)
    plt.close(fig)

    summary = {
        "target": target, "n_poses": len(poses), "n_atoms": len(names),
        "cutoff_A": cutoff, "n_clusters": n_clusters,
        "cluster_labels": labels.tolist(),
        "affinity_kcalmol": [None if np.isnan(a) else round(float(a), 3) for a in aff],
        "rmsd_to_best_A": [round(float(x), 3) for x in rmsd_to_best],
        "best_affinity": None if np.isnan(aff[0]) else round(float(aff[0]), 3),
        "figure": out_png.name,
    }
    print(f"  {target}: {len(poses)} poses, {len(names)} heavy atoms, "
          f"{n_clusters} clusters @ {cutoff} A  -> {out_png.name}")
    return summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", default=None, help="ex.: TcTR (default: todos)")
    ap.add_argument("--cutoff", type=float, default=2.0, help="corte de RMSD (A)")
    args = ap.parse_args()

    files = ([DOCK_DIR / f"result_{args.target}.pdbqt"] if args.target
             else sorted(DOCK_DIR.glob("result_*.pdbqt")))
    out = []
    print(f"Analisando {len(files)} arquivo(s) de docking (corte {args.cutoff} A):")
    for f in files:
        if f.exists():
            s = analyze_target(f, args.cutoff)
            if s:
                out.append(s)
        else:
            print(f"  ausente: {f.name}")
    (FIG_DIR / "docking_pose_rmsd_summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nResumo: {FIG_DIR / 'docking_pose_rmsd_summary.json'}")


if __name__ == "__main__":
    main()
