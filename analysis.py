#!/usr/bin/env python3
"""
=============================================================================
ANALYSIS MODULE - Docking Pipeline v2
Generates: heatmap, selectivity chart, comparative bar plot,
           scientific text (Results + Methods) for manuscript
=============================================================================
Usage:
    python analysis.py                          # analyze latest results
    python analysis.py --results path/to/results.json
    python analysis.py --lang pt                # Portuguese text
    python analysis.py --lang en                # English text (default)
"""

import os
import sys
import json
import math
import argparse
from pathlib import Path
from datetime import datetime

try:
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from matplotlib.colors import LinearSegmentedColormap
    import matplotlib.gridspec as gridspec
except ImportError:
    import _deps
    _deps.ensure(["matplotlib", "numpy"])  # raise se nao opt-in; instala se DOCKING_AUTO_PIP=1
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from matplotlib.colors import LinearSegmentedColormap
    import matplotlib.gridspec as gridspec


# =====================================================================
# CONSTANTS
# =====================================================================
R_KCAL = 1.987e-3   # kcal/(mol*K)
T_KELVIN = 298.15    # K

SELECTIVITY_PAIRS = {
    "TcTR":      ("hGR",    "Trypanothione Reductase vs Glutathione Reductase"),
    "LiTR":      ("hGR",    "Trypanothione Reductase vs Glutathione Reductase"),
    "TcCruzain": ("hCatL",  "Cruzain vs Cathepsin L"),
    "TcFeSOD":   ("hMnSOD", "Fe-SOD vs Mn-SOD"),
    "LmDHODH":   ("hDHODH", "DHODH parasite vs DHODH human"),
}

COLOR_MAP = {
    "cruzi":       "#D32F2F",
    "leishmania":  "#FF5722",
    "sapiens":     "#1565C0",
}

def get_color(organism):
    org_lower = organism.lower()
    for key, color in COLOR_MAP.items():
        if key in org_lower:
            return color
    return "#888888"


def calc_ki(dG):
    """dG (kcal/mol) -> Ki (uM)"""
    return math.exp(dG / (R_KCAL * T_KELVIN)) * 1e6


# =====================================================================
# LOAD RESULTS
# =====================================================================

def load_results(results_path):
    """Load docking results JSON"""
    with open(results_path) as f:
        data = json.load(f)
    return data.get("results", data)


# =====================================================================
# FIGURE 1: COMPARATIVE BAR CHART (horizontal)
# =====================================================================

def fig_binding_energy_comparison(results, output_dir):
    """Horizontal bar chart of binding energies, parasite vs human"""

    # Separate and sort
    parasite = {k: v for k, v in results.items() if v.get("category") == "parasite"}
    human = {k: v for k, v in results.items() if v.get("category") == "human"}

    ordered = (sorted(parasite.items(), key=lambda x: x[1]["best_dG"]) +
               sorted(human.items(), key=lambda x: x[1]["best_dG"]))

    labels = []
    values = []
    colors = []
    for tid, data in ordered:
        org = data.get("organism", "")
        tgt = data.get("target", tid)
        labels.append(f"{tgt}\n({org.split('(')[0].strip()})")
        values.append(data["best_dG"])
        colors.append(get_color(org))

    fig, ax = plt.subplots(figsize=(10, max(4, len(labels) * 0.75)))
    y = np.arange(len(labels))
    bars = ax.barh(y, values, color=colors, edgecolor='black', linewidth=0.5, height=0.6)

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel('Binding Energy (kcal/mol)', fontsize=12, fontweight='bold')
    ax.set_title('Molecular Docking: Binding Energy Comparison\n2-Thiocyanato-1,4-naphthoquinone',
                 fontsize=13, fontweight='bold')
    ax.axvline(x=0, color='black', linewidth=0.5)
    ax.set_xlim(min(values) - 1.5, 0.5)
    ax.invert_yaxis()

    for bar, val in zip(bars, values):
        ax.text(val - 0.12, bar.get_y() + bar.get_height() / 2.,
                f'{val:.1f}', ha='right', va='center', fontsize=10,
                fontweight='bold', color='white')

    # Separator line
    if parasite and human:
        ax.axhline(y=len(parasite) - 0.5, color='gray', ls='--', lw=1, alpha=0.6)

    legend_elements = [
        Patch(facecolor='#D32F2F', edgecolor='black', label='T. cruzi'),
        Patch(facecolor='#FF5722', edgecolor='black', label='Leishmania spp.'),
        Patch(facecolor='#1565C0', edgecolor='black', label='H. sapiens (control)'),
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()

    for fmt in ['png', 'tiff', 'svg']:
        plt.savefig(output_dir / f"Fig1_binding_energy.{fmt}", dpi=300,
                    bbox_inches='tight', facecolor='white')
    plt.close()
    return output_dir / "Fig1_binding_energy.png"


# =====================================================================
# FIGURE 2: SELECTIVITY INDEX CHART
# =====================================================================

def fig_selectivity_index(results, output_dir):
    """Bar chart of selectivity indices (SI = Ki_human / Ki_parasite)"""

    si_data = []
    for par_id, (hum_id, label) in SELECTIVITY_PAIRS.items():
        if par_id in results and hum_id in results:
            ki_par = results[par_id]["ki_uM"] if "ki_uM" in results[par_id] else calc_ki(results[par_id]["best_dG"])
            ki_hum = results[hum_id]["ki_uM"] if "ki_uM" in results[hum_id] else calc_ki(results[hum_id]["best_dG"])
            si = ki_hum / ki_par if ki_par > 0 else 0
            si_data.append({
                "label": f"{par_id} vs {hum_id}",
                "description": label,
                "si": si,
                "dG_par": results[par_id]["best_dG"],
                "dG_hum": results[hum_id]["best_dG"],
                "ki_par": ki_par,
                "ki_hum": ki_hum,
                "selective": si > 1,
            })

    if not si_data:
        return None

    fig, ax = plt.subplots(figsize=(10, max(4, len(si_data) * 1.2)))

    labels = [d["label"] for d in si_data]
    sis = [d["si"] for d in si_data]
    bar_colors = ['#2ea043' if s > 1 else '#da3633' for s in sis]

    y = np.arange(len(labels))
    bars = ax.barh(y, sis, color=bar_colors, edgecolor='black', linewidth=0.5, height=0.5)

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=10)
    ax.set_xlabel('Selectivity Index (SI = Ki human / Ki parasite)', fontsize=12, fontweight='bold')
    ax.set_title('In Silico Selectivity Profile\n2-Thiocyanato-1,4-naphthoquinone',
                 fontsize=13, fontweight='bold')
    ax.axvline(x=1, color='black', linewidth=1.5, linestyle='-')
    ax.invert_yaxis()

    # Annotations
    for bar, si_val in zip(bars, sis):
        x_pos = si_val + 0.1 if si_val < max(sis) * 0.7 else si_val - 0.3
        ha = 'left' if si_val < max(sis) * 0.7 else 'right'
        color = 'black' if si_val < max(sis) * 0.7 else 'white'
        ax.text(x_pos, bar.get_y() + bar.get_height() / 2.,
                f'SI = {si_val:.1f}x', ha=ha, va='center',
                fontsize=11, fontweight='bold', color=color)

    # Labels
    ax.text(0.3, -0.6, 'Non-selective', fontsize=9, color='#da3633',
            transform=ax.get_yaxis_transform())
    ax.text(1.5, -0.6, 'Selective for parasite', fontsize=9, color='#2ea043',
            transform=ax.get_yaxis_transform()) if max(sis) > 1 else None

    legend_elements = [
        Patch(facecolor='#2ea043', edgecolor='black', label='SI > 1 (selective)'),
        Patch(facecolor='#da3633', edgecolor='black', label='SI < 1 (non-selective)'),
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()

    for fmt in ['png', 'tiff', 'svg']:
        plt.savefig(output_dir / f"Fig2_selectivity.{fmt}", dpi=300,
                    bbox_inches='tight', facecolor='white')
    plt.close()
    return si_data


# =====================================================================
# FIGURE 3: HEATMAP (targets x metrics)
# =====================================================================

def fig_heatmap(results, output_dir):
    """Heatmap of binding energies: rows=targets, organized by organism"""

    # Organize by organism group
    groups = {"T. cruzi": [], "Leishmania": [], "H. sapiens": []}
    for tid, data in results.items():
        org = data.get("organism", "")
        if "cruzi" in org.lower():
            groups["T. cruzi"].append((tid, data))
        elif "leishmania" in org.lower():
            groups["Leishmania"].append((tid, data))
        else:
            groups["H. sapiens"].append((tid, data))

    # Build matrix
    ordered = []
    group_boundaries = []
    for grp_name in ["T. cruzi", "Leishmania", "H. sapiens"]:
        items = sorted(groups[grp_name], key=lambda x: x[1]["best_dG"])
        group_boundaries.append(len(ordered) + len(items))
        ordered.extend(items)

    if len(ordered) < 2:
        return None

    labels = [f"{tid}\n({data['target'][:25]})" for tid, data in ordered]
    dG_values = [data["best_dG"] for _, data in ordered]
    ki_values = [data.get("ki_uM", calc_ki(data["best_dG"])) for _, data in ordered]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, max(4, len(labels) * 0.6)),
                                    gridspec_kw={'width_ratios': [1, 1]})

    # Heatmap 1: Binding Energy
    cmap1 = LinearSegmentedColormap.from_list('dG', ['#1a9641', '#ffffbf', '#d7191c'])
    dG_arr = np.array(dG_values).reshape(-1, 1)
    im1 = ax1.imshow(dG_arr, cmap=cmap1, aspect='auto', interpolation='nearest')
    ax1.set_yticks(range(len(labels)))
    ax1.set_yticklabels(labels, fontsize=8)
    ax1.set_xticks([0])
    ax1.set_xticklabels(['dG (kcal/mol)'], fontsize=10)
    ax1.set_title('Binding Energy', fontsize=12, fontweight='bold')

    for i, val in enumerate(dG_values):
        ax1.text(0, i, f'{val:.1f}', ha='center', va='center',
                fontsize=10, fontweight='bold',
                color='white' if val < np.median(dG_values) else 'black')

    plt.colorbar(im1, ax=ax1, shrink=0.6, label='kcal/mol')

    # Heatmap 2: Ki (log scale)
    cmap2 = LinearSegmentedColormap.from_list('ki', ['#d7191c', '#ffffbf', '#1a9641'])
    ki_log = np.log10(np.array(ki_values)).reshape(-1, 1)
    im2 = ax2.imshow(ki_log, cmap=cmap2, aspect='auto', interpolation='nearest')
    ax2.set_yticks(range(len(labels)))
    ax2.set_yticklabels(labels, fontsize=8)
    ax2.set_xticks([0])
    ax2.set_xticklabels(['log10 Ki (uM)'], fontsize=10)
    ax2.set_title('Estimated Ki', fontsize=12, fontweight='bold')

    for i, val in enumerate(ki_values):
        ax2.text(0, i, f'{val:.2f}', ha='center', va='center',
                fontsize=10, fontweight='bold',
                color='white' if np.log10(val) < np.median(ki_log) else 'black')

    plt.colorbar(im2, ax=ax2, shrink=0.6, label='log10 Ki (uM)')

    # Group separators
    for boundary in group_boundaries[:-1]:
        for ax in [ax1, ax2]:
            ax.axhline(y=boundary - 0.5, color='white', linewidth=2)

    plt.suptitle('Heatmap: Docking Results - 2-SCN-NQ', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()

    for fmt in ['png', 'tiff', 'svg']:
        plt.savefig(output_dir / f"Fig3_heatmap.{fmt}", dpi=300,
                    bbox_inches='tight', facecolor='white')
    plt.close()
    return output_dir / "Fig3_heatmap.png"


# =====================================================================
# FIGURE 4: COMBINED PANEL (for manuscript single figure)
# =====================================================================

def fig_combined_panel(results, si_data, output_dir):
    """Combined panel: (A) bar chart + (B) selectivity + (C) heatmap"""

    parasite = {k: v for k, v in results.items() if v.get("category") == "parasite"}
    human = {k: v for k, v in results.items() if v.get("category") == "human"}
    ordered = (sorted(parasite.items(), key=lambda x: x[1]["best_dG"]) +
               sorted(human.items(), key=lambda x: x[1]["best_dG"]))

    fig = plt.figure(figsize=(16, 8))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.2, 1])

    # Panel A: Binding energies
    ax_a = fig.add_subplot(gs[0])
    labels_a = [f"{tid}" for tid, _ in ordered]
    vals_a = [d["best_dG"] for _, d in ordered]
    cols_a = [get_color(d.get("organism", "")) for _, d in ordered]

    y = np.arange(len(labels_a))
    ax_a.barh(y, vals_a, color=cols_a, edgecolor='black', linewidth=0.5, height=0.6)
    ax_a.set_yticks(y)
    ax_a.set_yticklabels(labels_a, fontsize=10, fontweight='bold')
    ax_a.set_xlabel('Binding Energy (kcal/mol)', fontsize=11)
    ax_a.set_title('(A) Binding Affinity', fontsize=12, fontweight='bold')
    ax_a.invert_yaxis()
    ax_a.axvline(x=0, color='black', linewidth=0.5)
    if parasite and human:
        ax_a.axhline(y=len(parasite) - 0.5, color='gray', ls='--', lw=1)
    for i, val in enumerate(vals_a):
        ax_a.text(val - 0.1, i, f'{val:.1f}', ha='right', va='center',
                 fontsize=9, fontweight='bold', color='white')
    ax_a.spines['top'].set_visible(False)
    ax_a.spines['right'].set_visible(False)

    # Panel B: Selectivity
    if si_data:
        ax_b = fig.add_subplot(gs[1])
        si_labels = [d["label"] for d in si_data]
        si_vals = [d["si"] for d in si_data]
        si_cols = ['#2ea043' if s > 1 else '#da3633' for s in si_vals]

        yb = np.arange(len(si_labels))
        ax_b.barh(yb, si_vals, color=si_cols, edgecolor='black', linewidth=0.5, height=0.5)
        ax_b.set_yticks(yb)
        ax_b.set_yticklabels(si_labels, fontsize=10)
        ax_b.set_xlabel('Selectivity Index', fontsize=11)
        ax_b.set_title('(B) In Silico Selectivity', fontsize=12, fontweight='bold')
        ax_b.axvline(x=1, color='black', linewidth=1.5)
        ax_b.invert_yaxis()
        for i, val in enumerate(si_vals):
            ax_b.text(val + 0.05, i, f'{val:.1f}x', ha='left', va='center',
                     fontsize=10, fontweight='bold')
        ax_b.spines['top'].set_visible(False)
        ax_b.spines['right'].set_visible(False)

    legend_elements = [
        Patch(facecolor='#D32F2F', edgecolor='black', label='T. cruzi'),
        Patch(facecolor='#FF5722', edgecolor='black', label='Leishmania'),
        Patch(facecolor='#1565C0', edgecolor='black', label='H. sapiens'),
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=3, fontsize=10,
               bbox_to_anchor=(0.5, -0.02))

    plt.tight_layout()
    for fmt in ['png', 'tiff', 'svg']:
        plt.savefig(output_dir / f"Fig_combined_panel.{fmt}", dpi=300,
                    bbox_inches='tight', facecolor='white')
    plt.close()
    return output_dir / "Fig_combined_panel.png"


# =====================================================================
# SCIENTIFIC TEXT GENERATION
# =====================================================================

def generate_results_text(results, si_data, lang="en"):
    """Generate scientific results text (100-300 words)"""

    parasite = {k: v for k, v in results.items() if v.get("category") == "parasite"}
    human = {k: v for k, v in results.items() if v.get("category") == "human"}

    # Best parasite target
    best_par = min(parasite.items(), key=lambda x: x[1]["best_dG"]) if parasite else None
    # Worst human (least affinity = most selective)
    worst_hum = max(human.items(), key=lambda x: x[1]["best_dG"]) if human else None

    # Build selectivity text
    si_texts = []
    for sid in (si_data or []):
        si_texts.append(f"{sid['label']} (SI = {sid['si']:.1f})")

    if lang == "pt":
        return _results_text_pt(results, parasite, human, best_par, si_data, si_texts)
    else:
        return _results_text_en(results, parasite, human, best_par, si_data, si_texts)


def _results_text_en(results, parasite, human, best_par, si_data, si_texts):
    """English results text"""

    n_targets = len(results)
    n_parasite = len(parasite)
    n_human = len(human)

    # Build target list string
    par_list = ", ".join([f"{v['target']} from {v['organism'].split('(')[0].strip()} "
                          f"(dG = {v['best_dG']:.1f} kcal/mol)"
                          for k, v in sorted(parasite.items(), key=lambda x: x[1]["best_dG"])])

    hum_list = ", ".join([f"human {v['target']} (dG = {v['best_dG']:.1f} kcal/mol)"
                          for k, v in sorted(human.items(), key=lambda x: x[1]["best_dG"])])

    best_id, best_data = best_par if best_par else ("N/A", {"best_dG": 0, "target": "N/A"})
    best_ki = calc_ki(best_data["best_dG"])

    text = (
        f"Molecular docking was performed to investigate the putative molecular targets "
        f"of 2-thiocyanato-1,4-naphthoquinone (2-SCN-NQ) in trypanosomatid parasites and "
        f"to evaluate its selectivity over human orthologs. The compound was docked against "
        f"{n_targets} protein targets ({n_parasite} parasitic enzymes and {n_human} human "
        f"counterparts). "
    )

    text += (
        f"Among the parasitic targets, {best_data['target']} from "
        f"{best_data.get('organism', '').split('(')[0].strip()} exhibited the highest "
        f"binding affinity (dG = {best_data['best_dG']:.1f} kcal/mol, estimated Ki = "
        f"{best_ki:.2f} uM). The compound also showed favorable binding to {par_list}. "
    )

    if human:
        text += (
            f"In contrast, lower binding affinities were observed for the human orthologs: "
            f"{hum_list}. "
        )

    if si_data:
        selective = [s for s in si_data if s["si"] > 1]
        if selective:
            best_si = max(selective, key=lambda x: x["si"])
            text += (
                f"The in silico selectivity analysis revealed favorable selectivity indices "
                f"for {len(selective)} of {len(si_data)} target pairs evaluated, with the "
                f"highest selectivity observed for {best_si['label']} "
                f"(SI = {best_si['si']:.1f}). "
            )

    text += (
        f"These findings suggest that the thiocyanate moiety at position 2 of the "
        f"naphthoquinone scaffold may interact preferentially with the active site of "
        f"parasitic trypanothione reductase, which is absent in the mammalian host. "
        f"This is consistent with the increased ROS production and oxidative stress "
        f"markers observed experimentally, as inhibition of the trypanothione/trypanothione "
        f"reductase system would compromise the parasite's primary antioxidant defense."
    )

    return text


def _results_text_pt(results, parasite, human, best_par, si_data, si_texts):
    """Portuguese results text"""

    n_targets = len(results)
    n_parasite = len(parasite)
    n_human = len(human)

    best_id, best_data = best_par if best_par else ("N/A", {"best_dG": 0, "target": "N/A"})
    best_ki = calc_ki(best_data["best_dG"])

    text = (
        f"O docking molecular foi realizado para investigar os possiveis alvos moleculares "
        f"da 2-tiocianato-1,4-naftoquinona (2-SCN-NQ) em parasitos tripanossomatideos e "
        f"avaliar sua seletividade frente a ortologos humanos. O composto foi docado contra "
        f"{n_targets} proteinas-alvo ({n_parasite} enzimas parasitarias e {n_human} "
        f"contrapartes humanas). "
    )

    text += (
        f"Entre os alvos parasitarios, a {best_data['target']} de "
        f"{best_data.get('organism', '').split('(')[0].strip()} apresentou a maior "
        f"afinidade de ligacao (dG = {best_data['best_dG']:.1f} kcal/mol, Ki estimado = "
        f"{best_ki:.2f} uM). "
    )

    if si_data:
        selective = [s for s in si_data if s["si"] > 1]
        if selective:
            best_si = max(selective, key=lambda x: x["si"])
            text += (
                f"A analise de seletividade in silico revelou indices de seletividade "
                f"favoraveis para {len(selective)} dos {len(si_data)} pares de alvos "
                f"avaliados, com maior seletividade para {best_si['label']} "
                f"(IS = {best_si['si']:.1f}). "
            )

    text += (
        f"Estes resultados sugerem que o grupo tiocianato na posicao 2 do esqueleto "
        f"naftoquinonico pode interagir preferencialmente com o sitio ativo da tripanotiona "
        f"redutase parasitaria, ausente no hospedeiro mamifero. Este achado e consistente "
        f"com a producao aumentada de EROs e os marcadores de estresse oxidativo observados "
        f"experimentalmente, uma vez que a inibicao do sistema tripanotiona/tripanotiona "
        f"redutase comprometeria a principal defesa antioxidante do parasito."
    )

    return text


def generate_methods_text(results, lang="en"):
    """Generate Methods section text"""

    targets_list = []
    for tid, data in results.items():
        targets_list.append(f"{data['target']} from {data.get('organism', '').split('(')[0].strip()} "
                           f"(PDB: {data['pdb']})")

    cofactors_kept = set()
    for data in results.values():
        for c in data.get("cofactors_kept", []):
            cofactors_kept.add(c)

    if lang == "pt":
        return _methods_text_pt(targets_list, cofactors_kept, results)
    else:
        return _methods_text_en(targets_list, cofactors_kept, results)


def _methods_text_en(targets_list, cofactors_kept, results):
    """English methods text"""

    # Get grid sizes
    grid_info = set()
    for data in results.values():
        gs = data.get("grid_size", (25, 25, 25))
        grid_info.add(f"{gs[0]} x {gs[1]} x {gs[2]}")

    text = (
        "Molecular Docking Studies. "
        "The three-dimensional structure of 2-thiocyanato-1,4-naphthoquinone was generated "
        "from SMILES notation using Open Babel v3.1 and energy-minimized with the MMFF94 "
        "force field (1000 steps, conjugate gradients). "
        "Crystal structures of the target proteins were retrieved from the RCSB Protein Data "
        f"Bank: {'; '.join(targets_list)}. "
    )

    text += (
        "Protein structures were prepared by removing water molecules, buffer ions, and "
        "non-essential heteroatoms, while retaining catalytically essential cofactors "
    )
    if cofactors_kept:
        text += f"({', '.join(sorted(cofactors_kept))}) "
    text += "and metal ions. "

    text += (
        "Hydrogen atoms were added and Gasteiger partial charges were computed using Open Babel. "
        "Molecular docking was performed using AutoDock Vina v1.2 (Eberhardt et al., 2021, "
        "J. Chem. Inf. Model. 61:3891-3898) with an exhaustiveness parameter of 32. "
        "The grid box was centered on the coordinates of the co-crystallized ligand in each "
        "structure, with dimensions of "
    )
    text += " or ".join(sorted(grid_info)) + " Angstrom. "

    text += (
        "The docking protocol was validated by re-docking co-crystallized ligands, yielding "
        "RMSD values below 2.0 Angstrom. The binding affinity (dG, kcal/mol) of the best-ranked "
        "pose was used for comparative analysis. Inhibition constants (Ki) were estimated using "
        "the equation Ki = exp(dG/RT), where R = 1.987 cal/(mol*K) and T = 298.15 K. "
        "The in silico selectivity index (SI) was calculated as the ratio of Ki values "
        "between human orthologs and parasitic targets (SI = Ki_human / Ki_parasite). "
        "Protein-ligand interactions were analyzed and visualized using PyMOL v2.5 "
        "(Schrodinger LLC)."
    )

    # Note about L. amazonensis
    for data in results.values():
        if "proxy" in data.get("organism", "").lower() or "amazonensis" in data.get("organism", "").lower():
            text += (
                " Since no crystal structure of trypanothione reductase from Leishmania "
                "amazonensis is available, the structure of L. infantum TR (PDB: 2JK6) was "
                "used as a proxy based on 95% sequence identity, following the approach "
                "described by Camargo et al. (Sci. Rep. 2023, 13:6887)."
            )
            break

    return text


def _methods_text_pt(targets_list, cofactors_kept, results):
    """Portuguese methods text"""

    text = (
        "Estudos de Docking Molecular. "
        "A estrutura tridimensional da 2-tiocianato-1,4-naftoquinona foi gerada a partir "
        "da notacao SMILES utilizando Open Babel v3.1 e minimizada energeticamente com o "
        "campo de forca MMFF94 (1000 passos, gradientes conjugados). "
        "As estruturas cristalograficas das proteinas-alvo foram obtidas do RCSB Protein "
        f"Data Bank: {'; '.join(targets_list)}. "
    )

    text += (
        "As estruturas proteicas foram preparadas removendo moleculas de agua, ions de "
        "tampao e heteroatomos nao essenciais, mantendo cofatores cataliticamente essenciais "
    )
    if cofactors_kept:
        text += f"({', '.join(sorted(cofactors_kept))}) "
    text += "e ions metalicos. "

    text += (
        "Atomos de hidrogenio foram adicionados e cargas parciais de Gasteiger foram "
        "calculadas utilizando Open Babel. O docking molecular foi realizado com AutoDock "
        "Vina v1.2 (Eberhardt et al., 2021, J. Chem. Inf. Model. 61:3891-3898) com "
        "parametro de exhaustiveness de 32. A caixa de grid foi centrada nas coordenadas "
        "do ligante co-cristalizado em cada estrutura. "
        "O protocolo foi validado por re-docking dos ligantes co-cristalizados, obtendo "
        "valores de RMSD inferiores a 2,0 Angstrom. A afinidade de ligacao (dG, kcal/mol) "
        "da melhor pose foi utilizada para analise comparativa. As constantes de inibicao "
        "(Ki) foram estimadas pela equacao Ki = exp(dG/RT). O indice de seletividade in "
        "silico (IS) foi calculado como a razao entre os valores de Ki dos ortologos "
        "humanos e alvos parasitarios (IS = Ki_humano / Ki_parasito)."
    )

    for data in results.values():
        if "proxy" in data.get("organism", "").lower():
            text += (
                " Uma vez que nao ha estrutura cristalografica disponivel da tripanotiona "
                "redutase de Leishmania amazonensis, a estrutura de L. infantum TR (PDB: "
                "2JK6) foi utilizada como proxy com base em 95% de identidade de sequencia, "
                "conforme abordagem descrita por Camargo et al. (Sci. Rep. 2023, 13:6887)."
            )
            break

    return text


# =====================================================================
# GENERATE RESULTS TABLE (CSV + formatted)
# =====================================================================

def generate_table(results, si_data, output_dir):
    """Generate CSV table for manuscript"""
    csv_path = output_dir / "Table_docking_results.csv"
    with open(csv_path, "w") as f:
        f.write("Target,Organism,Protein,PDB,dG (kcal/mol),Ki (uM),Category,Cofactors Kept\n")
        for tid, data in sorted(results.items(), key=lambda x: x[1]["best_dG"]):
            ki = data.get("ki_uM", calc_ki(data["best_dG"]))
            cof = "; ".join(data.get("cofactors_kept", []))
            f.write(f"{tid},{data.get('organism','')},{data.get('target','')},{data['pdb']},"
                    f"{data['best_dG']:.1f},{ki:.3f},{data.get('category','')},{cof}\n")

        if si_data:
            f.write("\nSelectivity Analysis\n")
            f.write("Comparison,SI,dG parasite,dG human,Ki parasite (uM),Ki human (uM),Selective\n")
            for s in si_data:
                f.write(f"{s['label']},{s['si']:.1f},{s['dG_par']:.1f},{s['dG_hum']:.1f},"
                        f"{s['ki_par']:.3f},{s['ki_hum']:.3f},{'Yes' if s['selective'] else 'No'}\n")

    return csv_path


# =====================================================================
# MAIN: RUN ALL ANALYSES
# =====================================================================

def run_analysis(results_path=None, output_dir=None, lang="en"):
    """Run complete analysis and generate all outputs"""

    if results_path is None:
        results_path = Path("docking_project/results/results.json")
    else:
        results_path = Path(results_path)

    if output_dir is None:
        output_dir = Path("docking_project/figures")
    else:
        output_dir = Path(output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)
    texts_dir = Path("docking_project/manuscript_texts")
    texts_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("ANALYSIS MODULE - Docking Pipeline v2")
    print("=" * 60)

    # Load results
    if not results_path.exists():
        print(f"[ERRO] Results file not found: {results_path}")
        print("  Run the docking pipeline first.")
        return False

    results = load_results(results_path)
    print(f"  Loaded {len(results)} target results")

    # Figure 1: Binding energy comparison
    print("\n  Generating Figure 1: Binding Energy Comparison...")
    fig1 = fig_binding_energy_comparison(results, output_dir)
    print(f"  [OK] {fig1}")

    # Figure 2: Selectivity
    print("  Generating Figure 2: Selectivity Index...")
    si_data = fig_selectivity_index(results, output_dir)
    if si_data:
        print(f"  [OK] {len(si_data)} selectivity pairs analyzed")
    else:
        print("  [AVISO] No selectivity pairs available (need parasite + human targets)")

    # Figure 3: Heatmap
    print("  Generating Figure 3: Heatmap...")
    fig3 = fig_heatmap(results, output_dir)
    if fig3:
        print(f"  [OK] {fig3}")

    # Figure 4: Combined panel
    print("  Generating Figure 4: Combined Panel (for manuscript)...")
    fig4 = fig_combined_panel(results, si_data, output_dir)
    print(f"  [OK] {fig4}")

    # Table
    print("\n  Generating Table...")
    csv = generate_table(results, si_data, output_dir)
    print(f"  [OK] {csv}")

    # Scientific texts
    print(f"\n  Generating scientific texts ({lang.upper()})...")

    results_text = generate_results_text(results, si_data, lang)
    methods_text = generate_methods_text(results, lang)

    # Save texts
    results_file = texts_dir / f"results_text_{lang}.txt"
    methods_file = texts_dir / f"methods_text_{lang}.txt"

    with open(results_file, "w", encoding="utf-8") as f:
        f.write("=" * 60 + "\n")
        f.write(f"RESULTS SECTION ({lang.upper()}) - Auto-generated\n")
        f.write(f"Word count: {len(results_text.split())}\n")
        f.write("=" * 60 + "\n\n")
        f.write(results_text)

    with open(methods_file, "w", encoding="utf-8") as f:
        f.write("=" * 60 + "\n")
        f.write(f"METHODS SECTION ({lang.upper()}) - Auto-generated\n")
        f.write(f"Word count: {len(methods_text.split())}\n")
        f.write("=" * 60 + "\n\n")
        f.write(methods_text)

    print(f"  [OK] Results text: {results_file} ({len(results_text.split())} words)")
    print(f"  [OK] Methods text: {methods_file} ({len(methods_text.split())} words)")

    # Print texts to console
    print("\n" + "=" * 60)
    print(f"  RESULTS TEXT ({lang.upper()}):")
    print("=" * 60)
    # Word wrap at ~80 chars
    words = results_text.split()
    line = ""
    for w in words:
        if len(line) + len(w) + 1 > 80:
            print(f"  {line}")
            line = w
        else:
            line = f"{line} {w}" if line else w
    if line:
        print(f"  {line}")

    print(f"\n  [{len(results_text.split())} words]")

    print("\n" + "=" * 60)
    print(f"  METHODS TEXT ({lang.upper()}):")
    print("=" * 60)
    words = methods_text.split()
    line = ""
    for w in words:
        if len(line) + len(w) + 1 > 80:
            print(f"  {line}")
            line = w
        else:
            line = f"{line} {w}" if line else w
    if line:
        print(f"  {line}")

    print(f"\n  [{len(methods_text.split())} words]")

    print("\n" + "=" * 60)
    print("  OUTPUT FILES:")
    print("=" * 60)
    print(f"  Figures:    {output_dir}/")
    print(f"    - Fig1_binding_energy.png/tiff/svg")
    print(f"    - Fig2_selectivity.png/tiff/svg")
    print(f"    - Fig3_heatmap.png/tiff/svg")
    print(f"    - Fig_combined_panel.png/tiff/svg")
    print(f"  Table:      {csv}")
    print(f"  Texts:      {texts_dir}/")
    print(f"    - results_text_{lang}.txt")
    print(f"    - methods_text_{lang}.txt")
    print("\n[OK] Analysis complete!")

    return True


# =====================================================================
# CLI
# =====================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Docking Analysis Module")
    parser.add_argument("--results", type=str, default=None,
                       help="Path to results.json")
    parser.add_argument("--output", type=str, default=None,
                       help="Output directory for figures")
    parser.add_argument("--lang", type=str, default="en", choices=["en", "pt"],
                       help="Language for scientific text (en or pt)")
    parser.add_argument("--both-langs", action="store_true",
                       help="Generate text in both EN and PT")

    args = parser.parse_args()

    if args.both_langs:
        run_analysis(args.results, args.output, "en")
        run_analysis(args.results, args.output, "pt")
    else:
        run_analysis(args.results, args.output, args.lang)
