#!/usr/bin/env python3
"""
=============================================================================
ADMET IN SILICO MODULE - Pharmacokinetic profiling
Calculates drug-likeness, absorption, distribution, metabolism and
toxicity properties for 2-SCN-NQ using RDKit + classical QSAR models.
=============================================================================
Usage:
    python admet_module.py
    python admet_module.py --smiles "O=C1C=C(SC#N)C(=O)c2ccccc21"
    python admet_module.py --output figures/
"""

import sys
import json
import math
import argparse
import threading
from pathlib import Path
from datetime import datetime

# ─── Dependencies ────────────────────────────────────────────────────────────
try:
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.patches import FancyBboxPatch
    from matplotlib.gridspec import GridSpec
except ImportError:
    import _deps
    _deps.ensure(["numpy", "matplotlib"])  # raise se nao opt-in; instala se DOCKING_AUTO_PIP=1
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.patches import FancyBboxPatch
    from matplotlib.gridspec import GridSpec

RDKIT_OK = False
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem, Draw
    from rdkit.Chem.Draw import rdMolDraw2D
    RDKIT_OK = True
except ImportError:
    pass

# ─── Compound defaults ───────────────────────────────────────────────────────
SMILES_DEFAULT  = "O=C1C=C(SC#N)C(=O)c2ccccc21"
COMPOUND_NAME   = "2-Thiocyanato-1,4-naphthoquinone"
COMPOUND_ABBR   = "2-SCN-NQ"
# Serializa a secao critica que muta/usa/restaura os globals COMPOUND_NAME/ABBR
# (run_admet) para evitar race condition sob Flask multithread.
_LOCK = threading.Lock()
COMPOUND_FORMULA = "C₁₁H₅NO₂S"

# ─── Lipinski / ADMET thresholds ─────────────────────────────────────────────
LIP = dict(MW_max=500, logP_max=5, HBD_max=5, HBA_max=10)
VEBER = dict(TPSA_max=140, RB_max=10)

# =============================================================================
# CORE CALCULATION ENGINE
# =============================================================================

def calculate_admet(smiles: str = SMILES_DEFAULT) -> dict:
    """Full ADMET profile. Returns a rich dict of all parameters."""

    if not RDKIT_OK:
        return {"error": "RDKit not installed. Run: pip install rdkit"}

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"error": f"Invalid SMILES: {smiles}"}

    # ── Physicochemical ────────────────────────────────────────────────────────
    mw          = Descriptors.MolWt(mol)
    exact_mw    = Descriptors.ExactMolWt(mol)
    logp        = Descriptors.MolLogP(mol)
    hbd         = rdMolDescriptors.CalcNumHBD(mol)
    hba         = rdMolDescriptors.CalcNumHBA(mol)
    tpsa        = Descriptors.TPSA(mol)
    rotbonds    = rdMolDescriptors.CalcNumRotatableBonds(mol)
    n_rings     = rdMolDescriptors.CalcNumRings(mol)
    arom_rings  = rdMolDescriptors.CalcNumAromaticRings(mol)
    n_heavy     = mol.GetNumHeavyAtoms()
    mol_ref     = Descriptors.MolMR(mol)
    formula     = rdMolDescriptors.CalcMolFormula(mol)
    fsp3        = rdMolDescriptors.CalcFractionCSP3(mol)
    stereo      = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    complexity  = Descriptors.BertzCT(mol)

    # ── Lipinski Rule of 5 ────────────────────────────────────────────────────
    lip_viol = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    lip_pass = lip_viol <= 1

    # ── Veber oral bioavailability rules ─────────────────────────────────────
    veber_pass = tpsa <= 140 and rotbonds <= 10

    # ── Ghose filter (drug-like) ─────────────────────────────────────────────
    ghose_pass = (160 <= mw <= 480) and (-0.4 <= logp <= 5.6) and (40 <= mol_ref <= 130) and (20 <= n_heavy <= 70)

    # ── Egan filter (oral bioavailability) ───────────────────────────────────
    egan_pass = (logp <= 5.88) and (tpsa <= 131.6)

    # ── Muegge filter (drug-like space) ──────────────────────────────────────
    muegge_pass = (200 <= mw <= 600) and (-2 <= logp <= 5) and (tpsa <= 150) and (n_rings <= 7) and (n_rings >= 1) and (rotbonds <= 15) and (hbd <= 5) and (hba <= 10)

    # ── Bioavailability score (0–1, based on Ertl & Murruzzu 2020) ────────────
    bio_score = _bioavailability_score(mw, logp, hbd, hba, tpsa, rotbonds)

    # ── Caco-2 permeability (Hou et al. 2004 model) ──────────────────────────
    # logPapp (x10⁻⁶ cm/s) = -0.0197·TPSA - 0.4496·nHBD + 0.5234·clogP - 2.7101
    log_papp = -0.0197 * tpsa - 0.4496 * hbd + 0.5234 * logp - 2.7101
    papp     = 10 ** log_papp
    caco2_class = "High" if log_papp > -5.15 else ("Medium" if log_papp > -5.85 else "Low")

    # ── PAMPA permeability (Bermejo et al. 2004 model) ────────────────────────
    log_pampa = -0.028 * tpsa + 0.397 * logp - 3.27
    pampa_class = "Good" if log_pampa > -6 else ("Moderate" if log_pampa > -7 else "Poor")

    # ── Intestinal absorption (%HIA, Zhao et al.) ────────────────────────────
    # %HIA ≈ -0.7416·TPSA + 97.22 (simplified for TPSA-based model)
    hia = max(0, min(100, -0.7416 * tpsa + 97.22))

    # ── BBB penetration (Clark 1999) ──────────────────────────────────────────
    log_bb = 0.0148 * logp - 0.0187 * tpsa + 0.16
    bbb_pass  = log_bb > -1.0

    # ── Water solubility ESOL (Delaney 2004) ─────────────────────────────────
    log_s = 0.16 - 0.63 * logp - 0.0062 * mw + 0.066 * rotbonds - 0.74 * arom_rings
    sol_mg_ml = (10 ** log_s) * mw  # mg/mL
    sol_class  = ("Very soluble" if log_s > -1 else
                  "Soluble"          if log_s > -2 else
                  "Moderately soluble" if log_s > -4 else
                  "Poorly soluble"   if log_s > -6 else
                  "Practically insoluble")

    # ── Plasma protein binding (Votano 2004 simple estimate) ─────────────────
    ppb_est = min(100, 82 + 5 * logp - 0.3 * tpsa)  # % bound
    ppb_class = "High (>90%)" if ppb_est > 90 else ("Moderate (70-90%)" if ppb_est > 70 else "Low (<70%)")

    # ── Volume of distribution estimate ──────────────────────────────────────
    vd_est = 10 ** (0.5 * logp - 0.01 * tpsa + 0.5)   # L/kg approx
    vd_class = "Low (<1 L/kg)" if vd_est < 1 else ("Moderate (1–3 L/kg)" if vd_est < 3 else "High (>3 L/kg)")

    # ── Half-life estimate (simplified) ──────────────────────────────────────
    # Based on Vd and typical clearance for similar compounds
    cl_est  = 15 * (mw / 400)  # mL/min/kg rough estimate
    t_half  = round((0.693 * vd_est * 1000) / (cl_est * 16.7), 1)  # hours
    t_half  = max(0.5, min(48, t_half))

    # ── CYP450 structural alerts ──────────────────────────────────────────────
    cyp_alerts = _check_cyp_alerts(mol)

    # ── hERG toxicity risk ────────────────────────────────────────────────────
    herg_risk = _herg_risk(logp, n_rings, hba)

    # ── PAINS (Pan-Assay Interference) filter ────────────────────────────────
    pains_alerts = _pains_filter(mol)

    # ── QED drug-likeness score (Bickerton et al. 2012) ──────────────────────
    qed = _calc_qed(mw, logp, hbd, hba, tpsa, rotbonds, arom_rings)

    # ── Synthetic Accessibility Score (Ertl & Schuffenhauer 2009) ────────────
    try:
        from rdkit.Chem import RDConfig
        import os
        sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
        from sascorer import calculateScore as sas_score
        sas = round(sas_score(mol), 2)
    except Exception:
        sas = round(1 + (complexity / 5000) * 9, 2)  # fallback estimate
    sas = max(1.0, min(10.0, sas))
    sas_class = ("Easy (SAS 1-3)" if sas <= 3 else
                 "Moderate (SAS 3-6)" if sas <= 6 else
                 "Difficult (SAS 6-10)")

    # ── AMES mutagenicity structural alert ────────────────────────────────────
    ames_alert = _ames_alert(mol)

    # ── Lead-like / Fragment-like ─────────────────────────────────────────────
    lead_like = (200 <= mw <= 350) and (-1 <= logp <= 3.5)
    frag_like = mw < 300 and logp < 3 and hbd <= 3 and hba <= 3

    return {
        "metadata": {
            "compound": COMPOUND_NAME,
            "abbreviation": COMPOUND_ABBR,
            "smiles": smiles,
            "formula": formula,
            "timestamp": datetime.now().isoformat(),
        },
        "physicochemical": {
            "MW":             round(mw,       2),
            "exact_MW":       round(exact_mw, 4),
            "logP":           round(logp,     2),
            "HBD":            hbd,
            "HBA":            hba,
            "TPSA":           round(tpsa,     1),
            "rotatable_bonds": rotbonds,
            "rings":          n_rings,
            "aromatic_rings": arom_rings,
            "heavy_atoms":    n_heavy,
            "Fsp3":           round(fsp3, 3),
            "molar_refraction": round(mol_ref, 2),
            "complexity":     round(complexity, 1),
            "stereocenters":  stereo,
        },
        "drug_likeness": {
            "lipinski": {
                "pass":       lip_pass,
                "violations": lip_viol,
                "MW_ok":      mw  <= 500,
                "logP_ok":    logp <= 5,
                "HBD_ok":     hbd  <= 5,
                "HBA_ok":     hba  <= 10,
            },
            "veber": {
                "pass":          veber_pass,
                "TPSA_ok":       tpsa     <= 140,
                "rotbonds_ok":   rotbonds <= 10,
            },
            "ghose":    {"pass": ghose_pass},
            "egan":     {"pass": egan_pass},
            "muegge":   {"pass": muegge_pass},
            "lead_like": lead_like,
            "frag_like": frag_like,
            "qed":       round(qed, 3),
            "sas":       sas,
            "sas_class": sas_class,
        },
        "absorption": {
            "bioavailability_score": round(bio_score, 2),
            "oral_class": ("Good" if bio_score > 0.55 else
                           "Moderate" if bio_score > 0.35 else "Poor"),
            "hia_percent":  round(hia, 1),
            "caco2_logPapp": round(log_papp, 2),
            "caco2_papp":   round(papp, 4),
            "caco2_class":  caco2_class,
            "pampa_logP":   round(log_pampa, 2),
            "pampa_class":  pampa_class,
        },
        "distribution": {
            "logBB":         round(log_bb, 2),
            "bbb_penetrant": bbb_pass,
            "ppb_percent":   round(ppb_est, 1),
            "ppb_class":     ppb_class,
            "vd_lkg":        round(vd_est, 2),
            "vd_class":      vd_class,
        },
        "metabolism": {
            "t_half_hours": t_half,
            "cyp_alerts":   cyp_alerts,
            "cyp_risk":     ("High" if len(cyp_alerts) > 1 else
                             "Medium" if len(cyp_alerts) == 1 else "Low"),
        },
        "excretion": {
            "logS_esol":    round(log_s, 2),
            "solubility_mg_mL": round(sol_mg_ml, 4),
            "solubility_class": sol_class,
        },
        "toxicity": {
            "herg_risk":    herg_risk,
            "pains_alerts": pains_alerts,
            "ames_alert":   ames_alert,
        },
    }


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def _bioavailability_score(mw, logp, hbd, hba, tpsa, rotbonds):
    """Composite bioavailability score 0–1"""
    s = 1.0
    if mw    > 500: s -= 0.20
    elif mw  > 400: s -= 0.10
    if logp  > 5:   s -= 0.20
    elif logp < -2: s -= 0.10
    if hbd   > 5:   s -= 0.20
    elif hbd > 3:   s -= 0.10
    if tpsa  > 140: s -= 0.30
    elif tpsa > 90: s -= 0.10
    if rotbonds > 10: s -= 0.10
    return round(max(0.0, min(1.0, s)), 2)


def _calc_qed(mw, logp, hbd, hba, tpsa, rotbonds, arom_rings):
    """Simplified QED score (Bickerton 2012)"""
    def d(x, a, b, c, dd):
        if x <= a or x >= dd: return 0.05
        if x <= b: return (x - a) / (b - a)
        if x <= c: return 1.0
        return (dd - x) / (dd - c)
    scores  = [d(mw, 0, 100, 400, 600),
               d(logp, -3, 0, 3.5, 7),
               d(hbd, 0, 0, 2, 6),
               d(hba, 0, 0, 6, 14),
               d(tpsa, 0, 0, 90, 180),
               d(rotbonds, 0, 0, 6, 15),
               d(arom_rings, 0, 0, 3, 6)]
    weights = [0.66, 0.56, 0.05, 0.39, 0.17, 0.15, 0.24]
    total_w = sum(weights)
    qed = math.exp(sum(w * math.log(s) for w, s in zip(weights, scores)) / total_w)
    return min(1.0, max(0.0, qed))


def _check_cyp_alerts(mol):
    alerts = []
    patterns = {
        "Naphthoquinone (CYP1A2/3A4 substrate)": "O=C1C=CC(=O)c2ccccc12",
        "Thiocyanate (CYP inhibitor)":            "[SX2][C]#[N]",
        "Thioether (CYP2C9 substrate)":           "[SX2]([#6])[#6]",
        "Quinone Michael acceptor":               "O=C1C=CC(=O)C=1",
    }
    for label, smarts in patterns.items():
        pat = Chem.MolFromSmarts(smarts)
        if pat and mol.HasSubstructMatch(pat):
            alerts.append(label)
    return alerts


def _herg_risk(logp, n_rings, hba):
    risk = (logp > 3) + (n_rings > 2) + (hba > 4)
    return "High" if risk >= 3 else ("Medium" if risk == 2 else "Low")


def _pains_filter(mol):
    """Check for common PAINS alerts"""
    alerts = []
    pains_smarts = {
        "Quinone":            "O=C1C=CC(=O)c2ccccc12",
        "Michael acceptor":   "C=CC(=O)",
        "Electrophilic S":    "[SX2][C]#[N]",
    }
    for name, smarts in pains_smarts.items():
        pat = Chem.MolFromSmarts(smarts)
        if pat and mol.HasSubstructMatch(pat):
            alerts.append(name)
    return alerts


def _ames_alert(mol):
    """Simple AMES mutagenicity structural alert"""
    ames_smarts = [
        "[N+](=O)[O-]",       # Nitro group
        "c1ccc2cc3ccccc3cc2c1",  # Polycyclic aromatic
        "[N]=[N+]=[N-]",       # Azide
        "C(=S)N",              # Thioamide
    ]
    for smarts in ames_smarts:
        pat = Chem.MolFromSmarts(smarts)
        if pat and mol.HasSubstructMatch(pat):
            return True
    return False


# =============================================================================
# FIGURE GENERATION
# =============================================================================

def generate_mol_image(smiles: str, output_path: Path) -> bool:
    """Render 2D structure image"""
    if not RDKIT_OK:
        return False
    try:
        mol = Chem.MolFromSmiles(smiles)
        mol = Draw.PrepareMolForDrawing(mol)
        drawer = rdMolDraw2D.MolDraw2DSVG(400, 300)
        drawer.drawOptions().addStereoAnnotation = True
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        output_path.write_text(drawer.GetDrawingText())
        return True
    except Exception as e:
        print(f"[WARN] Mol image: {e}")
        return False


def fig_radar_chart(admet: dict, output_dir: Path) -> Path:
    """
    Spider/radar chart of the 6 key ADMET dimensions.
    Normalized 0–1 (1 = ideal drug-like value).
    """
    pc  = admet["physicochemical"]
    ab  = admet["absorption"]
    di  = admet["distribution"]

    # Normalize each property to 0-1 (1 = best for oral drug)
    mw_score   = max(0, 1 - max(0, pc["MW"] - 160) / 340)           # best 160-500
    logp_score = max(0, 1 - abs(pc["logP"] - 2) / 3)                # best ≈ 2
    tpsa_score = max(0, 1 - max(0, pc["TPSA"] - 60) / 80)          # best <60
    hbd_score  = max(0, 1 - max(0, pc["HBD"] - 1) / 4)             # best ≤1
    bio_score  = ab["bioavailability_score"]
    caco_score = max(0, 1 - max(0, -ab["caco2_logPapp"] - 4.5) / 2) # best > -4.5

    categories = ['MW\nScore', 'LogP\nScore', 'TPSA\nScore',
                  'HBD\nScore', 'Bioavail.\nScore', 'Caco-2\nScore']
    values = [mw_score, logp_score, tpsa_score, hbd_score, bio_score, caco_score]

    N = len(categories)
    angles = [n / float(N) * 2 * math.pi for n in range(N)]
    angles += angles[:1]
    vals   = values + values[:1]

    fig, ax = plt.subplots(figsize=(6, 6), subplot_kw=dict(polar=True))
    ax.set_facecolor('#f8f9fa')
    fig.patch.set_facecolor('white')

    # Reference circles
    for r in [0.25, 0.5, 0.75, 1.0]:
        circle_pts = np.linspace(0, 2 * math.pi, 100)
        ax.plot(circle_pts, [r] * 100, '--', color='#cccccc', linewidth=0.7, alpha=0.6)

    # Ideal drug zone (>0.6 on all axes)
    ideal = [0.6] * N + [0.6]
    ax.fill(angles, ideal, alpha=0.10, color='#28a745')
    ax.plot(angles, ideal, '--', color='#28a745', linewidth=1, alpha=0.5, label='Drug-like zone')

    # Data
    ax.fill(angles, vals, alpha=0.25, color='#1f6feb')
    ax.plot(angles, vals, 'o-', color='#1f6feb', linewidth=2.5, markersize=7, label=COMPOUND_ABBR)

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories, fontsize=10, fontweight='bold')
    ax.set_ylim(0, 1)
    ax.set_yticks([])

    # Score annotations
    for angle, val, cat in zip(angles[:-1], values, categories):
        ax.annotate(f'{val:.2f}', xy=(angle, val),
                    xytext=(angle, val + 0.10),
                    ha='center', va='center',
                    fontsize=8, fontweight='bold', color='#1f6feb')

    ax.set_title(f'ADMET Radar Chart\n{COMPOUND_NAME}',
                 pad=20, fontsize=12, fontweight='bold')
    ax.legend(loc='upper right', bbox_to_anchor=(1.35, 1.12), fontsize=9)
    plt.tight_layout()

    out = output_dir / "ADMET_radar.png"
    plt.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    return out


def fig_lipinski_dashboard(admet: dict, output_dir: Path) -> Path:
    """
    Visual dashboard of Lipinski Rule of 5 + Veber rules.
    Traffic-light style gauges for each property.
    """
    pc   = admet["physicochemical"]
    lip  = admet["drug_likeness"]["lipinski"]
    veb  = admet["drug_likeness"]["veber"]

    # (label, value, limit, unit, is_below_good)
    props = [
        ("Molecular Weight",  pc["MW"],              500,  "Da",        True),
        ("cLogP",             pc["logP"],             5.0,  "",          True),
        ("H-Bond Donors",     pc["HBD"],              5,    "",          True),
        ("H-Bond Acceptors",  pc["HBA"],              10,   "",          True),
        ("TPSA",              pc["TPSA"],             140,  "Å²",        True),
        ("Rotatable Bonds",   pc["rotatable_bonds"],  10,   "",          True),
    ]

    fig, axes = plt.subplots(2, 3, figsize=(12, 7))
    fig.suptitle(f'Drug-Likeness Dashboard — {COMPOUND_NAME}',
                 fontsize=13, fontweight='bold', y=1.01)
    axes = axes.flatten()

    for ax, (label, val, limit, unit, below_good) in zip(axes, props):
        ok = val <= limit if below_good else val >= limit
        color   = '#28a745' if ok else '#dc3545'
        bg_col  = '#e8f5e9' if ok else '#fdecea'
        bar_col = color

        ax.set_facecolor(bg_col)
        pct = min(1.0, val / limit) if limit > 0 else 0
        ax.barh([0], [pct], color=bar_col, alpha=0.6, height=0.5)
        ax.barh([0], [1.0], color='#e0e0e0', alpha=0.3, height=0.5)
        ax.axvline(x=1.0, color='#888', linewidth=1.5, linestyle='--', alpha=0.8)

        lbl = f'{val:.1f}' if isinstance(val, float) else str(val)
        ax.text(0.5, 0, f'{lbl} {unit}', ha='center', va='center',
                fontsize=16, fontweight='bold', color=color, transform=ax.transAxes)
        ax.text(0.5, 0.75, label, ha='center', va='center',
                fontsize=10, fontweight='bold', transform=ax.transAxes)
        lim_text = f'Limit: ≤{limit} {unit}' if below_good else f'Min: ≥{limit} {unit}'
        ax.text(0.5, 0.15, lim_text, ha='center', va='center',
                fontsize=8, color='#555', transform=ax.transAxes)
        badge = '✓ PASS' if ok else '✗ FAIL'
        badge_col = '#28a745' if ok else '#dc3545'
        ax.text(0.5, -0.25, badge, ha='center', va='center',
                fontsize=12, fontweight='bold', color=badge_col, transform=ax.transAxes)
        ax.set_xlim(0, 1.2)
        ax.set_ylim(-0.4, 0.6)
        ax.axis('off')

    # Summary text
    violations = admet["drug_likeness"]["lipinski"]["violations"]
    lip_ok  = admet["drug_likeness"]["lipinski"]["pass"]
    veber_ok = admet["drug_likeness"]["veber"]["pass"]
    qed     = admet["drug_likeness"]["qed"]

    summary = (f"Lipinski: {'PASS ✓' if lip_ok else 'FAIL ✗'} ({violations} violations)  |  "
               f"Veber: {'PASS ✓' if veber_ok else 'FAIL ✗'}  |  QED = {qed:.3f}")
    fig.text(0.5, -0.01, summary, ha='center', fontsize=11, fontweight='bold',
             color='#28a745' if (lip_ok and veber_ok) else '#dc3545')

    plt.tight_layout()
    out = output_dir / "ADMET_lipinski.png"
    plt.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    return out


def fig_admet_summary(admet: dict, output_dir: Path) -> Path:
    """
    Comprehensive ADMET summary figure (4-panel):
    A) Bioavailability wheel  B) Absorption radar  C) Key properties table  D) SwissADME-style
    """
    pc = admet["physicochemical"]
    ab = admet["absorption"]
    di = admet["distribution"]
    me = admet["metabolism"]
    ex = admet["excretion"]
    to = admet["toxicity"]
    dl = admet["drug_likeness"]

    fig = plt.figure(figsize=(16, 10))
    gs  = GridSpec(2, 3, figure=fig, hspace=0.4, wspace=0.35)

    # ── Panel A: Bioavailability wheel ────────────────────────────────────────
    ax_a = fig.add_subplot(gs[0, 0])
    rules = {
        'Lipinski': admet["drug_likeness"]["lipinski"]["pass"],
        'Veber':    admet["drug_likeness"]["veber"]["pass"],
        'Ghose':    admet["drug_likeness"]["ghose"]["pass"],
        'Egan':     admet["drug_likeness"]["egan"]["pass"],
        'Muegge':   admet["drug_likeness"]["muegge"]["pass"],
    }
    n_pass = sum(rules.values())
    colors = ['#28a745' if v else '#dc3545' for v in rules.values()]
    wedge_colors = ['#28a745' if v else '#dc3545' for v in rules.values()]
    wedges, texts = ax_a.pie([1] * 5, labels=list(rules.keys()),
                              colors=wedge_colors, startangle=90,
                              wedgeprops=dict(width=0.6, edgecolor='white', linewidth=2))
    for t in texts:
        t.set_fontsize(9)
        t.set_fontweight('bold')
    ax_a.text(0, 0, f'{n_pass}/5\nPassed', ha='center', va='center',
              fontsize=14, fontweight='bold', color='#222')
    ax_a.set_title('(A) Drug-Likeness Filters', fontweight='bold', fontsize=11, pad=12)

    # ── Panel B: ADME properties bar chart ────────────────────────────────────
    ax_b = fig.add_subplot(gs[0, 1])
    props_b = {
        'Oral Bioavail.': ab["bioavailability_score"],
        'Caco-2 (norm)':  max(0, min(1, (ab["caco2_logPapp"] + 7) / 3)),
        'PAMPA (norm)':   max(0, min(1, (ab["pampa_logP"]    + 8) / 3)),
        'HIA/100':        ab["hia_percent"] / 100,
        'BBB (norm)':     max(0, min(1, (di["logBB"] + 2) / 3)),
    }
    bar_names = list(props_b.keys())
    bar_vals  = list(props_b.values())
    bar_cols  = ['#28a745' if v > 0.55 else ('#ffc107' if v > 0.35 else '#dc3545')
                 for v in bar_vals]
    bars = ax_b.barh(bar_names, bar_vals, color=bar_cols, edgecolor='white', linewidth=0.5)
    ax_b.axvline(x=0.55, color='#28a745', linestyle='--', linewidth=1.5, alpha=0.7)
    ax_b.set_xlim(0, 1.1)
    ax_b.set_xlabel('Normalized Score (0–1)', fontsize=9)
    ax_b.set_title('(B) Absorption & Distribution', fontweight='bold', fontsize=11)
    for bar, val in zip(bars, bar_vals):
        ax_b.text(val + 0.02, bar.get_y() + bar.get_height() / 2,
                  f'{val:.2f}', va='center', fontsize=9, fontweight='bold')
    ax_b.spines['top'].set_visible(False)
    ax_b.spines['right'].set_visible(False)

    # ── Panel C: Metabolism / Toxicity flags ─────────────────────────────────
    ax_c = fig.add_subplot(gs[0, 2])
    ax_c.axis('off')
    rows = [
        ['Property',          'Value',                        'Status'],
        ['CYP Risk',          me["cyp_risk"],                 '⚠' if me["cyp_risk"] != 'Low' else '✓'],
        ['hERG Risk',         to["herg_risk"],                '⚠' if to["herg_risk"] != 'Low' else '✓'],
        ['AMES Alert',        'Yes' if to["ames_alert"] else 'No', '⚠' if to["ames_alert"] else '✓'],
        ['PAINS Alerts',      str(len(to["pains_alerts"])),   '⚠' if to["pains_alerts"] else '✓'],
        ['T½ (est.)',         f'{me["t_half_hours"]:.1f} h',  '✓'],
        ['PPB',               di["ppb_class"][:10],           '—'],
        ['Vd (est.)',         f'{di["vd_lkg"]:.1f} L/kg',    '—'],
    ]
    colors_tbl = [['#e3f2fd', '#e3f2fd', '#e3f2fd']]
    for r in rows[1:]:
        status_col = ('#e8f5e9' if r[2] == '✓' else
                      '#fff8e1' if r[2] == '—' else '#fdecea')
        colors_tbl.append(['#fafafa', '#fafafa', status_col])

    tbl = ax_c.table(cellText=rows, cellLoc='center', loc='center',
                     cellColours=colors_tbl)
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1, 1.5)
    for (r, c), cell in tbl.get_celld().items():
        if r == 0:
            cell.set_text_props(fontweight='bold', color='white')
            cell.set_facecolor('#1565C0')
        cell.set_edgecolor('#dddddd')
    ax_c.set_title('(C) Metabolism & Safety', fontweight='bold', fontsize=11, pad=12)

    # ── Panel D: Physicochemical property table ────────────────────────────────
    ax_d = fig.add_subplot(gs[1, :])
    ax_d.axis('off')

    params = [
        ['Parameter', 'Value', 'Optimal Range', 'Assessment'],
        ['Molecular Weight (Da)', f'{pc["MW"]:.2f}', '160–500', '✓ PASS' if pc["MW"] <= 500 else '✗ FAIL'],
        ['cLogP', f'{pc["logP"]:.2f}', '0–5', '✓ PASS' if pc["logP"] <= 5 else '✗ FAIL'],
        ['TPSA (Å²)', f'{pc["TPSA"]:.1f}', '< 140', '✓ PASS' if pc["TPSA"] < 140 else '✗ FAIL'],
        ['H-Bond Donors', str(pc["HBD"]), '≤ 5', '✓ PASS' if pc["HBD"] <= 5 else '✗ FAIL'],
        ['H-Bond Acceptors', str(pc["HBA"]), '≤ 10', '✓ PASS' if pc["HBA"] <= 10 else '✗ FAIL'],
        ['Rotatable Bonds', str(pc["rotatable_bonds"]), '≤ 10', '✓ PASS' if pc["rotatable_bonds"] <= 10 else '✗ FAIL'],
        ['QED Drug-Likeness', f'{dl["qed"]:.3f}', '> 0.5 (good)', '✓ Good' if dl["qed"] > 0.5 else '▲ Moderate'],
        ['Water Solubility', ex["solubility_class"], 'Soluble or better', '—'],
        ['Oral Bioavailability', f'{ab["bioavailability_score"]:.2f}', '> 0.55', '✓ PASS' if ab["bioavailability_score"] > 0.55 else '▲ Moderate'],
        ['%HIA (estimated)', f'{ab["hia_percent"]:.1f}%', '> 70%', '✓ Good' if ab["hia_percent"] > 70 else '▲ Moderate'],
        ['Caco-2 Class', ab["caco2_class"], 'High or Medium', '✓ OK' if ab["caco2_class"] != 'Low' else '⚠ Low'],
    ]

    row_colors = [['#1565C0'] * 4]
    for p in params[1:]:
        if '✓' in p[3]:
            row_colors.append(['#fafafa', '#f0f7ff', '#f0f7ff', '#e8f5e9'])
        elif '✗' in p[3]:
            row_colors.append(['#fafafa', '#f0f7ff', '#f0f7ff', '#fdecea'])
        elif '⚠' in p[3]:
            row_colors.append(['#fafafa', '#f0f7ff', '#f0f7ff', '#fff8e1'])
        else:
            row_colors.append(['#fafafa', '#f0f7ff', '#f0f7ff', '#fafafa'])

    tbl_d = ax_d.table(cellText=params, cellLoc='center', loc='center',
                        cellColours=row_colors)
    tbl_d.auto_set_font_size(False)
    tbl_d.set_fontsize(10)
    tbl_d.scale(1, 1.6)
    for (r, c), cell in tbl_d.get_celld().items():
        if r == 0:
            cell.set_text_props(fontweight='bold', color='white', fontsize=10)
        cell.set_edgecolor('#dddddd')
    ax_d.set_title('(D) Full ADMET Parameter Table', fontweight='bold',
                   fontsize=11, pad=5)

    fig.suptitle(f'In Silico ADMET Profile — {COMPOUND_NAME} ({COMPOUND_ABBR})',
                 fontsize=14, fontweight='bold', y=1.01)

    out = output_dir / "ADMET_full_summary.png"
    plt.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    return out


# =============================================================================
# PAPER TEXT GENERATION
# =============================================================================

def generate_admet_text(admet: dict, lang: str = "en") -> str:
    """Generate manuscript-ready ADMET text (Results section)"""

    pc = admet["physicochemical"]
    ab = admet["absorption"]
    di = admet["distribution"]
    me = admet["metabolism"]
    ex = admet["excretion"]
    to = admet["toxicity"]
    dl = admet["drug_likeness"]

    if lang == "pt":
        return _admet_text_pt(pc, ab, di, me, ex, to, dl)
    return _admet_text_en(pc, ab, di, me, ex, to, dl)


def _admet_text_en(pc, ab, di, me, ex, to, dl):
    lip   = dl["lipinski"]
    veber = dl["veber"]
    n_rules_pass = sum([lip["pass"], veber["pass"], dl["ghose"]["pass"],
                        dl["egan"]["pass"], dl["muegge"]["pass"]])

    text = f"""In Silico ADMET Profiling

The in silico pharmacokinetic (ADMET) profile of {COMPOUND_NAME} ({COMPOUND_ABBR}) \
was evaluated using RDKit-based classical QSAR models. The compound has a molecular \
weight of {pc['MW']:.2f} Da, a calculated logP of {pc['logP']:.2f}, {pc['HBD']} \
hydrogen-bond donors, {pc['HBA']} hydrogen-bond acceptors, and a topological polar \
surface area (TPSA) of {pc['TPSA']:.1f} Å².

Drug-Likeness: {COMPOUND_ABBR} satisfied {n_rules_pass} of 5 drug-likeness filter \
sets assessed (Lipinski, Veber, Ghose, Egan, and Muegge). Lipinski's Rule of Five \
was {'met' if lip['pass'] else 'violated'} ({lip['violations']} violation{'s' if lip['violations'] != 1 else ''}). \
Veber's oral bioavailability criteria (TPSA ≤ 140 Å² and ≤ 10 rotatable bonds) were \
{'satisfied' if veber['pass'] else 'not satisfied'}. The QED drug-likeness score was \
{dl['qed']:.3f}, {'indicating good drug-like character' if dl['qed'] > 0.5 else 'indicating moderate drug-like character'}.

Absorption: The predicted oral bioavailability score was {ab['bioavailability_score']:.2f}, \
corresponding to {'good' if ab['oral_class'] == 'Good' else ab['oral_class'].lower()} \
oral absorption. The estimated Caco-2 apparent permeability coefficient \
(logPapp = {ab['caco2_logPapp']:.2f} × 10⁻⁶ cm/s) indicates {ab['caco2_class'].lower()} \
permeability across the intestinal epithelium. The predicted human intestinal absorption \
(%HIA) was {ab['hia_percent']:.1f}%.

Distribution: The predicted blood–brain barrier (BBB) penetration index (logBB = \
{di['logBB']:.2f}) suggests the compound is {'likely' if di['bbb_penetrant'] else 'unlikely'} \
to cross the BBB, which is {'advantageous for avoiding CNS side effects' if not di['bbb_penetrant'] else 'a consideration for CNS-related toxicity'}. \
Estimated plasma protein binding was {di['ppb_percent']:.1f}% ({di['ppb_class']}).

Metabolism and Excretion: Structural analysis identified the following potential CYP450 \
interaction alerts: {', '.join(me['cyp_alerts']) if me['cyp_alerts'] else 'none'}. \
The estimated metabolic half-life is {me['t_half_hours']:.1f} hours. Water solubility \
was estimated as {ex['solubility_class'].lower()} \
(logS = {ex['logS_esol']:.2f}, ≈ {ex['solubility_mg_mL']:.4f} mg/mL).

Safety: The hERG cardiac toxicity risk was predicted as {to['herg_risk'].lower()}. \
{'No AMES mutagenicity structural alerts were detected.' if not to['ames_alert'] else 'One or more AMES mutagenicity structural alerts were detected; experimental validation is recommended.'} \
Collectively, the ADMET profile of {COMPOUND_ABBR} suggests \
{'a favorable pharmacokinetic profile for oral administration' if ab['bioavailability_score'] > 0.55 else 'moderate pharmacokinetic properties that may require formulation optimization for oral delivery'}, \
consistent with its evaluation as an oral antiparasitic candidate."""

    return text


def _admet_text_pt(pc, ab, di, me, ex, to, dl):
    lip   = dl["lipinski"]
    veber = dl["veber"]
    n_rules_pass = sum([lip["pass"], veber["pass"], dl["ghose"]["pass"],
                        dl["egan"]["pass"], dl["muegge"]["pass"]])

    text = f"""Perfil ADMET In Silico

O perfil farmacocinético in silico (ADMET) de {COMPOUND_NAME} ({COMPOUND_ABBR}) \
foi avaliado utilizando modelos QSAR clássicos baseados em RDKit. O composto possui \
peso molecular de {pc['MW']:.2f} Da, logP calculado de {pc['logP']:.2f}, {pc['HBD']} \
doadores de ligação de hidrogênio, {pc['HBA']} aceptores de ligação de hidrogênio \
e área de superfície polar topológica (TPSA) de {pc['TPSA']:.1f} Å².

Drug-likeness: A {COMPOUND_ABBR} satisfez {n_rules_pass} de 5 filtros de \
drug-likeness avaliados (Lipinski, Veber, Ghose, Egan e Muegge). A Regra dos \
Cinco de Lipinski foi {'cumprida' if lip['pass'] else 'violada'} ({lip['violations']} \
violação{'ões' if lip['violations'] != 1 else ''}). Os critérios de biodisponibilidade \
oral de Veber (TPSA ≤ 140 Å² e ≤ 10 ligações rotáveis) foram \
{'satisfeitos' if veber['pass'] else 'não satisfeitos'}. A pontuação QED de \
drug-likeness foi {dl['qed']:.3f}.

Absorção: O escore de biodisponibilidade oral previsto foi de {ab['bioavailability_score']:.2f}, \
correspondendo a absorção oral {'boa' if ab['oral_class'] == 'Good' else ab['oral_class'].lower()}. \
O coeficiente de permeabilidade aparente Caco-2 estimado (logPapp = {ab['caco2_logPapp']:.2f}) \
indica permeabilidade {ab['caco2_class'].lower()} através do epitélio intestinal. \
A absorção intestinal humana prevista (%HIA) foi de {ab['hia_percent']:.1f}%.

Distribuição: O índice de penetração na barreira hematoencefálica (BHE) previsto \
(logBB = {di['logBB']:.2f}) sugere que o composto é \
{'provável' if di['bbb_penetrant'] else 'improvável'} de atravessar a BHE. \
A ligação proteica plasmática estimada foi de {di['ppb_percent']:.1f}% ({di['ppb_class']}).

Metabolismo e Excreção: A análise estrutural identificou os seguintes alertas \
potenciais de interação com CYP450: {', '.join(me['cyp_alerts']) if me['cyp_alerts'] else 'nenhum'}. \
A meia-vida metabólica estimada é de {me['t_half_hours']:.1f} horas. A solubilidade \
em água foi estimada como {ex['solubility_class'].lower()} (logS = {ex['logS_esol']:.2f}).

Segurança: O risco de toxicidade cardíaca (hERG) foi previsto como \
{to['herg_risk'].lower()}. {'Nenhum alerta de mutagenicidade AMES foi detectado.' if not to['ames_alert'] else 'Um ou mais alertas estruturais de mutagenicidade AMES foram detectados; validação experimental é recomendada.'} \
O perfil ADMET geral da {COMPOUND_ABBR} sugere \
{'um perfil farmacocinético favorável para administração oral' if ab['bioavailability_score'] > 0.55 else 'propriedades farmacocinéticas moderadas que podem requerer otimização de formulação para via oral'}."""

    return text


# =============================================================================
# RUN MODULE STANDALONE
# =============================================================================

def run_admet(smiles: str = SMILES_DEFAULT, output_dir: str = "docking_project/figures",
              lang: str = "en",
              compound_name: str = None,
              compound_abbr: str = None) -> dict:
    """Main entry point: calculate ADMET and generate all figures + text"""

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Override module-level name with user-provided compound name
    global COMPOUND_NAME, COMPOUND_ABBR
    # Serializa a regiao de mutacao/uso/restauracao dos globals para evitar
    # entrelacamento de rotulos sob Flask multithread (race condition).
    with _LOCK:
        _orig_name = COMPOUND_NAME
        _orig_abbr = COMPOUND_ABBR
        if compound_name:
            COMPOUND_NAME = compound_name
            COMPOUND_ABBR = compound_abbr or compound_name
        try:
            print(f"\n{'='*60}")
            print(f"  ADMET Analysis: {COMPOUND_NAME}")
            print(f"  SMILES: {smiles}")
            print(f"{'='*60}\n")

            admet = calculate_admet(smiles)
            if "error" in admet:
                print(f"[ERRO] {admet['error']}")
                return admet

            # Print summary
            pc = admet["physicochemical"]
            dl = admet["drug_likeness"]
            ab = admet["absorption"]

            print(f"  MW    = {pc['MW']:.2f} Da       | logP = {pc['logP']:.2f}")
            print(f"  HBD   = {pc['HBD']}              | HBA  = {pc['HBA']}")
            print(f"  TPSA  = {pc['TPSA']:.1f} Å²     | RB   = {pc['rotatable_bonds']}")
            print(f"  QED   = {dl['qed']:.3f}          | SAS  = {dl['sas']:.2f}")
            lip_ok = dl["lipinski"]["pass"]
            print(f"  Lipinski: {'PASS ✓' if lip_ok else 'FAIL ✗'}  |  Veber: {'PASS ✓' if dl['veber']['pass'] else 'FAIL ✗'}")
            print(f"  Bioavailability: {ab['bioavailability_score']:.2f}  |  Caco-2: {ab['caco2_class']}")

            # Generate figures
            print("\n  Generating ADMET figures...")
            try:
                f1 = fig_radar_chart(admet, out_dir)
                print(f"  [OK] {f1.name}")
                f2 = fig_lipinski_dashboard(admet, out_dir)
                print(f"  [OK] {f2.name}")
                f3 = fig_admet_summary(admet, out_dir)
                print(f"  [OK] {f3.name}")
            except Exception as e:
                print(f"  [WARN] Figure generation error: {e}")

            # Generate SVG molecular structure
            safe_abbr = COMPOUND_ABBR.replace(" ", "_").replace("-", "").replace(",", "")
            mol_svg = out_dir / f"{safe_abbr}_structure.svg"
            if generate_mol_image(smiles, mol_svg):
                print(f"  [OK] {mol_svg.name}")

            # Generate text
            text = generate_admet_text(admet, lang=lang)
            text_file = out_dir / f"admet_text_{lang}.txt"
            text_file.write_text(text, encoding='utf-8')
            print(f"  [OK] {text_file.name}")

            # Save JSON
            json_file = out_dir / "admet_results.json"
            json_file.write_text(json.dumps(admet, indent=2, ensure_ascii=False), encoding='utf-8')
            print(f"  [OK] {json_file.name}")

            print(f"\n  All ADMET outputs saved to: {out_dir}")
            return admet
        finally:
            # Restore original global names to avoid side effects across calls
            COMPOUND_NAME = _orig_name
            COMPOUND_ABBR = _orig_abbr


# =============================================================================
# CLI
# =============================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ADMET In Silico Analysis")
    parser.add_argument("--smiles",  default=SMILES_DEFAULT)
    parser.add_argument("--output",  default="docking_project/figures")
    parser.add_argument("--lang",    default="en", choices=["en", "pt"])
    args = parser.parse_args()
    run_admet(args.smiles, args.output, args.lang)
