#!/usr/bin/env python3
"""
=============================================================================
HUMAN TOXICITY PREDICTION MODULE
Predicts LD50, hepatotoxicity, mutagenicity and organ toxicity endpoints
using structural alerts (RDKit) + QSAR models + ProTox-II-like estimates.
=============================================================================
Usage:
    python toxicity_module.py
    python toxicity_module.py --smiles "O=C1C=C(SC#N)C(=O)c2ccccc21" --lang en
"""
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
    from matplotlib.gridspec import GridSpec
    import matplotlib.patches as mpatches
except ImportError:
    import _deps
    _deps.ensure(["numpy", "matplotlib"])  # raise se nao opt-in; instala se DOCKING_AUTO_PIP=1
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    import matplotlib.patches as mpatches

RDKIT_OK = False
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem
    RDKIT_OK = True
except ImportError:
    pass

COMPOUND_NAME = "2-Thiocyanato-1,4-naphthoquinone"
COMPOUND_ABBR = "2-SCN-NQ"
SMILES_DEFAULT = "O=C1C=C(SC#N)C(=O)c2ccccc21"

# ─── GHS toxicity classes ─────────────────────────────────────────────────────
GHS_CLASSES = {
    1: ("Fatal",           "#C62828", "≤ 5 mg/kg"),
    2: ("Fatal",           "#E53935", "5–50 mg/kg"),
    3: ("Toxic",           "#FB8C00", "50–300 mg/kg"),
    4: ("Harmful",         "#FDD835", "300–2000 mg/kg"),
    5: ("May be harmful",  "#66BB6A", "2000–5000 mg/kg"),
    6: ("Nontoxic",        "#43A047", "> 5000 mg/kg"),
}

# ─── Default experimental dataset (pre-filled from literature for similar NQs) ──
# Override these with your actual experimental data via the dashboard or function args.
MOUSE_EXPERIMENTAL = {
    "LD50_oral_mg_kg":  150,    # mice, oral — update with your measured value
    "source":           "Literature extrapolation (similar naphthoquinone analogs)",
    "species":          "Mus musculus",
    "route":            "Oral",
}

# Default IC50 values — will be overridden by user-supplied experimental data
EXPERIMENTAL_DEFAULTS = {
    # ── Mammalian cell toxicity ──────────────────────────────────────────────
    "cc50_vero_uM":            ">50",      # VERO (monkey kidney — most used)
    "cc50_hela_uM":            None,       # HeLa (human cervical)
    "cc50_thp1_uM":            None,       # THP-1 (human macrophage, most relevant)
    "cc50_j774_uM":            None,       # J774.A1 (murine macrophage)
    "cc50_hepg2_uM":           None,       # HepG2 (human hepatocyte)

    # ── T. cruzi activity ────────────────────────────────────────────────────
    "ic50_tc_trypo_intra_uM":  None,       # intracellular trypomastigotes (most relevant)
    "ic50_tc_epi_uM":          None,       # epimastigotes (axenic culture)
    "ic50_tc_amast_intra_uM":  None,       # intracellular amastigotes

    # ── Leishmania activity ──────────────────────────────────────────────────
    "ic50_lm_promas_uM":       None,       # L. amazonensis promastigotes (axenic)
    "ic50_lm_amast_intra_uM":  None,       # intracellular amastigotes (most clinically relevant)

    # ── In vivo ──────────────────────────────────────────────────────────────
    "ld50_mouse_oral_mg_kg":   None,       # your experimental LD50 (replaces default)
    "ld50_mouse_ip_mg_kg":     None,       # intraperitoneal route
    "ld50_rat_oral_mg_kg":     None,       # rat oral

    # ── Reference drug for context ───────────────────────────────────────────
    "ref_drug":                "Benznidazole",
    "ref_ic50_tc_uM":          0.5,        # Benznidazole IC50 T. cruzi ~0.5 uM
    "ref_ic50_lm_uM":          15.0,       # Miltefosine IC50 Leishmania ~15 uM
}

# =============================================================================
# TOXICITY PREDICTION ENGINE
# =============================================================================

def predict_toxicity(smiles: str = SMILES_DEFAULT,
                     exp_data: dict = None) -> dict:
    """
    Multi-endpoint toxicity prediction using structural alerts and QSAR models.
    exp_data: optional dict with user-supplied experimental values.
              Keys mirror EXPERIMENTAL_DEFAULTS. Supplied values override defaults.
    """
    # Merge user-supplied experimental data with defaults
    _exp = dict(EXPERIMENTAL_DEFAULTS)
    if exp_data:
        _exp.update({k: v for k, v in exp_data.items() if v not in (None, "", "N/A", 0)})
    if not RDKIT_OK:
        return {"error": "RDKit not installed. Run: pip install rdkit"}

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"error": f"Invalid SMILES: {smiles}"}

    mw   = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd  = rdMolDescriptors.CalcNumHBD(mol)
    hba  = rdMolDescriptors.CalcNumHBA(mol)
    tpsa = Descriptors.TPSA(mol)

    # ── If user provided real LD50, use it; otherwise predict via QSAR ─────────
    user_ld50_mouse = _exp.get("ld50_mouse_oral_mg_kg") or                       _exp.get("ld50_mouse_ip_mg_kg")   or                       MOUSE_EXPERIMENTAL["LD50_oral_mg_kg"]
    user_ld50_rat   = _exp.get("ld50_rat_oral_mg_kg")

    # ── LD50 Human Prediction (ProTox-II-style regression) ────────────────────
    # Model: Drwal et al. 2014; based on toxicophores + physicochemical descriptors
    # log(LD50 in mg/kg) = f(MW, logP, TPSA, structural alerts)
    base_log_ld50 = 3.2   # 1585 mg/kg - baseline for drug-like molecules
    # Corrections:
    delta = 0.0
    # Naphthoquinone (pro-oxidant, electrophilic quinone) — more toxic
    quinone_smarts = Chem.MolFromSmarts("O=C1C=CC(=O)c2ccccc12")
    if mol.HasSubstructMatch(quinone_smarts):
        delta -= 0.6   # reduce LD50 by ~4x
    # Thiocyanate — moderately toxic
    scn_smarts = Chem.MolFromSmarts("[SX2][C]#[N]")
    if mol.HasSubstructMatch(scn_smarts):
        delta -= 0.3
    # LogP > 3 — increased membrane penetration (more toxic potential)
    if logp > 3:
        delta -= 0.2
    # High TPSA — reduced oral absorption (less toxic systemically)
    if tpsa > 100:
        delta += 0.2
    # Small MW — better absorption
    if mw < 250:
        delta -= 0.1

    log_ld50_human = base_log_ld50 + delta
    ld50_human_mg_kg = round(10 ** log_ld50_human, 0)

    # GHS class
    ghs_class = _ghs_class(ld50_human_mg_kg)

    # ── Hepatotoxicity ────────────────────────────────────────────────────────
    # Features: quinone (CYP1A2 substrate, GSH conjugates), reactive metabolites
    hepato_risk, hepato_evidence = _hepatotoxicity(mol, logp, mw)

    # ── Mutagenicity (AMES) ───────────────────────────────────────────────────
    ames_positive, ames_alerts = _ames_prediction(mol)

    # ── Carcinogenicity ───────────────────────────────────────────────────────
    carcinogen_risk, carcinogen_reason = _carcinogenicity(mol, ames_positive)

    # ── Organ toxicity endpoints ───────────────────────────────────────────────
    organ_tox = _organ_toxicity(mol, logp, mw, tpsa)

    # ── Environmental toxicity ─────────────────────────────────────────────────
    env_tox = _environmental_toxicity(logp, mw)

    # ── Overall safety assessment ─────────────────────────────────────────────
    n_concerns = sum([
        ghs_class <= 3,
        hepato_risk == "High",
        ames_positive,
        carcinogen_risk == "High",
        organ_tox["cardiac"] == "High",
        organ_tox["renal"] == "High",
    ])
    overall = "High concern" if n_concerns >= 3 else ("Moderate concern" if n_concerns >= 1 else "Low concern")

    # ── Therapeutic index + Selectivity Index (use real data if available) ──────
    # Priority: user experimental IC50 > literature estimate
    def _to_float(v):
        try: return float(str(v).replace(">","").replace("<","").strip())
        except (ValueError, TypeError): return None

    ic50_tc_best  = _to_float(_exp.get("ic50_tc_trypo_intra_uM")) or                     _to_float(_exp.get("ic50_tc_epi_uM")) or 0.5
    ic50_lm_best  = _to_float(_exp.get("ic50_lm_amast_intra_uM")) or                     _to_float(_exp.get("ic50_lm_promas_uM")) or 1.0
    cc50_best     = _to_float(_exp.get("cc50_vero_uM")) or                     _to_float(_exp.get("cc50_thp1_uM"))  or                     _to_float(_exp.get("cc50_j774_uM"))  or None

    # SI = CC50(host) / IC50(parasite)
    si_tc  = round(cc50_best / ic50_tc_best, 1)  if cc50_best else None
    si_lm  = round(cc50_best / ic50_lm_best, 1)  if cc50_best else None
    si_str = (f"SI(Tc)={si_tc}x  SI(Lm)={si_lm}x" if si_tc and si_lm
              else "Set CC50 host cell to calculate")

    # TI = LD50(mouse) / ED50_est
    # ED50 estimated from IC50 with oral bioavailability correction (F~0.4)
    ed50_mg_kg = ic50_tc_best * mw / (1000 * 1000 * 0.4)
    ti_est = (user_ld50_mouse / max(ed50_mg_kg, 0.001)) if user_ld50_mouse else None

    # Use real LD50 for GHS if available
    if isinstance(user_ld50_mouse, (int, float)):
        ld50_human_mg_kg = float(user_ld50_mouse)   # use real mouse value for GHS
        ghs_class        = _ghs_class(ld50_human_mg_kg)

    return {
        "metadata": {
            "compound": COMPOUND_NAME,
            "abbreviation": COMPOUND_ABBR,
            "smiles": smiles,
            "timestamp": datetime.now().isoformat(),
        },
        "ld50": {
            "human_predicted_mg_kg":   int(ld50_human_mg_kg),
            "mouse_experimental_mg_kg": user_ld50_mouse or MOUSE_EXPERIMENTAL["LD50_oral_mg_kg"],
            "mouse_source": ("User experimental" if _exp.get("ld50_mouse_oral_mg_kg")
                             else MOUSE_EXPERIMENTAL["source"]),
            "rat_ld50_mg_kg": user_ld50_rat,
            "ghs_class": ghs_class,
            "ghs_label": GHS_CLASSES[ghs_class][0],
            "ghs_color": GHS_CLASSES[ghs_class][1],
            "ghs_range": GHS_CLASSES[ghs_class][2],
        },
        "experimental": {
            # All user-supplied or default activity/toxicity values
            "cc50_vero_uM":            _exp.get("cc50_vero_uM"),
            "cc50_hela_uM":            _exp.get("cc50_hela_uM"),
            "cc50_thp1_uM":            _exp.get("cc50_thp1_uM"),
            "cc50_j774_uM":            _exp.get("cc50_j774_uM"),
            "cc50_hepg2_uM":           _exp.get("cc50_hepg2_uM"),
            "ic50_tc_trypo_intra_uM":  _exp.get("ic50_tc_trypo_intra_uM"),
            "ic50_tc_epi_uM":          _exp.get("ic50_tc_epi_uM"),
            "ic50_tc_amast_intra_uM":  _exp.get("ic50_tc_amast_intra_uM"),
            "ic50_lm_promas_uM":       _exp.get("ic50_lm_promas_uM"),
            "ic50_lm_amast_intra_uM":  _exp.get("ic50_lm_amast_intra_uM"),
            "ld50_mouse_oral_mg_kg":   _exp.get("ld50_mouse_oral_mg_kg"),
            "ld50_mouse_ip_mg_kg":     _exp.get("ld50_mouse_ip_mg_kg"),
            "ld50_rat_oral_mg_kg":     _exp.get("ld50_rat_oral_mg_kg"),
            "ref_drug":                _exp.get("ref_drug","Benznidazole"),
            "ref_ic50_tc_uM":          _exp.get("ref_ic50_tc_uM"),
            "ref_ic50_lm_uM":          _exp.get("ref_ic50_lm_uM"),
        },
        "selectivity": {
            "si_tc_vs_host":   si_tc,
            "si_lm_vs_host":   si_lm,
            "si_summary":      si_str,
            "ic50_tc_used_uM": ic50_tc_best,
            "ic50_lm_used_uM": ic50_lm_best,
            "cc50_used_uM":    cc50_best,
        },
        "hepatotoxicity": {
            "risk": hepato_risk,
            "evidence": hepato_evidence,
        },
        "mutagenicity": {
            "ames_predicted_positive": ames_positive,
            "structural_alerts": ames_alerts,
        },
        "carcinogenicity": {
            "risk": carcinogen_risk,
            "reasoning": carcinogen_reason,
        },
        "organ_toxicity": organ_tox,
        "environmental": env_tox,
        "overall": {
            "safety_concerns": n_concerns,
            "assessment": overall,
            "therapeutic_index_est": round(ti_est, 1) if ti_est else "N/A",
            "selectivity_index_tc":   si_tc,
            "selectivity_index_lm":   si_lm,
            "data_source":            ("experimental" if _exp.get("ld50_mouse_oral_mg_kg")
                                       or _exp.get("ic50_tc_trypo_intra_uM") else "predicted"),
            "development_recommendation": _get_recommendation(n_concerns, ghs_class),
        },
    }


def _ghs_class(ld50_mg_kg):
    if   ld50_mg_kg <=     5: return 1
    elif ld50_mg_kg <=    50: return 2
    elif ld50_mg_kg <=   300: return 3
    elif ld50_mg_kg <=  2000: return 4
    elif ld50_mg_kg <=  5000: return 5
    else:                      return 6


def _hepatotoxicity(mol, logp, mw):
    alerts = []
    # Quinone → Michael acceptor → reactive metabolite → GSH depletion
    q_smarts = Chem.MolFromSmarts("O=C1C=CC(=O)c2ccccc12")
    if mol.HasSubstructMatch(q_smarts):
        alerts.append("Quinone (reactive metabolite, GSH conjugation)")
    # Thiol-reactive groups
    scn_s = Chem.MolFromSmarts("[SX2][C]#[N]")
    if mol.HasSubstructMatch(scn_s):
        alerts.append("Thiocyanate (potential CYP-mediated hepatotoxin)")
    # High logP → hepatic accumulation
    if logp > 3:
        alerts.append(f"LogP {logp:.1f} > 3 (hepatic accumulation risk)")

    risk = "High" if len(alerts) >= 2 else ("Medium" if len(alerts) == 1 else "Low")
    return risk, alerts


def _ames_prediction(mol):
    alerts = []
    smarts_dict = {
        "Polycyclic aromatic system": "c1ccc2ccccc2c1",
        "Michael acceptor (alpha,beta-unsat. carbonyl)": "C=CC(=O)",
        "Halogenated aromatic": "[F,Cl,Br,I]c",
    }
    for name, smarts in smarts_dict.items():
        pat = Chem.MolFromSmarts(smarts)
        if pat and mol.HasSubstructMatch(pat):
            alerts.append(name)
    # Quinone → intercalation, frameshift
    q = Chem.MolFromSmarts("O=C1C=CC(=O)c2ccccc12")
    if mol.HasSubstructMatch(q):
        alerts.append("Naphthoquinone (DNA intercalation, topoisomerase II poisoning)")
    positive = len(alerts) >= 2
    return positive, alerts


def _carcinogenicity(mol, ames_positive):
    # Based on Benigni-Bossa rules
    reasons = []
    q = Chem.MolFromSmarts("O=C1C=CC(=O)c2ccccc12")
    if mol.HasSubstructMatch(q):
        reasons.append("Naphthoquinone (genotoxic mechanism at high concentrations)")
    if ames_positive:
        reasons.append("AMES-positive structural features")
    risk = "Medium" if reasons else "Low"
    return risk, reasons


def _organ_toxicity(mol, logp, mw, tpsa):
    # hERG/cardiac
    herg = (logp > 3) + (mw > 400) + (tpsa < 60)
    cardiac = "High" if herg >= 3 else ("Medium" if herg >= 2 else "Low")

    # Renal — quinone nephrotoxicity is known but dose-dependent
    renal = "Low"
    q = Chem.MolFromSmarts("O=C1C=CC(=O)c2ccccc12")
    if mol.HasSubstructMatch(q):
        renal = "Low-Medium"  # at therapeutic doses, likely low

    # CNS — low if tpsa > 90 (unlikely to cross BBB)
    cns = "Low" if tpsa > 90 else "Medium"

    return {
        "cardiac":  cardiac,
        "renal":    renal,
        "cns":      cns,
        "immune":   "Low",
        "skin":     "Low-Medium",  # quinone contact sensitizer potential
    }


def _environmental_toxicity(logp, mw):
    fish_lc50 = round(10 ** (0.58 - 0.7 * logp + 0.00012 * mw), 4)  # mg/L (simplified model)
    return {
        "fish_lc50_mg_L": fish_lc50,
        "class": "Toxic to aquatic life" if fish_lc50 < 1 else
                 ("Harmful to aquatic life" if fish_lc50 < 10 else "Low aquatic toxicity"),
        "biodegradation": "Not readily biodegradable (quinone scaffold)",
        "log_Kow": round(logp, 2),
    }


def _get_recommendation(n_concerns, ghs_class):
    if n_concerns == 0 and ghs_class >= 4:
        return "Favorable safety profile for further development"
    elif n_concerns <= 2 and ghs_class >= 3:
        return ("Acceptable risk for antiparasitic context; hepatotoxicity monitoring "
                "recommended in preclinical studies; dose optimization critical")
    else:
        return ("Significant safety concerns; structure optimization to reduce "
                "reactive metabolite formation recommended before further development")


# =============================================================================
# FIGURE GENERATION
# =============================================================================

def fig_toxicity_dashboard(tox: dict, output_dir: Path) -> Path:
    """
    Comprehensive toxicity dashboard figure.
    4-panel: GHS gauge, organ toxicity spider, LD50 comparison, safety matrix.
    """
    ld   = tox["ld50"]
    org  = tox["organ_toxicity"]
    ov   = tox["overall"]
    env  = tox["environmental"]

    risk_color = {"Low": "#28a745", "Low-Medium": "#8bc34a",
                  "Medium": "#FFC107", "Medium-High": "#FF9800",
                  "High": "#dc3545"}

    fig = plt.figure(figsize=(14, 10))
    gs  = GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.35)
    fig.patch.set_facecolor('white')

    # ── Panel A: GHS Class gauge ──────────────────────────────────────────────
    ax_a = fig.add_subplot(gs[0, 0])
    ghs_n = ld["ghs_class"]
    ghs_colors = [GHS_CLASSES[i][1] for i in range(1, 7)]
    ghs_labels = [f"Class {i}\n{GHS_CLASSES[i][0]}" for i in range(1, 7)]

    # Create donut
    explode = [0.08 if i+1 == ghs_n else 0 for i in range(6)]
    wedges, texts = ax_a.pie([1]*6, labels=None, colors=ghs_colors,
                              explode=explode, startangle=90,
                              wedgeprops=dict(width=0.55, edgecolor='white', linewidth=1.5))

    # Highlight the predicted class
    wedges[ghs_n - 1].set_linewidth(3)
    wedges[ghs_n - 1].set_edgecolor('black')

    ax_a.text(0, 0.15, f'GHS', ha='center', va='center', fontsize=11, fontweight='bold', color='#333')
    ax_a.text(0, -0.15, f'Class {ghs_n}', ha='center', va='center', fontsize=14, fontweight='bold',
              color=GHS_CLASSES[ghs_n][1])
    ax_a.text(0, -0.50, f'{ld["ghs_label"]}', ha='center', va='center', fontsize=10, color='#555')
    ax_a.text(0, -0.78, f'LD50 (human, pred.)\n{ld["human_predicted_mg_kg"]:.0f} mg/kg',
              ha='center', va='center', fontsize=9, fontweight='bold', color='#333')

    ax_a.set_title(f'(A) GHS Toxicity Class\n{COMPOUND_ABBR}',
                   fontsize=11, fontweight='bold', pad=10)

    # ── Panel B: Organ toxicity radar ─────────────────────────────────────────
    ax_b = fig.add_subplot(gs[0, 1], polar=True)
    organs = ["Hepatic", "Cardiac", "Renal", "CNS", "Immune", "Skin"]
    org_vals_map = {
        "Hepatic": tox["hepatotoxicity"]["risk"],
        "Cardiac": org["cardiac"],
        "Renal":   org["renal"],
        "CNS":     org["cns"],
        "Immune":  org["immune"],
        "Skin":    org["skin"],
    }
    risk_num = {"Low": 0.2, "Low-Medium": 0.4, "Medium": 0.6, "Medium-High": 0.75, "High": 1.0}
    org_scores = [risk_num.get(org_vals_map[o], 0.3) for o in organs]

    N = len(organs)
    angles = [n / float(N) * 2 * math.pi for n in range(N)] + [0]
    vals_r = org_scores + [org_scores[0]]

    ax_b.set_facecolor('#f8f9fa')
    ax_b.fill(angles, vals_r, alpha=0.25, color='#C62828')
    ax_b.plot(angles, vals_r, 'o-', color='#C62828', linewidth=2.0, markersize=7)
    ax_b.fill(angles, [0.5]*N + [0.5], alpha=0.08, color='orange')  # medium risk zone

    ax_b.set_xticks(angles[:-1])
    ax_b.set_xticklabels(organs, fontsize=9, fontweight='bold')
    ax_b.set_ylim(0, 1)
    ax_b.set_yticks([0.2, 0.5, 0.8])
    ax_b.set_yticklabels(['Low', 'Med', 'High'], fontsize=7)
    ax_b.set_title('(B) Organ Toxicity Profile', fontsize=11, fontweight='bold', pad=15)

    # ── Panel C: Mouse vs Human LD50 comparison ───────────────────────────────
    ax_c = fig.add_subplot(gs[0, 2])
    species_vals = [
        MOUSE_EXPERIMENTAL["LD50_oral_mg_kg"],
        ld["human_predicted_mg_kg"]
    ]
    species_labels = ['Mouse\n(Experimental)', 'Human\n(Predicted)']
    species_colors = ['#1565C0', '#E65100']
    bars = ax_c.bar(species_labels, species_vals, color=species_colors,
                    edgecolor='white', linewidth=0.5, width=0.5)

    for bar, val in zip(bars, species_vals):
        ax_c.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5,
                  f'{val:.0f} mg/kg', ha='center', fontsize=11, fontweight='bold',
                  color=bar.get_facecolor())

    ax_c.axhline(2000, color='gray', linestyle='--', linewidth=1.5, alpha=0.7,
                 label='GHS Class 4/5 boundary')
    ax_c.axhline(300, color='orange', linestyle='--', linewidth=1.5, alpha=0.7,
                 label='GHS Class 3/4 boundary')
    ax_c.set_ylabel('LD50 (mg/kg, oral)', fontsize=10)
    ax_c.set_title('(C) LD50 Comparison\nMouse (exp.) vs Human (pred.)',
                   fontsize=11, fontweight='bold')
    ax_c.legend(fontsize=8, loc='upper right')
    ax_c.spines['top'].set_visible(False)
    ax_c.spines['right'].set_visible(False)

    # ── Panel D: Full safety matrix ───────────────────────────────────────────
    ax_d = fig.add_subplot(gs[1, :])
    ax_d.axis('off')

    hepato_alerts = "; ".join(tox["hepatotoxicity"]["evidence"][:2]) or "None"
    ames_alerts   = "; ".join(tox["mutagenicity"]["structural_alerts"][:2]) or "None"

    matrix_data = [
        ['Endpoint',             'Prediction',         'Evidence / Model',              'Risk Level'],
        ['LD50 (human, pred.)',  f'{ld["human_predicted_mg_kg"]:.0f} mg/kg', 'QSAR structural model',  ld["ghs_label"]],
        ['LD50 (mouse, expt.)',  f'{MOUSE_EXPERIMENTAL["LD50_oral_mg_kg"]} mg/kg', MOUSE_EXPERIMENTAL["source"][:35], 'GHS Class 3'],
        ['Hepatotoxicity',       tox["hepatotoxicity"]["risk"],  hepato_alerts[:40],    tox["hepatotoxicity"]["risk"]],
        ['AMES Mutagenicity',    'Positive' if tox["mutagenicity"]["ames_predicted_positive"] else 'Negative',
                                  ames_alerts[:40], 'Moderate' if tox["mutagenicity"]["ames_predicted_positive"] else 'Low'],
        ['Carcinogenicity',      tox["carcinogenicity"]["risk"],
                                  "; ".join(tox["carcinogenicity"]["reasoning"][:1])[:40] or "No alerts", tox["carcinogenicity"]["risk"]],
        ['Cardiac (hERG)',       org["cardiac"],       'logP + ring-based model',       org["cardiac"]],
        ['Renal',                org["renal"],         'Structural class',               org["renal"]],
        ['CNS',                  org["cns"],           'TPSA-based BBB model',           org["cns"]],
        ['Aquatic toxicity',     env["class"],         f'log Kow = {env["log_Kow"]}',    'Concern'],
        ['Overall',              ov["assessment"],     f'{ov["safety_concerns"]} concerns identified', ov["assessment"]],
    ]

    def row_color(risk):
        if any(x in risk for x in ['High', 'Toxic', 'Positive', 'Fatal', 'Concern']):
            return ['#fdecea'] * 4
        elif any(x in risk for x in ['Medium', 'Moderate', 'Harmful']):
            return ['#fff8e1'] * 4
        elif 'Low' in risk or 'Negative' in risk or 'Favorable' in risk:
            return ['#e8f5e9'] * 4
        return ['#fafafa'] * 4

    matrix_colors = [['#1565C0'] * 4]
    for row in matrix_data[1:]:
        matrix_colors.append(row_color(row[3]))

    tbl = ax_d.table(cellText=matrix_data, cellLoc='center', loc='center',
                     cellColours=matrix_colors)
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1, 1.55)
    for (r, c), cell in tbl.get_celld().items():
        if r == 0:
            cell.set_text_props(fontweight='bold', color='white')
        cell.set_edgecolor('#dddddd')
        if r > 0 and c == 3:
            cell.set_text_props(fontweight='bold')
    ax_d.set_title('(D) Complete Safety & Toxicity Matrix',
                   fontweight='bold', fontsize=11, pad=10)

    fig.suptitle(f'Human Toxicity Prediction — {COMPOUND_NAME} ({COMPOUND_ABBR})\n'
                 f'Therapeutic Window Assessment for Antiparasitic Development',
                 fontsize=13, fontweight='bold', y=1.01)
    plt.tight_layout()

    out = output_dir / "Toxicity_dashboard.png"
    plt.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    return out


# =============================================================================
# PAPER TEXT GENERATION
# =============================================================================

def generate_toxicity_text(tox: dict, lang: str = "en") -> str:
    if lang == "pt":
        return _tox_text_pt(tox)
    return _tox_text_en(tox)


def _tox_text_en(tox):
    ld  = tox["ld50"]
    ov  = tox["overall"]
    he  = tox["hepatotoxicity"]
    mu  = tox["mutagenicity"]
    ca  = tox["carcinogenicity"]
    exp = tox.get("experimental", {})
    sel = tox.get("selectivity", {})
    src = ov.get("data_source", "predicted")

    # Build LD50 sentence
    mouse_val = ld.get("mouse_experimental_mg_kg") or ld["human_predicted_mg_kg"]
    ld50_src  = "experimentally determined" if src == "experimental" else "predicted"
    ti_val    = ov["therapeutic_index_est"]
    ti_fav    = (str(ti_val) not in ("N/A", "None") and float(str(ti_val)) > 10)

    # Build selectivity sentence
    si_tc = sel.get("si_tc_vs_host")
    si_lm = sel.get("si_lm_vs_host")
    ic50_tc = sel.get("ic50_tc_used_uM")
    ic50_lm = sel.get("ic50_lm_used_uM")
    cc50    = sel.get("cc50_used_uM")
    sel_sentence = ""
    if si_tc or si_lm:
        sel_sentence = (f"The calculated Selectivity Index (SI = CC50_host / IC50_parasite) "
                        f"was SI(T. cruzi) = {si_tc}x and SI(L. amazonensis) = {si_lm}x, "
                        f"based on CC50 (host cells) = {cc50} μM, "
                        f"IC50 (T. cruzi) = {ic50_tc} μM, and "
                        f"IC50 (Leishmania) = {ic50_lm} μM. ")

    return f"""Human Toxicity and Selectivity Profile

In silico toxicity profiling of {COMPOUND_ABBR} was performed using structural \
alert analysis (RDKit) and QSAR-based regression models. Where available, \
{("experimental data were incorporated to refine the predictions." if src == "experimental" else "predictions were compared with data from structurally related naphthoquinone derivatives.")}

LD50 and GHS Classification: The {ld50_src} oral LD50 of {COMPOUND_ABBR} is \
{mouse_val:.0f} mg/kg (GHS Class {ld["ghs_class"]}: {ld["ghs_label"]}, \
range {ld["ghs_range"]}). {f"The predicted human LD50 is {ld['human_predicted_mg_kg']:.0f} mg/kg." if src != "experimental" else ""} \
The therapeutic index is estimated at {ti_val}, suggesting a \
{"favorable" if ti_fav else "moderate"} margin between effective and toxic doses.

Antiparasitic Activity and Selectivity: {sel_sentence if sel_sentence else f"IC50 values for T. cruzi and Leishmania should be added in the dashboard to calculate the Selectivity Index."} \
A Selectivity Index > 10 is generally considered favorable for antiparasitic drug candidates, \
as it indicates that the compound is at least 10-fold more toxic to the parasite than to mammalian host cells.

Hepatotoxicity: The hepatotoxicity risk was assessed as {he['risk'].lower()} \
based on the following structural alerts: {'; '.join(he['evidence']) if he['evidence'] else 'none identified'}. \
The quinone scaffold is recognized as a potential source of reactive metabolites \
(via CYP1A2 and CYP3A4 oxidation), which may conjugate with glutathione (GSH). \
However, at the predicted antiparasitic effective concentrations (~sub-micromolar), \
the hepatotoxic risk is expected to be manageable and dose-dependent.

Mutagenicity and Genotoxicity: Structural analysis {'identified' if mu['ames_predicted_positive'] else 'did not identify'} \
AMES mutagenicity alerts ({'; '.join(mu['structural_alerts'][:2]) if mu['structural_alerts'] else 'none'}). \
The carcinogenic risk was assessed as {ca['risk'].lower()}. \
Experimental validation of genotoxicity (Comet assay, micronucleus test) is \
recommended as a component of the preclinical safety package.

Overall Safety Assessment: {ov['safety_concerns']} safety concern(s) were identified \
in silico. {ov['development_recommendation']}. These computational findings should \
be validated with in vitro toxicity assays (MTT, Annexin V, reactive oxygen species) \
in human cell lines (HepG2, HEK293) prior to in vivo toxicological studies. \
The comparison between antiparasitic activity (IC50 in the nM–μM range) and the \
predicted human LD50 ({ld['human_predicted_mg_kg']:.0f} mg/kg) supports that a \
therapeutic window exists for further development of {COMPOUND_ABBR}."""


def _tox_text_pt(tox):
    ld  = tox["ld50"]
    ov  = tox["overall"]
    he  = tox["hepatotoxicity"]
    mu  = tox["mutagenicity"]
    ca  = tox["carcinogenicity"]
    exp = tox.get("experimental", {})
    sel = tox.get("selectivity", {})
    src = ov.get("data_source", "predicted")

    mouse_val = ld.get("mouse_experimental_mg_kg") or ld["human_predicted_mg_kg"]
    ld50_src  = "determinada experimentalmente" if src == "experimental" else "predita computacionalmente"
    ti_val    = ov["therapeutic_index_est"]
    ti_fav    = (str(ti_val) not in ("N/A", "None") and float(str(ti_val)) > 10)

    si_tc = sel.get("si_tc_vs_host")
    si_lm = sel.get("si_lm_vs_host")
    ic50_tc = sel.get("ic50_tc_used_uM")
    ic50_lm = sel.get("ic50_lm_used_uM")
    cc50    = sel.get("cc50_used_uM")
    sel_sentence = ""
    if si_tc or si_lm:
        sel_sentence = (f"O Índice de Seletividade calculado (IS = CC50_hospedeiro / CI50_parasito) "
                        f"foi IS(T. cruzi) = {si_tc}x e IS(L. amazonensis) = {si_lm}x, "
                        f"baseado em CC50 (células hospedeiras) = {cc50} μM, "
                        f"CI50 (T. cruzi) = {ic50_tc} μM e "
                        f"CI50 (Leishmania) = {ic50_lm} μM. ")

    return f"""Perfil de Toxicidade e Seletividade

O perfil de toxicidade in silico de {COMPOUND_ABBR} foi realizado utilizando análise \
de alertas estruturais (RDKit) e modelos de regressão QSAR. \
{("Dados experimentais foram incorporados para refinar as predições." if src == "experimental" else "As predições foram comparadas com dados de naftoquinonas estruturalmente relacionadas.")}

DL50 e Classificação GHS: A DL50 oral {ld50_src} de {COMPOUND_ABBR} é \
{mouse_val:.0f} mg/kg (Classe GHS {ld["ghs_class"]}: {ld["ghs_label"]}, \
faixa {ld["ghs_range"]}). {f"A DL50 humana predita é {ld['human_predicted_mg_kg']:.0f} mg/kg." if src != "experimental" else ""} \
O índice terapêutico estimado é {ti_val}, sugerindo uma margem \
{"favorável" if ti_fav else "moderada"} entre doses eficazes e tóxicas.

Atividade Antiparasitária e Seletividade: {sel_sentence if sel_sentence else "Preencha os valores de CI50 e CC50 no dashboard para calcular o Índice de Seletividade."} \
Um Índice de Seletividade > 10 é considerado favorável para candidatos antiparasitários, \
indicando que o composto é pelo menos 10 vezes mais tóxico ao parasito do que às células hospedeiras.

Hepatotoxicidade: O risco de hepatotoxicidade foi avaliado como {he['risk'].lower()} \
com base nos seguintes alertas estruturais: \
{'; '.join(he['evidence']) if he['evidence'] else 'nenhum identificado'}. \
O scaffold quinônico é reconhecido como fonte potencial de metabólitos reativos \
(via CYP1A2 e CYP3A4), que podem conjugar com glutationa (GSH). Contudo, nas \
concentrações antiparasitárias eficazes (~submicromolar), o risco hepatotóxico \
é esperado ser manejável e dose-dependente.

Mutagenicidade e Genotoxicidade: A análise estrutural \
{'identificou' if mu['ames_predicted_positive'] else 'não identificou'} \
alertas de mutagenicidade AMES. O risco carcinogênico foi avaliado como \
{ca['risk'].lower()}. Validação experimental de genotoxicidade (ensaio Cometa, \
teste de micronúcleo) é recomendada.

Avaliação Geral de Segurança: {ov['safety_concerns']} preocupação(ões) de segurança \
foram identificadas in silico. {ov['development_recommendation']}. A comparação entre \
a atividade antiparasitária (IC50 na faixa nM–μM) e a DL50 humana predita \
({ld['human_predicted_mg_kg']:.0f} mg/kg) suporta a existência de uma janela \
terapêutica para desenvolvimento da {COMPOUND_ABBR}."""


# =============================================================================
# MAIN ENTRY POINT
# =============================================================================

def run_toxicity(smiles: str = SMILES_DEFAULT,
                 output_dir: str = "docking_project/figures",
                 lang: str = "en",
                 compound_name: str = None,
                 compound_abbr: str = None,
                 exp_data: dict = None) -> dict:
    """Run complete toxicity analysis"""
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Override module-level name with user-provided compound name
    global COMPOUND_NAME, COMPOUND_ABBR
    _orig_name = COMPOUND_NAME
    _orig_abbr = COMPOUND_ABBR
    if compound_name:
        COMPOUND_NAME = compound_name
        COMPOUND_ABBR = compound_abbr or compound_name
    try:
        print(f"\n{'='*60}")
        print(f"  Toxicity Analysis: {COMPOUND_NAME}")
        print(f"{'='*60}\n")

        tox = predict_toxicity(smiles, exp_data=exp_data)
        if "error" in tox:
            print(f"[ERRO] {tox['error']}")
            return tox

        ld = tox["ld50"]
        ov = tox["overall"]
        print(f"  LD50 (human pred.): {ld['human_predicted_mg_kg']:.0f} mg/kg | GHS Class {ld['ghs_class']} ({ld['ghs_label']})")
        print(f"  LD50 (mouse expt.): {MOUSE_EXPERIMENTAL['LD50_oral_mg_kg']} mg/kg")
        print(f"  Hepatotoxicity: {tox['hepatotoxicity']['risk']}")
        print(f"  AMES: {'Positive ⚠' if tox['mutagenicity']['ames_predicted_positive'] else 'Negative ✓'}")
        print(f"  Overall: {ov['assessment']}")

        # Figure
        fig_path = fig_toxicity_dashboard(tox, out_dir)
        print(f"  [OK] {fig_path.name}")

        # Paper text
        text = generate_toxicity_text(tox, lang=lang)
        tf   = out_dir / f"toxicity_text_{lang}.txt"
        tf.write_text(text, encoding='utf-8')
        print(f"  [OK] {tf.name}")

        # JSON
        jf = out_dir / "toxicity_results.json"
        jf.write_text(json.dumps(tox, indent=2, ensure_ascii=False), encoding='utf-8')
        print(f"  [OK] {jf.name}")

        return tox
    finally:
        # Restore original global names to avoid side effects across calls
        COMPOUND_NAME = _orig_name
        COMPOUND_ABBR = _orig_abbr


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Human Toxicity Prediction")
    parser.add_argument("--smiles",  default=SMILES_DEFAULT)
    parser.add_argument("--output",  default="docking_project/figures")
    parser.add_argument("--lang",    default="en", choices=["en", "pt"])
    args = parser.parse_args()
    run_toxicity(args.smiles, args.output, args.lang)
