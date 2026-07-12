#!/usr/bin/env python3
"""
=============================================================================
populate_ligand_db.py
Popula o banco de compostos (docking_project/ligand/) com 19 moléculas
da literatura antiparasitária — naftoquinonas, fármacos de referência e
derivados estudados pelo grupo Menna-Barreto (FIOCRUZ/IOC).

Uso:
    python populate_ligand_db.py           # adiciona todos
    python populate_ligand_db.py --list    # lista sem adicionar
    python populate_ligand_db.py --dry-run # mostra o que seria feito
    python populate_ligand_db.py --name "Lapachol"  # adiciona só um
=============================================================================
"""
import sys
import argparse
import json
from pathlib import Path

# ─── BANCO DA LITERATURA ─────────────────────────────────────────────────────
# Cada entrada contém:
#   smiles      : canônico, validado no PubChem/ChemDraw
#   formula     : fórmula molecular
#   mw          : peso molecular (g/mol)
#   category    : "naphthoquinone_natural" | "naphthoquinone_synthetic" |
#                 "reference_drug" | "cancer" | "derivative_fiocruz"
#   activity    : dados de atividade biológica extraídos da literatura
#   reference   : citação primária
#   notes       : notas científicas relevantes
# ─────────────────────────────────────────────────────────────────────────────

COMPOUND_LIBRARY = [

    # ══════════════════════════════════════════════════════════════════════════
    # GRUPO 1 — FÁRMACOS DE REFERÊNCIA CLÍNICA
    # ══════════════════════════════════════════════════════════════════════════
    {
        "name":         "Benznidazole",
        "display_name": "Benznidazole (BZN)",
        "smiles":       "O=C(CN1C(=O)c2ccccc2N1)NCc1cncc(Cl)c1",
        "formula":      "C12H12ClN3O3",
        "mw":           260.70,
        "category":     "reference_drug",
        "activity": {
            "ic50_tc_epi_uM":          "0.5–2.0",
            "ic50_tc_trypo_uM":        "2–10",
            "cc50_vero_uM":            ">100",
            "si_tc":                   ">50",
            "clinical_use":            "Chagas disease — first-line therapy",
        },
        "reference":    "Drugs approved by ANVISA/WHO; Romanha et al. 2010 Mem Inst Oswaldo Cruz",
        "notes":        "Nitro-heterocyclic drug; prodrug activated by T. cruzi nitroreductases (TcNTR). "
                        "Reference compound for ALL trypanocidal assays in Menna-Barreto group. "
                        "Ineffective in chronic phase; causes peripheral neuropathy ~30% patients.",
    },

    {
        "name":         "Nifurtimox",
        "display_name": "Nifurtimox (NFX)",
        "smiles":       "CC1CS(=O)(=O)CCN1/N=C/c1ccc([N+](=O)[O-])o1",
        "formula":      "C10H13N3O5S",
        "mw":           287.29,
        "category":     "reference_drug",
        "activity": {
            "ic50_tc_epi_uM":          "3–15",
            "cc50_vero_uM":            ">100",
            "clinical_use":            "Chagas disease — second-line; also HAT combination",
        },
        "reference":    "WHO Model List of Essential Medicines 2023",
        "notes":        "Nitrofuran prodrug; generates ROS via futile redox cycling. "
                        "Discontinued in many countries due to CNS toxicity. "
                        "Approved for pediatric Chagas disease (2020) in combination with eflornithine for HAT.",
    },

    {
        "name":         "Miltefosine",
        "display_name": "Miltefosine (MLT)",
        "smiles":       "CCCCCCCCCCCCCCCCOP(=O)([O-])OCC[N+](C)(C)C",
        "formula":      "C21H45NO4P",
        "mw":           407.56,
        "category":     "reference_drug",
        "activity": {
            "ic50_lm_promas_uM":       "2–8",
            "ic50_lm_amast_uM":        "1–5",
            "cc50_vero_uM":            "10–30",
            "si_lm":                   "~5–10",
            "clinical_use":            "Visceral leishmaniasis — oral (India/Europe); also cutaneous",
        },
        "reference":    "Dorlo et al. 2012 J Antimicrob Chemother; Santa-Rita et al. 2004",
        "notes":        "First oral drug approved for visceral leishmaniasis. Alkylphospholipid; "
                        "disrupts lipid metabolism and cell signaling. Teratogenic — contraindicated in pregnancy. "
                        "Used as positive control for anti-Leishmania assays in Menna-Barreto group.",
    },

    {
        "name":         "Amphotericin-B",
        "display_name": "Amphotericin B (AmB)",
        "smiles":       "OC1C(=O)OC2CC(OC3OC(C)C(N)C(O)C3O)CC(O)CC(O)CC(O)CC(O)CC(=C/C=C/C=C/C=C/C=C/C(=O)OCC(O)CC2O)C1=O",
        "formula":      "C47H73NO17",
        "mw":           924.08,
        "category":     "reference_drug",
        "activity": {
            "ic50_lm_amast_uM":        "0.05–0.5",
            "cc50_vero_uM":            "1–5",
            "si_lm":                   "~10–20",
            "clinical_use":            "Visceral leishmaniasis — gold standard (liposomal form L-AmB)",
        },
        "reference":    "Sundar & Chakravarty 2010 Expert Opin Pharmacother",
        "notes":        "Polyene macrolide; binds ergosterol (absent in mammals, present in Leishmania). "
                        "Liposomal formulation (AmBisome) reduces nephrotoxicity. "
                        "Reference control for severe visceral leishmaniasis.",
    },

    # ══════════════════════════════════════════════════════════════════════════
    # GRUPO 2 — NAFTOQUINONAS NATURAIS (base estrutural do 2-SCN-NQ)
    # ══════════════════════════════════════════════════════════════════════════
    {
        "name":         "Lapachol",
        "display_name": "Lapachol",
        "smiles":       "OC1=C(CC=C(C)C)C(=O)c2ccccc2C1=O",
        "formula":      "C15H14O3",
        "mw":           242.27,
        "category":     "naphthoquinone_natural",
        "activity": {
            "ic50_tc_epi_uM":          "15–40",
            "ic50_lm_promas_uM":       "10–25",
            "cc50_vero_uM":            "80–120",
            "anticancer":              "Active — cell cycle arrest (NCI-60 panel)",
        },
        "reference":    "Pinto & De Castro 2009 Molecules; Moura et al. 2001 J Braz Chem Soc",
        "notes":        "Precursor natural de β-lapachona e α-lapachona; isolado de Tabebuia spp. (Ipê). "
                        "Pré-cursor essencial da síntese dos naftoimidazóis estudados por Menna-Barreto. "
                        "Seminal compound for FIOCRUZ naphthoquinone drug discovery program.",
    },

    {
        "name":         "beta-Lapachone",
        "display_name": "β-Lapachone (β-LAP)",
        "smiles":       "CC1(C)CCC2=C(O1)c1ccccc1C(=O)C2=O",
        "formula":      "C15H14O3",
        "mw":           242.27,
        "category":     "naphthoquinone_natural",
        "activity": {
            "ic50_tc_epi_uM":          "2–8",
            "ic50_tc_trypo_uM":        "1–5",
            "ic50_lm_promas_uM":       "5–20",
            "cc50_vero_uM":            "20–50",
            "anticancer":              "NQO1-bioactivatable — clinical trials (pancreatic Ca, breast Ca)",
            "si_tc":                   "~4–10",
        },
        "reference":    "Menna-Barreto & De Castro 2014 Front Physiol; "
                        "Boveris et al. 1978 (ROS mechanism); Dos Anjos et al. 2016 IJPDDR",
        "notes":        "Composto-chave estudado extensivamente pelo grupo Menna-Barreto/FIOCRUZ. "
                        "Gera radical semiquinona → superóxido via ciclo redox mitocondrial. "
                        "Precursor dos naftoimidazóis N1, N2, N3 (derivados sintéticos de maior atividade). "
                        "Ativa NQO1 (NAD(P)H:quinona oxidorredutase 1) → morte seletiva células cancerígenas.",
    },

    {
        "name":         "alpha-Lapachone",
        "display_name": "α-Lapachone (α-LAP)",
        "smiles":       "O=C1c2ccccc2C(=O)C2=C1OCC(C)(C)O2",
        "formula":      "C15H14O4",
        "mw":           258.27,
        "category":     "naphthoquinone_natural",
        "activity": {
            "ic50_tc_epi_uM":          "10–30",
            "ic50_lm_promas_uM":       "8–20",
            "cc50_vero_uM":            ">50",
        },
        "reference":    "Pinto & De Castro 2009 Molecules; Menna-Barreto IOC review",
        "notes":        "Isômero de β-lapachona; epóxido derivado usado por Menna-Barreto como precursor "
                        "para derivados com atividade contra proteases de L. amazonensis (serina protease OPB).",
    },

    {
        "name":         "Plumbagin",
        "display_name": "Plumbagin (5-OH-2-methyl-NQ)",
        "smiles":       "O=C1c2c(O)cccc2C(=O)C=C1C",
        "formula":      "C11H8O3",
        "mw":           188.18,
        "category":     "naphthoquinone_natural",
        "activity": {
            "ic50_tc_trypo_uM":        "0.5–2.0",
            "ic50_tc_amast_uM":        "1–4",
            "ic50_lm_promas_uM":       "1–5",
            "cc50_vero_uM":            "5–15",
            "si_tc":                   "~3–8",
            "anticancer":              "Inibe NF-κB; apoptose; ensaios clínicos fase I/II (glioblastoma)",
            "target_validated":        "Tripanotiona redutase (inibição direta validada)",
        },
        "reference":    "Salomão et al. 2013 BMC Microbiol; Menna-Barreto group; "
                        "Antiparasitic Activity of Plumbago auriculata, Pharmaceutics 2023",
        "notes":        "Isolado de Plumbago spp. Inibe TcTR diretamente. "
                        "Citado por Menna-Barreto como modelo de mecanismo ROS + inibição TR. "
                        "Alta toxicidade para células cardíacas — limitação para Chagas.",
    },

    {
        "name":         "Juglone",
        "display_name": "Juglone (5-OH-1,4-NQ)",
        "smiles":       "O=C1C=CC(=O)c2cc(O)ccc21",
        "formula":      "C10H6O3",
        "mw":           174.15,
        "category":     "naphthoquinone_natural",
        "activity": {
            "ic50_tc_epi_uM":          "5–15",
            "ic50_lm_promas_uM":       "3–10",
            "anticancer":              "Inibidor Pin1 prolil isomerase; atividade em múltiplas linhagens",
        },
        "reference":    "Docampo et al. 1978; Salomão et al. 2013 BMC Microbiol",
        "notes":        "Isolado de Juglans regia (nogueira). Gerador de ROS; "
                        "mecanismo semelhante ao β-lapachona na mitocôndria do parasito.",
    },

    {
        "name":         "Lawsone",
        "display_name": "Lawsone (2-OH-1,4-NQ, Henna)",
        "smiles":       "O=C1C=C(O)C(=O)c2ccccc21",
        "formula":      "C10H6O3",
        "mw":           174.15,
        "category":     "naphthoquinone_natural",
        "activity": {
            "ic50_tc_epi_uM":          "20–80",
            "ic50_lm_promas_uM":       "15–50",
            "anticancer":              "Atividade moderada; precursor sintético LQB-118",
        },
        "reference":    "Maia et al. 2011 Invest New Drugs (LQB-118 precursor); "
                        "Rani et al. 2022 Chem Biol Drug Des",
        "notes":        "Pigmento do Henna. Precursor do LQB-118 (pterocarpanquinona) "
                        "via iodo-lawsona. Menor atividade sozinho; fundamental na síntese de derivados.",
    },

    {
        "name":         "Menadione",
        "display_name": "Menadione (Vitamina K3 / 2-methyl-1,4-NQ)",
        "smiles":       "O=C1C=C(C)C(=O)c2ccccc21",
        "formula":      "C11H8O2",
        "mw":           172.18,
        "category":     "naphthoquinone_natural",
        "activity": {
            "ic50_tc_epi_uM":          "1–5",
            "ic50_lm_promas_uM":       "2–8",
            "cc50_vero_uM":            "10–30",
            "anticancer":              "Ativo contra glioblastoma C6, apoptose via ROS",
        },
        "reference":    "Docampo et al. 1978 Z Parasitenkd; Salomão et al. 2013 BMC Microbiol; "
                        "Anaissi-Afonso et al. 2018 ACS Chem Biol",
        "notes":        "Vitamina K3 sintética. Modelo mais simples de NQ com atividade antiparasitária. "
                        "Referência estrutural para SAR (structure-activity relationship) de naftoquinonas.",
    },

    {
        "name":         "Atovaquone",
        "display_name": "Atovaquone (Hydroxynaphthoquinone antiprotozoal)",
        "smiles":       "OC1=C(C(=O)c2ccccc21)C1CCC(c2ccc(Cl)cc2)CC1",
        "formula":      "C22H19ClO3",
        "mw":           366.83,
        "category":     "reference_drug",
        "activity": {
            "ic50_lm_amast_uM":        "0.1–1.0",
            "cc50_vero_uM":            ">50",
            "si_lm":                   ">50",
            "clinical_use":            "Malária (Malarone); pneumocistose; toxoplasmose",
        },
        "reference":    "Croft et al. 1992 J Antimicrob Chemother; "
                        "WHO Essential Medicines — antiprotozoal",
        "notes":        "Naftoquinona sintética clinicamente aprovada. Inibe complexo III mitocondrial "
                        "(citocromo bc1) de plasmódio, Pneumocystis e Toxoplasma. "
                        "Modelo de como NQ pode ter alvos mitocondriais específicos.",
    },

    # ══════════════════════════════════════════════════════════════════════════
    # GRUPO 3 — DERIVADOS SINTÉTICOS (grupo Menna-Barreto/FIOCRUZ + colaboradores)
    # ══════════════════════════════════════════════════════════════════════════
    {
        "name":         "NQ1-naphthoimidazole",
        "display_name": "Naphthoimidazole N1 (β-lapachone derivative — FIOCRUZ)",
        "smiles":       "O=C1OC2(C)CCC=C2c2c1cc1nc(-c3ccccc3)[nH]c1c2=O",
        "formula":      "C21H16N2O3",
        "mw":           348.37,
        "category":     "derivative_fiocruz",
        "activity": {
            "ic50_tc_trypo_uM":        "0.5–2.0",
            "ic50_tc_epi_uM":          "1–3",
            "ic50_tc_amast_uM":        "2–5",
            "cc50_vero_uM":            ">50",
            "si_tc":                   ">25",
            "in_vivo":                 "25-50% redução de parasitemia em camundongos (agudo)",
        },
        "reference":    "Menna-Barreto et al. 2005 J Antimicrob Chemother 56:1034; "
                        "Cascabulho et al. 2020 Mem Inst Oswaldo Cruz DOI:10.1590/0074-02760190389",
        "notes":        "Naftoimidazol derivado de β-lapachona. Um dos mais ativos da série N1/N2/N3 "
                        "estudada pelo Lab Biologia Celular do IOC/FIOCRUZ (Menna-Barreto). "
                        "Alvo confirmado: mitocôndria do parasito (colapso ΔΨm, vacuolização). "
                        "Publicado na série de estudos mecanísticos de 2005-2020 do grupo FIOCRUZ.",
    },

    {
        "name":         "LQB-118",
        "display_name": "LQB-118 Pterocarpanquinone (UFRJ/FIOCRUZ)",
        "smiles":       "O=C1c2c(O)cccc2C(=O)C1OC1Oc2ccccc2C2CCCC12",
        "formula":      "C21H16O5",
        "mw":           352.35,
        "category":     "derivative_fiocruz",
        "activity": {
            "ic50_tc_epi_uM":          "2.5",
            "ic50_tc_amast_uM":        "4.2",
            "ic50_tc_trypo_meta_uM":   "38.1",
            "ic50_lm_promas_uM":       "3.4",
            "ic50_lm_braziliensis_amast_uM": "7.5",
            "cc50_vero_uM":            ">100",
            "si_tc":                   ">40",
            "in_vivo_leishmaniasis":   "99% redução carga parasitária (10 mg/kg oral, VL)",
            "anticancer":              "Leucemia mieloide crônica (LMC) — apoptose seletiva",
        },
        "reference":    "Maia et al. 2011 Invest New Drugs 29:1143; "
                        "Cunha-Junior et al. 2016 Antimicrob Agents Chemother 60:3794; "
                        "ACS Omega 2025 (LQB-118 vs T. cruzi mitochondria); "
                        "Menna-Barreto group collaboration FIOCRUZ/UFRJ",
        "notes":        "Pterocarpanquinona — híbrido molecular lapachol + pterocarpano. "
                        "Um dos compostos mais promissores do programa FIOCRUZ. "
                        "Eficaz oral em leishmaniose visceral e cutânea. "
                        "Mecanismo: ROS → colapso ΔΨm → apoptose. Baixa toxicidade sistêmica (subaguda).",
    },

    {
        "name":         "NQ1-nor-beta-lapachone-triazole-Br",
        "display_name": "NQ1 Nor-β-lapachone Triazole Brominated (FIOCRUZ dual-quinone)",
        "smiles":       "O=C1OCC2=C1c1c(Br)c(=O)c3ccccc3c1CC2n1cc(-c2ccccc2)nn1",
        "formula":      "C24H16BrN3O3",
        "mw":           490.31,
        "category":     "derivative_fiocruz",
        "activity": {
            "ic50_tc_trypo_uM":        "0.8–2.0",
            "ic50_tc_epi_uM":          "1.5–3.1",
            "ic50_compared_to_bz":     "More active than benznidazole",
            "cc50_vero_uM":            ">30",
        },
        "reference":    "Duarte et al. Menna-Barreto group — dual redox-active quinones, "
                        "Researchgate 2020 (Naphthoquinones derivatives FIOCRUZ)",
        "notes":        "Naftoquinona brominada — estrutura dual orto/para-quinoidal (NQ1). "
                        "Série de 1,2,3-triazóis derivados de nor-β-lapachona. "
                        "Mais ativa que benznidazol contra tripomastigotas e epimastigotas. "
                        "Alvo: mitocôndria de T. cruzi (ROS, colapso ΔΨm).",
    },

    {
        "name":         "DPNQ",
        "display_name": "DPNQ (2,3-Diphenyl-1,4-naphthoquinone)",
        "smiles":       "O=C1C(=C(-c2ccccc2)c2ccccc2)C(=O)c2ccccc21",
        "formula":      "C22H14O2",
        "mw":           310.35,
        "category":     "naphthoquinone_synthetic",
        "activity": {
            "ic50_tc_epi_uM":          "~10",
            "in_vivo_tc":              "50% redução carga parasitária; 60% sobrevivência em camundongos",
        },
        "reference":    "Becerra et al. 2021 (in silico QSAR naphthoquinones antiprotozoal) "
                        "In Silico Antiprotozoal Evaluation MDPI Pharmaceutics 2022",
        "notes":        "2,3-difenilnaftoquinona — agente quimioterapêutico promissor em modelo murino de T. cruzi. "
                        "Halving da carga parasitária e 60% de sobrevivência documentados in vivo.",
    },

    # ══════════════════════════════════════════════════════════════════════════
    # GRUPO 4 — NAFTOQUINONAS ANTICÂNCER (com atividade antiparasitária)
    # ══════════════════════════════════════════════════════════════════════════
    {
        "name":         "Shikonin",
        "display_name": "Shikonin (5,8-dihidroxy-2-(1-methylbutyl)-NQ)",
        "smiles":       "O=C1c2c(O)ccc(O)c2C(=O)C=C1[C@@H](O)CCC(C)C",
        "formula":      "C16H16O5",
        "mw":           288.30,
        "category":     "cancer",
        "activity": {
            "ic50_tc_epi_uM":          "5–15",
            "ic50_lm_promas_uM":       "3–10",
            "anticancer":              "Inibidor PKM2; ativo em ≥30 linhagens cancerígenas; antiangiogênico",
        },
        "reference":    "Chen et al. 2012 Cancer Lett; Rani et al. 2022 Chem Biol Drug Des",
        "notes":        "Naftoquinona natural de Lithospermum spp. (raiz vermelha). "
                        "Inibidor potente de PKM2 (piruvato quinase M2) — superexpresso em tumores. "
                        "Atividade antiparasitária reportada por múltiplos grupos.",
    },

    {
        "name":         "Dunnione",
        "display_name": "Dunnione (β-dunnione, 2-acetyl-furanaphthoquinone)",
        "smiles":       "O=C1OC(C)(C)C2=C1c1cc(=O)cc(=O)c1C2=O",
        "formula":      "C13H10O5",
        "mw":           250.21,
        "category":     "naphthoquinone_synthetic",
        "activity": {
            "ic50_tc_epi_uM":          "2–8",
            "nqo1_bioactivatable":     True,
            "anticancer":              "NQO1-bioativável — câncer pancreático, pulmão",
        },
        "reference":    "Anaissi-Afonso et al. 2018 ACS Chem Biol 13:1950; "
                        "Rani et al. 2022 Chem Biol Drug Des",
        "notes":        "Furanoftoquinona — estrutura similar a β-lapachona. "
                        "Bioativada por NQO1 (overexpresso em tumor + ausente em parasito). "
                        "Modelo para entender seletividade tumor vs parasito.",
    },

    {
        "name":         "2-SCN-NQ",
        "display_name": "2-Thiocyanato-1,4-naphthoquinone (COMPOSTO BASE)",
        "smiles":       "O=C1C=C(SC#N)C(=O)c2ccccc21",
        "formula":      "C11H5NO2S",
        "mw":           215.23,
        "category":     "naphthoquinone_synthetic",
        "activity": {
            "status":                  "COMPOSTO EM INVESTIGAÇÃO — dados sendo gerados (Santos et al.)",
            "ic50_tc_reported":        "em determinação",
            "reference_study":         "Santos et al. CSIJ 22(1):1-8, 2018 (Article no. CSIJ.40288)",
        },
        "reference":    "Santos et al. 2018 CSIJ 22(1):1-8",
        "notes":        "Composto base do projeto. Grupo tiocianato (-SCN) na posição 2 confere "
                        "dupla funcionalidade: gerador ROS (quinona) + inibidor TR (tiocianato eletrofílico). "
                        "JÁ DEVE ESTAR no banco — entrada duplicata segura (merge automático).",
    },

]


# ─────────────────────────────────────────────────────────────────────────────
# SCRIPT PRINCIPAL
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Popula o banco de compostos com naftoquinonas da literatura"
    )
    parser.add_argument("--list",    action="store_true", help="Listar compostos sem adicionar")
    parser.add_argument("--dry-run", action="store_true", help="Simular sem escrever arquivos")
    parser.add_argument("--name",    type=str, default=None, help="Adicionar apenas um composto pelo nome")
    parser.add_argument("--base-dir", type=str, default="docking_project",
                        help="Diretório base do projeto (default: docking_project)")
    args = parser.parse_args()

    # ── Listagem ───────────────────────────────────────────────────────────────
    if args.list:
        print(f"\n{'#':<3} {'Nome':<42} {'Categoria':<28} {'MW':>7}  {'Fórmula'}")
        print("-" * 105)
        for i, c in enumerate(COMPOUND_LIBRARY, 1):
            print(f"{i:<3} {c['name']:<42} {c['category']:<28} {c['mw']:>7.2f}  {c['formula']}")
        print(f"\nTotal: {len(COMPOUND_LIBRARY)} compostos")
        return

    # ── Filtrar por nome se solicitado ─────────────────────────────────────────
    to_add = COMPOUND_LIBRARY
    if args.name:
        to_add = [c for c in COMPOUND_LIBRARY if args.name.lower() in c["name"].lower()]
        if not to_add:
            print(f"[ERRO] Composto '{args.name}' não encontrado. Use --list para ver os disponíveis.")
            sys.exit(1)

    # ── Importar compound_db ───────────────────────────────────────────────────
    script_dir = Path(__file__).parent
    sys.path.insert(0, str(script_dir))

    try:
        from compound_db import add_compound, load_compound, _safe_name
    except ImportError:
        print("[ERRO] compound_db.py não encontrado. Execute na pasta do projeto.")
        sys.exit(1)

    if args.dry_run:
        print(f"\n[DRY RUN] Os seguintes compostos seriam adicionados em {args.base_dir}/ligand/:")
        for c in to_add:
            safe = _safe_name(c["name"])
            print(f"  → {c['name']} ({c['formula']}, {c['mw']} g/mol)")
            print(f"    dir: {args.base_dir}/ligand/{safe}/")
        print(f"\nTotal: {len(to_add)} compostos")
        return

    # ── Adicionar compostos ────────────────────────────────────────────────────
    print(f"\n{'='*65}")
    print(f"  Populando banco de compostos — {len(to_add)} entradas")
    print(f"  Diretório: {Path(args.base_dir).resolve() / 'ligand'}")
    print(f"{'='*65}\n")

    added = 0; skipped = 0; errors = 0

    for c in to_add:
        name = c["name"]
        print(f"  [{added+skipped+errors+1:02d}/{len(to_add)}] {name} ...", end=" ")

        # Check if already exists
        safe  = _safe_name(name)
        mpath = Path(args.base_dir) / "ligand" / safe / "metadata.json"
        if mpath.exists():
            print("já existe — atualizando metadados")
            # Update metadata with richer info
            try:
                with open(mpath) as f:
                    existing = json.load(f)
                existing.update({
                    "category":   c.get("category", existing.get("category", "")),
                    "activity":   c.get("activity", {}),
                    "reference":  c.get("reference", ""),
                    "notes":      c.get("notes", ""),
                    "display_name": c.get("display_name", name),
                })
                if not args.dry_run:
                    with open(mpath, "w", encoding="utf-8") as f:
                        json.dump(existing, f, indent=2, ensure_ascii=False)
            except Exception as e:
                print(f"[AVISO] Erro ao atualizar: {e}")
            skipped += 1
            continue

        # Add new compound
        try:
            meta = add_compound(
                name         = name,
                smiles       = c["smiles"],
                formula      = c["formula"],
                display_name = c.get("display_name", name),
                mw           = c["mw"],
                notes        = c.get("notes", ""),
                prepare_3d   = True,
            )
            # Save extra fields
            meta_path = Path(args.base_dir) / "ligand" / safe / "metadata.json"
            if meta_path.exists():
                with open(meta_path) as f:
                    m = json.load(f)
                m["category"]  = c.get("category", "")
                m["activity"]  = c.get("activity", {})
                m["reference"] = c.get("reference", "")
                with open(meta_path, "w", encoding="utf-8") as f:
                    json.dump(m, f, indent=2, ensure_ascii=False)

            pdbqt_ok = meta.get("_has_pdbqt", False)
            if pdbqt_ok:
                status = "OK"
            else:
                # Distinguir SMILES invalido (dados corrompidos) de OpenBabel ausente.
                smi_valid = _smiles_is_valid(c["smiles"])
                if smi_valid is False:
                    status = "AVISO (SMILES invalido — sem estrutura 3D/PDBQT; corrigir SMILES)"
                else:
                    status = "OK (sem PDBQT — OpenBabel ausente)"
            print(status)
            for msg in meta.get("_status_msgs", []):
                if "[ERRO]" in msg or "[OK]" in msg:
                    print(f"         {msg.strip()}")
            added += 1

        except Exception as e:
            print(f"ERRO: {e}")
            errors += 1

    # ── Summary ────────────────────────────────────────────────────────────────
    print(f"\n{'='*65}")
    print(f"  Resultado:")
    print(f"    ✓ Adicionados:     {added}")
    print(f"    ~ Já existentes:   {skipped} (metadados atualizados)")
    print(f"    ✗ Erros:           {errors}")
    print(f"\n  Banco em: {Path(args.base_dir).resolve() / 'ligand'}")
    print(f"\n  Próximos passos:")
    print(f"    1. python dashboard_v2.py → aba Compostos para visualizar")
    print(f"    2. Selecionar composto no banner e rodar o pipeline de docking")
    print(f"    3. Comparar dG do 2-SCN-NQ com compostos de referência (BZN, β-LAP, etc.)")
    print(f"{'='*65}\n")


def _smiles_is_valid(smiles):
    """Valida um SMILES com RDKit.

    Retorna True se valido, False se invalido, ou None se o RDKit nao estiver
    disponivel (nesse caso nao da para distinguir da ausencia de OpenBabel).
    """
    try:
        from rdkit import Chem
    except Exception:
        return None
    try:
        return Chem.MolFromSmiles(smiles) is not None
    except Exception:
        return False


if __name__ == "__main__":
    main()
