# -*- coding: utf-8 -*-
"""
tools/ingest_campaign.py -- Ingestao READ-ONLY da campanha SUB04 (Fase 7.1).

SUB04 = 4 actinas de Trypanosoma cruzi (Dm28c: TcAct1-4) + actina beta humana
(Hs_ACTB, regua) x 4 drogas do citoesqueleto (citocalasina D, jasplaquinolida,
latrunculina A, latrunculina B), a partir do trabalho ja feito e documentado em
Rho_OCRL/in_silico/drogas_citoesqueleto/ -- 21 scripts rodados a mao, narrativa
em .project-memory/, nunca antes conectados a esta plataforma.

Por que SUB04 e nao SUB02: subprojeto/subprojeto_config.json ja reserva SUB02
(otimizacao dual de derivados NQ) e SUB03 (NQ + beta-ciclodextrina) desde
2026-06-19 -- antes deste plano existir. O PLANO.md de origem ja se chama
"SUB04" e pressupoe SUB01/02/03 ocupados. Confirmado com o autor 2026-08-09.

Principio: NUNCA fingir execucao propria. Todo run ingerido leva
provenance.mode="ingested", reproducible=false e original_script apontando pro
script de origem. A plataforma LE os JSONs de onde estao; o Rho_OCRL continua
dono dos dados -- nada e movido, nada e re-rodado. O que este script GRAVA
(quando --write) fica inteiramente dentro de docking_suite_platform/registry/
e docking_suite_platform/docking_project/SUB04/ -- nunca escreve no Rho_OCRL.

dry_run e o DEFAULT (mandado pelo plano, Fase 7.1): sem --write, so imprime o
relatorio do que seria gerado. --write grava de verdade e roda
registry.validate_all() no final.

Peculiaridades preservadas como dado (nao codigo), por decisao D-0009:
- citocalasina D tem Kd experimental que varia ~3 ordens de grandeza por
  ESTADO da actina (G-monomero 4.6 uM vs extremidade barbada do filamento
  ~2 nM) -- por isso experimental_kd[] carrega 'state' obrigatorio.
- jasplaquinolida so faz sentido em F-actina (assembly_state=filament);
  latrunculina/citocalasina primariamente em G-actina (assembly_state=monomer).
- CY9/LAR/LAB sao macrociclos RIGIDOS com 11-13 torsoes ativas -- ABAIXO do
  gate padrao da plataforma (GATE_MIN_TORSIONS=15 em runner.py). Nao
  contornado aqui: fica registrado em notes para quem for rodar docking real.
- TcAct2 e TcAct3 falham no veredito geometrico do Boltz-2 p/ citocalasina D
  ("outro_sitio") -- affinity_estimator_reliable=False universal (a cabeca de
  afinidade do Boltz-2 inverte a ordem experimental, ate fator 6x medido).
- O sitio da latrunculina NAO tem tabela de conservacao por residuo na fonte
  (so citocalasina_D/_Fform, jasplaquinolida, faloidina tem) -- registrado
  como lacuna real, nao preenchido por invencao.
- Faloidina e referencia de sitio (controle positivo), NUNCA docada -- sem
  entrada em ligands/.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import registry  # noqa: E402

CAMPAIGN_ID = "SUB04"

SOURCE_ROOT = Path(os.environ.get(
    "SUB04_SOURCE_ROOT",
    r"C:\Users\luizd\OneDrive\Data_Science\Projetos\Rho_OCRL\in_silico\drogas_citoesqueleto"
))

TARGET_IDS = ["TcAct1", "TcAct2", "TcAct3", "TcAct4", "Hs_ACTB"]
LIGAND_CCDS = ["CY9", "9UE", "LAR", "LAB"]

# ccd -> id do registry (minusculo/hifen, como registry._safe_id produziria)
LIGAND_ID = {"CY9": "citocalasina-d", "9UE": "jasplaquinolida",
             "LAR": "latrunculina-a", "LAB": "latrunculina-b"}

# ccd -> chave em receptors_qc.json[alvo]['sites'] / binding_sites.json (site
# 'primario' de cada ligante -- CY9 tem 2, ver SITE_KEYS_BY_CCD abaixo).
SITE_KEYS_BY_CCD = {
    "CY9": ["citocalasina_D", "citocalasina_D_Fform"],
    "9UE": ["jasplaquinolida"],
    "LAR": [],   # sem tabela de conservacao propria na fonte (ver docstring)
    "LAB": [],
}

# box_site usado em docking_prospective.json (2 caixas: citocalasina, latrunculina)
BOX_SITE_BY_CCD = {"CY9": "citocalasina", "LAR": "latrunculina", "LAB": "latrunculina"}

ORGANISM = {"TcAct1": "Trypanosoma cruzi", "TcAct2": "Trypanosoma cruzi",
            "TcAct3": "Trypanosoma cruzi", "TcAct4": "Trypanosoma cruzi",
            "Hs_ACTB": "Homo sapiens"}
CATEGORY = {"TcAct1": "parasite", "TcAct2": "parasite", "TcAct3": "parasite",
            "TcAct4": "parasite", "Hs_ACTB": "human"}


# ======================================================================
# carga dos JSONs de origem (leitura pura, nunca escreve em Rho_OCRL)
# ======================================================================

def _load(relpath):
    p = SOURCE_ROOT / relpath
    if not p.exists():
        raise FileNotFoundError(f"esperado em Rho_OCRL, nao encontrado: {p}")
    return json.loads(p.read_text(encoding="utf-8"))


def load_source():
    return {
        "receptors_qc": _load("analysis/receptors_qc.json"),
        "binding_sites": _load("sites/binding_sites.json"),
        "redock_gate": _load("analysis/redock_gate.json"),
        "docking_prospective": _load("analysis/docking_prospective.json"),
        "boltz_eval": _load("analysis/boltz_eval.json"),
        "factin_complexes": _load("analysis/factin_complexes.json"),
        "filament_assembly": _load("analysis/filament_assembly.json"),
        "cytod_controls": _load("analysis/cytod_controls.json"),
        "calibracao_afinidade": _load("analysis/calibracao_afinidade.json"),
        "blind_test_8trn": _load("analysis/blind_test_8trn.json"),
        "gpu_check": _load("analysis/gpu_check.json"),
        "ligands": _load("ligands/ligands.json"),
    }


# ======================================================================
# helpers
# ======================================================================

def _divergent_residues(site_block, target_id):
    """Residuos onde este alvo diverge do consenso ('ref') no sitio -- diff
    mecanico sobre dado ja carregado, nunca inventado."""
    out = []
    for chain, residues in (site_block.get("chains") or {}).items():
        for r in residues:
            val = r.get(target_id)
            if val is not None and val != r.get("ref"):
                out.append({"chain": chain, "resnum": r["resnum"],
                            "ref": r.get("ref"), "ref_resname": r.get("ref_resname"),
                            "target_residue": val})
    return out


def _boltz_job(src, label):
    for item in src["boltz_eval"]:
        if item.get("label") == label:
            return item
    return None


# ======================================================================
# construtores: alvo, ligante, par, campanha
# ======================================================================

def build_target(tid, src):
    qc = src["receptors_qc"][tid]
    sites = {}

    # box_size_A: a fonte (Rho_OCRL) NAO computou/gravou uma caixa de docking
    # explicita em nenhum dos 15 JSONs -- ela usou outro metodo (Boltz-2 co-
    # dobramento + Vina com box centrada no co-cristal via script proprio, sem
    # persistir as dimensoes). Nao inventamos um numero "medido": usamos o
    # PLACEHOLDER ja convencionado na plataforma (25 A, mesmo valor do TcTR em
    # SUB01) com box_policy sinalizando que precisa ser recalculada a partir do
    # ligante real antes de qualquer execucao nativa (Fase 7.2) -- nunca usar
    # para docking de verdade sem revisar.
    PLACEHOLDER_BOX = [25.0, 25.0, 25.0]

    for site_key in ("citocalasina_D", "citocalasina_D_Fform", "jasplaquinolida", "faloidina"):
        qc_site = qc["sites"].get(site_key, {})
        bs = src["binding_sites"].get(site_key, {})
        divergences = _divergent_residues(bs, tid) if bs else []
        sites[site_key] = {
            "definition": "derived_from_reference_complex",
            "reference_pdb": bs.get("pdb"),
            "reference_ligand_resname": bs.get("ligand"),
            "box_size_A": PLACEHOLDER_BOX,
            "box_policy": "xtal_bbox_plus_padding",
            "cutoff_A": bs.get("cutoff_A"),
            "n_contacts": bs.get("n_contacts"),
            "n_positions_mapped": qc_site.get("n_mapeados"),
            "n_positions_total": qc_site.get("n_total"),
            "plddt_medio": qc_site.get("plddt_medio"),
            "plddt_min": qc_site.get("plddt_min"),
            "divergent_residues_vs_consensus": divergences,
            "notes": (f"Sitio derivado de contatos <= {bs.get('cutoff_A')} A do "
                     f"ligante no co-cristal {bs.get('pdb')}. "
                     f"{len(divergences)} residuo(s) divergem do consenso multi-especie "
                     f"neste alvo. box_size_A e PLACEHOLDER (convencao da plataforma, "
                     f"nao medido na fonte) -- recalcular a partir do ligante real "
                     f"antes de rodar docking de verdade (Fase 7.2)."
                     if bs else
                     "Sem tabela de conservacao por residuo na fonte para este sitio. "
                     "box_size_A e PLACEHOLDER, nao medido."),
        }

    # sitio da latrunculina: SEM tabela de conservacao na fonte (gap real, nao
    # preenchido por invencao). Referencia = 2Q0U (latrunculina B, melhor
    # resolucao: 1.45 A vs 2.0 A do 1ESV/latrunculina A) -- ambos citados na nota.
    sites["latrunculina"] = {
        "definition": "derived_from_reference_complex",
        "reference_pdb": "2Q0U",
        "reference_ligand_resname": "LAB",
        "box_size_A": PLACEHOLDER_BOX,
        "box_policy": "xtal_bbox_plus_padding",
        "notes": "Fenda do nucleotideo. A fonte NAO computou tabela de conservacao "
                 "por residuo nem pLDDT por sitio para este bolso (diferente dos "
                 "outros 3 sitios, todos com binding_sites.json/receptors_qc.json "
                 "dedicados). Referencia = 2Q0U (latrunculina B, melhor resolucao); "
                 "1ESV (latrunculina A, 2.0 A) e co-cristal alternativo -- ver "
                 "pairs/*__latrunculina-a.json e *__latrunculina-b.json. box_size_A "
                 "e PLACEHOLDER, nao medido na fonte.",
    }

    divergencias_modelo = qc.get("divergencias_dm28c_vs_modelo", [])
    notes = (
        f"Modelo AlphaFold, identidade Dm28c x CL Brener = "
        f"{qc.get('identity_dm28c_vs_clb')}%. "
        + (f"{len(divergencias_modelo)} substituicao(oes) do modelo (baseado em CL "
           f"Brener) vs a cepa real Dm28c: {', '.join(divergencias_modelo)} -- "
           f"caveat de fidelidade da cepa, nao de conservacao de sitio."
           if divergencias_modelo else
           "Sem divergencias conhecidas entre a cepa Dm28c e o genoma-base CL Brener.")
    )

    return {
        "schema_version": "2.0",
        "id": tid,
        "aliases": [qc["uniprot"]],
        "name": ("Actina beta humana (referencia)" if tid == "Hs_ACTB"
                 else f"Actina {tid[-1]} de Trypanosoma cruzi (cepa Dm28c)"),
        "organism": ORGANISM[tid],
        "category": CATEGORY[tid],
        "ortholog_of": None,
        "structure": {
            "source": "alphafold",
            "pdb_id": None,
            "assembly": None,
            "n_protomers_required": 1,   # base monomerica; filamento e peculiaridade do PAR
            "local_file": str((SOURCE_ROOT / "structures" / qc["file"])),
            "confidence": {
                "global": qc["plddt_global"],
                "site": {k: v.get("plddt_medio") for k, v in qc["sites"].items()},
            },
            "resolution_A": None,
            "cofactor_copies": {},
        },
        "chemistry": {
            "keep_cofactors": [],
            "essential_cofactors": [],
            "cofactor_source": "none",
            "nonstandard_residues": [],
            "protonation_pH": 7.4,
        },
        "sites": sites,
        "default_site": "citocalasina_D",
        "notes": notes,
        "literature": "",
        "provenance": {
            "mode": "ingested",
            "reproducible": False,
            "generated_from": "Rho_OCRL/in_silico/drogas_citoesqueleto/analysis/receptors_qc.json "
                              "+ sites/binding_sites.json",
            "generator": "tools/ingest_campaign.py",
            "original_script": "prep_docking_receptors.py (Rho_OCRL, nao versionado aqui)",
        },
        "n_residues": qc.get("n_res"),
        "uniprot": qc.get("uniprot"),
    }


def build_ligand(ccd, src):
    lig = src["ligands"][ccd]
    mc = lig.get("macrocycle") or {}
    torsions = lig.get("pdbqt_active_torsions")
    lig_id = LIGAND_ID[ccd]

    gate_note = ""
    if torsions is not None and torsions < 15:
        gate_note = (f" ATENCAO: {torsions} torsoes ativas -- ABAIXO do gate padrao "
                     f"da plataforma (GATE_MIN_TORSIONS=15 em runner.py). Macrociclo "
                     f"rigido; o gate generico nao se aplica sem ajuste explicito.")

    return {
        "schema_version": "2.0",
        "id": lig_id,
        "name": lig.get("label"),
        "display_name": lig.get("name_pdb"),
        "smiles": lig.get("smiles"),
        "formula": (lig.get("formula") or "").replace(" ", ""),
        "mw_g_mol": lig.get("mw"),
        "heavy_atoms": None,
        "charge": 0,
        "class": "macrociclo natural (droga do citoesqueleto)",
        "chemotype_tags": ["cytoskeleton_drug", "macrocycle"],
        "topology": {
            "is_macrocycle": bool(mc.get("is_macrocycle")),
            "largest_ring_size": mc.get("largest_ring"),
            "ring_sizes": mc.get("ring_sizes"),
            "n_rings": mc.get("n_rings"),
            "expected_active_torsions": torsions,
        },
        "prep": {
            "optimization": f"MMFF94 ({lig.get('n_confs')} confs; melhor "
                            f"{lig.get('mmff_best_kcal')} kcal/mol)",
            "protonation": None,
            "pdbqt": {"tool": "meeko", "rigid_macrocycles": False,
                      "min_ring_size": 7, "max_ring_size": 33,
                      "notes": lig.get("pdbqt_note")},
            "gates": [],  # gate generico nao se aplica -- ver notes
        },
        "notes": f"Papel no estudo: {lig.get('role', '')}.{gate_note}",
        "ccd": ccd,
        "reference_pdb": lig.get("ref_pdb"),
        "reference_resolution_A": lig.get("ref_resolution_A"),
        "provenance": {
            "mode": "ingested",
            "reproducible": False,
            "generated_from": "Rho_OCRL/in_silico/drogas_citoesqueleto/ligands/ligands.json",
            "generator": "tools/ingest_campaign.py",
        },
    }


def build_pair(tid, ccd, src):
    lig_key = LIGAND_ID[ccd]
    pair_id = f"{CAMPAIGN_ID}__{tid}__{lig_key}"
    lig = src["ligands"][ccd]

    # ---- assembly_state: jasplaquinolida exige filamento; os demais, monomero
    if ccd == "9UE":
        fc = next((x for x in src["factin_complexes"] if x["receptor"] == tid), None)
        assembly = {
            "required": "filament",
            "min_protomers": 3,   # n_copias_9UE: protomeros de fato contatados pelo ligante
            "site_topology": "interprotomer (juncao de 3 protomeros)",
            "rationale": "Jasplaquinolida liga na juncao de multiplos protomeros da "
                         "F-actina; nao existe pose monomerica com sentido biologico "
                         "(PLANO.md da fonte).",
            "built_by": "transferencia estrutural / threading (build_factin_complexes.py, "
                        "Rho_OCRL) -- NAO e predicao de novo.",
        }
        if fc:
            assembly["n_protomeros_no_modelo"] = fc.get("n_protomeros")
            assembly["rmsd_colocacao_medio_A"] = fc.get("rmsd_colocacao_medio_A")
            assembly["choques_interprotomero_lt_2A"] = fc.get("choques_interprotomero_lt_2A")
    else:
        assembly = {
            "required": "monomer",
            "min_protomers": 1,
            "site_topology": "intraprotomer",
            "rationale": "Docking prospectivo e redocking rodados no monomero "
                         "(G-actina). Ver validation.experimental_kd para o "
                         "efeito de ate 3 ordens de grandeza do estado "
                         "(monomero x extremidade barbada do filamento).",
        }

    # ---- experimental_kd: multi-estado, de calibracao_afinidade.escala_experimental
    # 'value'+'unit' e 'state' sao exigidos por registry.validate_all(): comparar
    # afinidade sem o estado do alvo e erro de categoria (a razao de ser deste campo).
    exp_kd = []
    for row in src["calibracao_afinidade"].get("escala_experimental", []):
        if row.get("ccd") == ccd:
            exp_kd.append({
                "value": row.get("kd_uM"),
                "unit": "uM",
                "type": row.get("tipo"),
                "state": row.get("estado_actina"),
                "comparable_to_monomer_prediction": row.get("comparavel_ao_monomero"),
                "pmid": row.get("pmid"),
                "doi": row.get("doi"),
                "citation": row.get("citacao"),
            })

    # ---- validation: redock gate (protocolo, nao por-alvo) + boltz por-alvo
    redock_entries = [(k, v) for k, v in src["redock_gate"].items() if v.get("ccd") == ccd]
    validation = {"redock_gate": [
        {"reference_pdb": v.get("pdb"), "rmsd_best_A": v.get("rmsd_best"),
         "rmsd_top1_A": v.get("rmsd_top1"), "passed": v.get("passou"),
         "rmsd_function": "menor RMSD entre 20 poses vs cristal (equivalente a "
                          "CalcRMS -- nunca GetBestRMS, que alinha antes de medir "
                          "e mascara pose errada)"}
        for _, v in redock_entries
    ]}

    boltz = _boltz_job(src, f"{tid}_x_{ccd}")
    if boltz:
        validation["boltz2_geometry"] = {
            "verdict": (boltz.get("geometria") or {}).get("veredito"),
            "ptm": (boltz.get("confianca") or {}).get("ptm"),
            "iptm": (boltz.get("confianca") or {}).get("iptm"),
            "complex_plddt": (boltz.get("confianca") or {}).get("complex_plddt"),
        }

    # ---- docking prospectivo Vina (site nativo + cruzado)
    docking = {}
    for key, entry in src["docking_prospective"].items():
        t, c, box = key.split("__")
        if t == tid and c == ccd:
            docking[box] = {"is_native_site": entry.get("is_native_site"),
                            "dG_top1_kcal_mol": entry.get("dG_top1")}

    # ---- controles G/F para citocalasina D (cytod_controls.json)
    state_controls = [
        {"state": c["estado"], "reference_pdb": c["cristal"],
         "kd_experimental_uM": c["kd_experimental_uM"],
         "rmsd_superposicao_A": c["rmsd_superposicao_A"],
         "recuperacao_sitio_pct": c["recuperacao_sitio_pct"]}
        for c in src["cytod_controls"] if c["receptor"] == tid and ccd == "CY9"
    ]

    # ---- treinamento/vazamento: registry.validate_all() exige 'held_out' explicito
    # se o bloco existir. O teste cego FORMAL na fonte e para TgAct1/8TRN (T. gondii),
    # um alvo diferente destes 5 -- nao temos held_out firme PARA ESTE par
    # especificamente, entao NAO populamos validation.training_leakage aqui (evita
    # afirmar o que nao sabemos). O contexto fica em status.note, sem gate de schema.
    leakage_note = None
    if ccd == "9UE":
        bt = src["blind_test_8trn"]
        leakage_note = (f"Teste cego formal de vazamento de treino do Boltz-2 existe "
                        f"para {bt.get('alvo')}/8TRN (T. gondii, NAO este alvo) -- "
                        f"alcance_da_aplicacao varia 40-74% conforme o TcAct "
                        f"(corte_treino_boltz2={bt.get('corte_treino_boltz2')}). "
                        f"Ver blind_test_8trn.json na fonte para o detalhe por alvo; "
                        f"nao ha held_out determinado diretamente para este par.")

    engine_primary = "vina" if ccd != "9UE" else "structural_transfer_threading"

    return {
        "schema_version": "2.0",
        "campaign": CAMPAIGN_ID,
        "target": tid,
        "ligand": lig_key,
        "site": "citocalasina_D" if ccd == "CY9" else
                ("jasplaquinolida" if ccd == "9UE" else "latrunculina"),
        "role": "primary_target" if tid.startswith("TcAct") else "human_reference",

        "assembly_state": assembly,

        "cofactor_requirement": {
            "required": [],
            "rationale": "Sem cofator obrigatorio identificado na fonte para este "
                         "sistema (actina + droga do citoesqueleto). Estado "
                         "nucleotidico (ATP/ADP) da actina e variavel relevante "
                         "para afinidade (ver experimental_kd[].state), mas nao "
                         "modelado como cofator explicito nesta ingestao.",
        },

        "experimental_kd": exp_kd,
        "_experimental_kd_note": "Kd de literatura (PubMed), state obrigatorio por "
                                 "peculiaridade cientifica real: a afinidade varia ate "
                                 "3 ordens de grandeza com o estado da actina.",

        "state_controls": state_controls,   # so preenchido p/ CY9 (cytod_controls.json)

        "engines": {
            "primary_pose": engine_primary,
            "secondary_pose": "boltz2",
            "affinity_estimator_reliable": False,
            "affinity_estimator_note": "calibracao_afinidade.json: a cabeca de afinidade "
                                       "do Boltz-2 inverte a ordem experimental em pelo "
                                       "menos 1 caso medido (fator de erro ate 6x p/ "
                                       "latrunculina A). Score/afinidade do Boltz-2 NAO "
                                       "confiavel para ranking absoluto neste sistema.",
        },

        "validation": validation,
        "docking_prospective": docking,

        "protocols": {"docking": "vina_exh32"} if ccd != "9UE" else {},

        "status": {"docking": "done" if ccd != "9UE" else "n/a (filamento, sem docking monomerico)",
                   "structure_prediction": "done",
                   "note": ("Ingerido de Rho_OCRL/in_silico/drogas_citoesqueleto. "
                           "Execucao nativa pela plataforma = Fase 7.2, pendente."
                           + (f" {leakage_note}" if leakage_note else ""))},

        "provenance": {
            "mode": "ingested",
            "reproducible": False,
            "generated_from": "Rho_OCRL/in_silico/drogas_citoesqueleto/analysis/*.json",
            "generator": "tools/ingest_campaign.py",
            "original_scripts": ["prep_docking_receptors.py", "redock_gate.py",
                                 "run_boltz.py", "build_factin_complexes.py"],
        },
        "_pair_id": pair_id,
    }


def build_campaign(targets, ligands, pairs, src):
    return {
        "schema_version": "2.0",
        "id": CAMPAIGN_ID,
        "sigla": "TcAct-citoesqueleto",
        "programa": "Rho_OCRL-actina",
        "nome_especifico": "TcActins (Dm28c) x Citocalasina D / Jasplaquinolida / "
                           "Latrunculina A / Latrunculina B -- selecao de paraloga "
                           "e estado oligomerico como determinantes de ligacao",
        "titulo_paper_atual": "TcActins (Dm28c) x Cytochalasin D / Jasplakinolide / "
                              "Latrunculin A / Latrunculin B -- paralog selection "
                              "and oligomeric state as binding determinants",
        "status": "ingerido (dados de trabalho; manuscrito em preparacao no Rho_OCRL)",
        "fase_atual": "Ingestao read-only (Fase 7.1). Execucao nativa pela plataforma "
                      "= Fase 7.2, ainda nao feita.",
        "periodico_alvo": None,
        "criado": "2026-08-09",

        "targets": TARGET_IDS,
        "primary_targets": ["TcAct1", "TcAct2", "TcAct3", "TcAct4"],
        "counter_targets": [],
        "human_reference_targets": ["Hs_ACTB"],

        "ligands": [LIGAND_ID[c] for c in LIGAND_CCDS],
        "default_ligand": LIGAND_ID["CY9"],

        "protocols": {"docking": "vina_exh32"},

        "paths": {"data_root": None,
                  "_note": "null = usa o default (docking_project/SUB04 dentro da "
                           "plataforma). Os dados brutos de origem continuam em "
                           f"{SOURCE_ROOT} -- nunca movidos."},

        "ui": {
            "header": "Actinas de Trypanosoma cruzi x drogas do citoesqueleto",
            "subtitle": "4 paralogas (Dm28c) vs actina beta humana -- "
                        "conservacao de sitio como determinante de seletividade",
            "recommended_targets": ["TcAct1", "TcAct3", "Hs_ACTB"],
            "organism_badge_field": "category",
        },

        "provenance": {
            "mode": "ingested",
            "reproducible": False,
            "source_root": str(SOURCE_ROOT),
            "generator": "tools/ingest_campaign.py",
            "note": "Ingestao read-only da Fase 7.1. A plataforma le os dados de onde "
                    "estao; o projeto Rho_OCRL continua dono deles. Ver "
                    ".project-memory/PROJECT_STATE.md e PLANO.md em Rho_OCRL para a "
                    "narrativa cientifica completa (nao reescrita aqui).",
        },
    }


# ======================================================================
# relatorio (sempre impresso, dry-run ou nao) + escrita
# ======================================================================

def report(targets, ligands, pairs, campaign):
    print("=" * 72)
    print(f"INGESTAO {CAMPAIGN_ID} -- relatorio (fonte: {SOURCE_ROOT})")
    print("=" * 72)
    print(f"\nAlvos ({len(targets)}):")
    for t in targets:
        conf = t["structure"]["confidence"]["global"]
        print(f"  {t['id']:10s} pLDDT={conf:5.1f}  {t['organism']:20s} {t['category']}")
    print(f"\nLigantes ({len(ligands)}):")
    for l in ligands:
        torsions = l["topology"]["expected_active_torsions"]
        flag = " <<< abaixo do gate (15)" if torsions and torsions < 15 else ""
        print(f"  {l['id']:18s} {l['ccd']:5s} torsoes={torsions}{flag}")
    print(f"\nPares ({len(pairs)}):")
    for p in pairs:
        assembly = p["assembly_state"]["required"]
        n_kd = len(p["experimental_kd"])
        boltz_v = (p["validation"].get("boltz2_geometry") or {}).get("verdict", "-")
        print(f"  {p['target']:10s} x {p['ligand']:18s} assembly={assembly:8s} "
              f"kd_states={n_kd}  boltz2={boltz_v}")
    print()


def write_all(targets, ligands, pairs, campaign, out_dir):
    out_dir = Path(out_dir)
    written = []

    def _write(subdir, id_, data):
        p = out_dir / subdir / f"{id_}.json"
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        written.append(p)

    _write("campaigns", CAMPAIGN_ID, campaign)
    for t in targets:
        _write("targets", t["id"], t)
    for l in ligands:
        _write("ligands", l["id"], l)
    for p in pairs:
        pid = p.pop("_pair_id")
        _write("pairs", pid, p)

    print(f"\n{len(written)} arquivo(s) gravado(s) em {out_dir}\n")
    return written


# ======================================================================
# CLI
# ======================================================================

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true",
                    help="grava de verdade em registry/. Sem isso: dry-run (so relatorio).")
    ap.add_argument("--out", default=str(ROOT / "registry"),
                    help="raiz do registry a gravar (default: registry/ do repo).")
    args = ap.parse_args(argv)

    src = load_source()

    targets = [build_target(t, src) for t in TARGET_IDS]
    ligands = [build_ligand(c, src) for c in LIGAND_CCDS]
    pairs = [build_pair(t, c, src) for t in TARGET_IDS for c in LIGAND_CCDS]
    campaign = build_campaign(targets, ligands, pairs, src)

    report(targets, ligands, pairs, campaign)

    if not args.write:
        print("DRY-RUN (default): nada foi gravado. Use --write para persistir em "
              f"{args.out}/{{campaigns,targets,ligands,pairs}}/.")
        return 0

    write_all(targets, ligands, pairs, campaign, args.out)

    registry.clear_cache()
    errs = registry.validate_all()
    if errs:
        print(f"[AVISO] {len(errs)} erro(s) de validacao apos a escrita:")
        for e in errs:
            print(f"  - {e}")
        return 1
    print("[ok] registry.validate_all() limpo apos a ingestao.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
