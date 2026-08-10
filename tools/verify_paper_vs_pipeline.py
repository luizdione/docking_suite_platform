# -*- coding: utf-8 -*-
"""
tools/verify_paper_vs_pipeline.py -- Fase 8.1: diverge PAPER (export_web.py, dados
do manuscrito submetido) vs os resultados REAIS do pipeline (docking_results.json,
admet_profile.json, binding-mode da MD).

READ-ONLY: nunca escreve em export_web.py, no registry, nem republica o viewer.
Objetivo desta fase e so ENUMERAR as divergencias -- resolve-las (ou decidir que
sao esperadas/documentadas) e trabalho humano, antes de qualquer troca de fonte
(Fase 8.2). Gate: "resolver cada uma antes de trocar a fonte" (plano, Fase 8).

Fontes comparadas:
- Docking : results_io.read_results('SUB01', ligante) -- pipeline real.
- ADMET   : <data_root>/results/<ligante>/admet_profile.json -- pipeline real.
- MD      : binding_mode_<alvo>_<modo>[_r2|_r3].json no kit de figuras do paper
            (NQ_KIT_REDOCK ou o default hardcoded em dashboard_v2.py) -- estes SAO
            dados reais de MD (nao reinventados aqui), so NAO passam por
            results_io/registry ainda -- por isso o script le do mesmo lugar que
            a rota /binding_mode_status do dashboard ja le.

Uso: python tools/verify_paper_vs_pipeline.py [--campaign SUB01] [--ligand 2-SCN-NQ]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import paths          # noqa: E402
import results_io     # noqa: E402
import export_web     # noqa: E402

TOL_DG = 0.05          # kcal/mol -- casas decimais do PAPER, tolerancia de arredondamento
TOL_KI = 0.05          # uM relativo
TOL_MD_A = 0.15        # Angstrom -- ruido normal de MD entre reamostragens/versao do OpenMM
TOL_ADMET = 0.05        # relativo, campos numericos


# ======================================================================
# comparadores genericos
# ======================================================================

class Divergence:
    def __init__(self, section, key, field, paper_val, pipeline_val, note=""):
        self.section, self.key, self.field = section, key, field
        self.paper_val, self.pipeline_val, self.note = paper_val, pipeline_val, note

    def __str__(self):
        return (f"[{self.section}] {self.key}.{self.field}: "
                f"PAPER={self.paper_val!r}  pipeline={self.pipeline_val!r}"
                + (f"  ({self.note})" if self.note else ""))


def _num_diverges(a, b, tol, rel=False):
    if a is None or b is None:
        return a != b
    try:
        a, b = float(a), float(b)
    except (TypeError, ValueError):
        return a != b
    if rel:
        denom = max(abs(a), abs(b), 1e-9)
        return abs(a - b) / denom > tol
    return abs(a - b) > tol


# ======================================================================
# 1) DOCKING
# ======================================================================

def verify_docking(campaign, ligand):
    divs = []
    live = results_io.read_results(campaign, ligand)
    if live is None:
        divs.append(Divergence("docking", "*", "*", "existe", None,
                               "results_io.read_results() nao achou nenhum resultado -- "
                               "nada para comparar"))
        return divs, {}

    paper_by_code = {row["code"]: row for row in export_web.PAPER["docking"]}
    live_codes = {k for k in live if not k.startswith("_")}
    paper_codes = set(paper_by_code)

    only_in_paper = paper_codes - live_codes
    only_in_pipeline = live_codes - paper_codes
    if only_in_paper:
        divs.append(Divergence("docking", "*", "alvos", sorted(only_in_paper), None,
                               "codigo(s) do PAPER sem correspondente EXATO no pipeline -- "
                               "conferir se e so diferenca de nome (ex.: LdDHODH vs LmDHODH)"))
    if only_in_pipeline:
        divs.append(Divergence("docking", "*", "alvos", None, sorted(only_in_pipeline),
                               "alvo(s) no pipeline sem entrada no PAPER"))

    for code in sorted(paper_codes & live_codes):
        p, l = paper_by_code[code], live[code]
        if _num_diverges(p.get("dg"), l.get("best_dG"), TOL_DG):
            divs.append(Divergence("docking", code, "dG_kcal_mol", p.get("dg"), l.get("best_dG")))
        if _num_diverges(p.get("ki_uM"), l.get("ki_uM"), TOL_KI, rel=True):
            divs.append(Divergence("docking", code, "ki_uM", p.get("ki_uM"), l.get("ki_uM")))
        if p.get("pdb") and l.get("pdb") and p["pdb"] != l["pdb"]:
            divs.append(Divergence("docking", code, "pdb", p.get("pdb"), l.get("pdb")))
        p_cof = set(p.get("cofactors") or [])
        l_cof = set(l.get("cofactors_kept") or [])
        if p_cof != l_cof:
            divs.append(Divergence("docking", code, "cofactors", sorted(p_cof), sorted(l_cof)))

    return divs, {"live_source": live.get("_source_path"), "n_targets_compared": len(paper_codes & live_codes)}


# ======================================================================
# 2) ADMET
# ======================================================================

# PAPER usa nomes proprios; admet_profile.json (pipeline) usa outros. Mapa
# EXPLICITO em vez de comparar por posicao -- assim um campo faltante aparece
# como divergencia, nao e silenciosamente pulado.
ADMET_FIELD_MAP = {
    "MW": "MW", "cLogP": "LogP", "TPSA": "TPSA", "HBD": "HBD", "HBA": "HBA",
    "RB": "RotBonds", "QED": "QED", "bioavailability": "bioavailability_score",
    "hia": "hia_percent", "logBB": None, "ppb": None, "vd": None,
    "ames": "ames_alert", "pains": None, "t_half": "t_half_hours",
    "caco2": "Caco2_class", "bbb": "bbb_penetrant", "herg": "herg_risk",
}


def _admet_candidates(campaign, ligand):
    """Cascata de candidatos, igual em espirito a paths.candidate_results_paths():
    paths.data_root() depende de NQ_DATA_DIR estar setado NO PROCESSO ATUAL para
    resolver a raiz legada (E:/data/docking) -- sem isso cai no layout novo, que
    para o SUB01 ainda nao tem o arquivo. Tentamos os dois em vez de confiar em
    uma unica resolucao dependente de ambiente."""
    seen, out = set(), []
    for base in (paths.data_root(campaign), Path(os.environ.get("NQ_DATA_DIR", r"E:\data\docking"))):
        p = base / "results" / ligand / "admet_profile.json"
        key = str(p).lower()
        if key not in seen:
            seen.add(key)
            out.append(p)
    return out


def verify_admet(campaign, ligand):
    divs = []
    admet_path = next((p for p in _admet_candidates(campaign, ligand) if p.exists()), None)
    if admet_path is None:
        tried = ", ".join(str(p) for p in _admet_candidates(campaign, ligand))
        divs.append(Divergence("admet", "*", "*", "existe", None,
                               f"nao encontrado em nenhum candidato: {tried}"))
        return divs, {}
    live = json.loads(admet_path.read_text(encoding="utf-8"))
    paper = export_web.PAPER["admet"]

    for paper_key, live_key in ADMET_FIELD_MAP.items():
        if live_key is None:
            divs.append(Divergence("admet", "*", paper_key, paper.get(paper_key), "(sem mapa)",
                                   "campo do PAPER sem equivalente identificado no "
                                   "admet_profile.json -- conferir manualmente"))
            continue
        pv, lv = paper.get(paper_key), live.get(live_key)
        if isinstance(pv, (int, float)) and isinstance(lv, (int, float)):
            if _num_diverges(pv, lv, TOL_ADMET, rel=True):
                divs.append(Divergence("admet", "*", paper_key, pv, lv))
        elif str(pv).strip().lower() != str(lv).strip().lower():
            divs.append(Divergence("admet", "*", paper_key, pv, lv))

    return divs, {"live_source": str(admet_path)}


# ======================================================================
# 3) MD -- binding mode (lig-FAD / lig-ditiol), n=3 onde aplicavel
# ======================================================================

def _kit_dir():
    kit = os.environ.get("NQ_KIT_REDOCK") or (
        r"C:\Users\luizd\OneDrive\Data_Science\Projetos\2-SNC-NQ"
        r"\montagem_resultados\fig13_3d_kit\redock")
    return Path(kit)


# alvo -> (tag da replica 1, tags das replicas extras) -- os nomes de arquivo
# NAO sao uniformes na fonte (TcTR usa sufixo _50ns na r1; LiTR/hGR nao usam
# sufixo na r1). Mapear explicitamente em vez de adivinhar padrao.
MD_REPLICA_TAGS = {
    "TcTR": ["TcTR_subversive_50ns", "TcTR_subversive_r2", "TcTR_subversive_r3"],
    "LiTR": ["LiTR_subversive", "LiTR_subversive_r2", "LiTR_subversive_r3"],
    "hGR":  ["hGR_subversive", "hGR_subversive_r2", "hGR_subversive_r3"],
}
MD_INHIBITOR_TAGS = {"TcTR": "TcTR_inhibitor", "LiTR": "LiTR_inhibitor"}


def _load_binding_mode(kit, tag):
    p = kit / f"binding_mode_{tag}.json"
    if not p.exists():
        return None
    return json.loads(p.read_text(encoding="utf-8"))


def verify_md(campaign, ligand):
    divs = []
    kit = _kit_dir()
    if not kit.exists():
        divs.append(Divergence("md", "*", "*", "existe", None, f"kit nao encontrado em {kit}"))
        return divs, {}

    paper_subversive = {row["target"]: row for row in export_web.PAPER["md"]["subversive"]}
    n3_meta = {}
    for target, tags in MD_REPLICA_TAGS.items():
        replicas = [_load_binding_mode(kit, t) for t in tags]
        found = [r for r in replicas if r is not None]
        missing = [t for t, r in zip(tags, replicas) if r is None]
        if missing:
            divs.append(Divergence("md", target, "replicas", f"{len(tags)} esperadas",
                                   f"{len(found)} achadas", f"faltando: {missing}"))
        if not found:
            continue
        fad_means = [r["lig_fad_A"]["mean"] for r in found if "lig_fad_A" in r]
        dithiol_means = [r["lig_dithiol_A"]["mean"] for r in found if "lig_dithiol_A" in r]
        n3_meta[target] = {"n_replicas_achadas": len(found),
                           "n_frames_total": sum(r.get("n_frames", 0) for r in found)}
        if fad_means:
            mean_of_means = sum(fad_means) / len(fad_means)
            paper_row = paper_subversive.get(target)
            if paper_row and _num_diverges(paper_row.get("fad"), mean_of_means, TOL_MD_A):
                divs.append(Divergence("md", target, "fad_A (media das replicas)",
                                       paper_row.get("fad"), round(mean_of_means, 3)))
        if dithiol_means:
            n3_meta.setdefault(target, {})["dithiol_A_mean"] = round(
                sum(dithiol_means) / len(dithiol_means), 3)
        # dithiol_pct (% de frames com contato) NAO e recomputavel a partir do
        # resumo disponivel (so min/mean/max/first/last por replica, sem
        # timeseries) -- reportado como lacuna, nao inventado.
        if target in paper_subversive and paper_subversive[target].get("dithiol_pct") is not None:
            divs.append(Divergence("md", target, "dithiol_pct", paper_subversive[target]["dithiol_pct"],
                                   "NAO RECOMPUTAVEL",
                                   "resumo binding_mode_*.json so tem min/mean/max/first/last, "
                                   "sem serie temporal -- precisaria do CSV/DCD bruto p/ recalcular"))

    paper_inhibitor = {row["target"]: row for row in export_web.PAPER["md"]["inhibitor"]}
    for target, tag in MD_INHIBITOR_TAGS.items():
        r = _load_binding_mode(kit, tag)
        if r is None:
            divs.append(Divergence("md", target, "inhibitor", "existe", None, f"tag {tag} nao achada"))
            continue
        paper_row = paper_inhibitor.get(target, {})
        if "lig_fad_A" in r and _num_diverges(paper_row.get("fad"), r["lig_fad_A"]["mean"], TOL_MD_A):
            divs.append(Divergence("md", target, "fad_A (inhibitor)",
                                   paper_row.get("fad"), round(r["lig_fad_A"]["mean"], 3)))

    return divs, {"kit_dir": str(kit), "n3": n3_meta}


# ======================================================================
# relatorio
# ======================================================================

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--campaign", default="SUB01")
    ap.add_argument("--ligand", default="2-SCN-NQ")
    args = ap.parse_args(argv)

    print("=" * 72)
    print(f"verify_paper_vs_pipeline.py -- campanha={args.campaign} ligante={args.ligand}")
    print("READ-ONLY: nao altera export_web.py, registry, nem republica o viewer.")
    print("=" * 72)

    all_divs = []
    for name, fn in (("DOCKING", verify_docking), ("ADMET", verify_admet), ("MD", verify_md)):
        print(f"\n--- {name} ---")
        divs, meta = fn(args.campaign, args.ligand)
        for k, v in meta.items():
            print(f"  ({k}: {v})")
        if not divs:
            print("  [ok] sem divergencias.")
        for d in divs:
            print(f"  - {d}")
        all_divs.extend(divs)

    print("\n" + "=" * 72)
    print(f"TOTAL: {len(all_divs)} divergencia(s)/lacuna(s) encontrada(s).")
    print("Fase 8.2 (trocar a fonte do viewer p/ o pipeline) NAO deve prosseguir")
    print("ate cada uma destas ser resolvida ou explicitamente aceita como esperada.")
    print("=" * 72)
    return 1 if all_divs else 0


if __name__ == "__main__":
    raise SystemExit(main())
