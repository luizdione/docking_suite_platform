# -*- coding: utf-8 -*-
"""tools/validate_derived_sites.py -- gate GEOMETRICO (sem docking, sem GPU) dos
sitios 'derived_from_reference_complex' de uma campanha.

Por que geometrico e nao redock rigido: os alvos do SUB04 sao modelos AlphaFold
*apo*. Redock rigido do ligante nesses modelos NAO reproduz a pose cristalografica
(bolso parcialmente fechado; o Vina otimiza energia e foge do bolso) -- validado
empiricamente (ate o caso de superposicao 0.70 A da ~8 A de RMSD). O gate honesto:
a pose de referencia, transformada p/ o frame do alvo, cai num BOLSO real do modelo?
(fracao de atomos do ligante em contato < 4.5 A com a proteina, e imersao media).
Discrimina bolso presente x ausente -- ex.: jasplaquinolida (liga FILAMENTO) da
imersao fraca num modelo MONOMERICO, biologicamente correto.

DRY-RUN por padrao (so relata). --write persiste box_center_A (target_frame),
box_size_A e as metricas nos sitios com veredito 'center_ok', com proveniencia
honesta. Sitios sem veredito 'center_ok' seguem SEM box_center_A -> fail-loud no
runner (nunca centro-de-massa em silencio).

Uso:
    python tools/validate_derived_sites.py --campaign SUB04           # relatorio
    python tools/validate_derived_sites.py --campaign SUB04 --write   # persiste os OK
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import registry            # noqa: E402
import site_derivation     # noqa: E402

DL_DIR = ROOT / ".pdb_cache"


def _iter_derived(campaign_id):
    for target in registry.list_targets(campaign_id):
        tid = target.get("id")
        cif = registry.resolve_local_file((target.get("structure") or {}).get("local_file"))
        for sname, site in (target.get("sites") or {}).items():
            if site.get("definition") == "derived_from_reference_complex":
                yield tid, sname, site, cif


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--campaign", default="SUB04")
    ap.add_argument("--write", action="store_true",
                    help="persiste box_center_A dos sitios com veredito center_ok")
    args = ap.parse_args()
    DL_DIR.mkdir(exist_ok=True)

    rows, ok, local, weak, blocked = [], 0, 0, 0, 0
    print(f"{'ALVO':9}{'SITIO':22}{'ref':6}{'supRMSD':>8}{'imersao':>9}{'contato':>8}  veredito")
    for tid, sname, site, cif in _iter_derived(args.campaign):
        rec = {"target": tid, "site": sname, "ref": site.get("reference_pdb")}
        if not cif or not Path(cif).exists():
            rec["verdict"] = "blocked_no_cif"; blocked += 1; rows.append(rec)
            print(f"{tid:9}{sname:22}{str(rec['ref']):6}{'':>8}{'':>9}{'':>8}  BLOQUEADO: .cif ausente")
            continue
        try:
            g = site_derivation.geometric_validation(
                cif, site["reference_pdb"], site.get("reference_ligand_resname", ""), DL_DIR)
        except Exception as exc:  # noqa: BLE001
            rec["verdict"] = f"blocked:{type(exc).__name__}"; blocked += 1; rows.append(rec)
            print(f"{tid:9}{sname:22}{str(rec['ref']):6}  BLOQUEADO: {type(exc).__name__}: {str(exc)[:45]}")
            continue
        rec.update(g)
        v = g.get("verdict")
        if v == "center_ok":
            ok += 1
        elif v == "center_ok_local":
            local += 1
        elif v == "weak_immersion":
            weak += 1
        else:
            blocked += 1
        rows.append(rec)
        if "frac_contact" in g:
            print(f"{tid:9}{sname:22}{str(rec['ref']):6}{g['superposition_rmsd_A']:>8}"
                  f"{g['immersion_mean_A']:>8}A{g['frac_contact']:>8}  {v}")
        else:
            print(f"{tid:9}{sname:22}{str(rec['ref']):6}{g.get('superposition_rmsd_A',''):>8}"
                  f"{'':>9}{'':>8}  {v}")

    print(f"\ncenter_ok: {ok}   center_ok_local: {local}   imersao_fraca: {weak}   "
          f"bloqueado: {blocked}   total: {len(rows)}")

    if args.write:
        written = 0
        for r in rows:
            if r.get("verdict") not in ("center_ok", "center_ok_local"):
                continue
            is_local = r["verdict"] == "center_ok_local"
            path = registry.REGISTRY_DIR / "targets" / f"{r['target']}.json"
            obj = json.loads(path.read_text(encoding="utf-8"))
            s = obj["sites"][r["site"]]
            s["box_center_A"] = r["box_center_A_target_frame"]
            s["box_size_A"] = r["box_size_A"]
            s["box_center_frame"] = "target"
            s["superposition_rmsd_A"] = r["superposition_rmsd_A"]
            s["n_ca_matched"] = r["n_ca_matched"]
            s["superposition_mode"] = r.get("superposition_mode", "global")
            gv = {
                "frac_contact_4.5A": r["frac_contact"],
                "immersion_mean_A": r["immersion_mean_A"],
                "immersion_min_A": r["immersion_min_A"],
                "method": "site_derivation.geometric_validation",
            }
            if is_local:
                gv["superposition_mode"] = "local_pocket"
                gv["local_rmsd_A"] = r.get("local_rmsd_A")
                gv["n_ca_local"] = r.get("n_ca_local")
                gv["note"] = (
                    "RESGATADO por superposicao LOCAL do bolso: a dobra GLOBAL diverge "
                    f"(RMSD global {r['superposition_rmsd_A']} A > 3 A, paralogo com "
                    "insercoes fora do sitio), mas o bolso e localmente conservado "
                    f"(RMSD local {r.get('local_rmsd_A')} A sobre {r.get('n_ca_local')} "
                    "CAs) e a pose de referencia continua IMERSA nos atomos reais do "
                    "modelo (target_frame, receptor AlphaFold intocado). Redock rigido "
                    "NAO reproduz pose em modelo apo -- para pose exata usar docking "
                    "flexivel/template holo.")
            else:
                gv["superposition_mode"] = "global"
                gv["note"] = (
                    "Centro por superposicao estrutural GLOBAL (target_frame, receptor "
                    "AlphaFold intocado); validado GEOMETRICAMENTE (pose de referencia "
                    "enterrada no bolso). Redock rigido NAO reproduz a pose em modelo "
                    "AlphaFold apo (limitacao documentada) -- para pose exata usar "
                    "docking flexivel/template holo.")
            s["geometric_validation"] = gv
            path.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
            written += 1
        print(f"\nPERSISTIDO: {written} sitios (box_center_A target_frame; {local} por "
              f"superposicao local). Rode 'python registry.py' e os gates. "
              f"Os demais seguem fail-loud.")
    else:
        print("\nDRY-RUN: nada escrito. Use --write para persistir os center_ok.")


if __name__ == "__main__":
    main()
