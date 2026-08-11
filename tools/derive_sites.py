# -*- coding: utf-8 -*-
"""tools/derive_sites.py -- deriva (por superposicao) o centro dos sitios
'derived_from_reference_complex' de uma campanha e RELATA. Por padrao NAO escreve
nada (dry-run): imprime uma tabela para revisao humana. So com --write persiste
box_center_A/box_size_A/superposition_rmsd_A/n_ca_matched no registry, e mesmo
assim SO deveria ser usado depois que um redock gate confirmar a convencao.

Convencao: por padrao grava o centro no referencial do ALVO (target_frame,
receptor AlphaFold intocado). --crystal-frame grava o do cristal (exige receptor
movido, prepared_receptor_local_file -- nao gerado por este script).

Uso:
    python tools/derive_sites.py --campaign SUB04            # dry-run (relatorio)
    python tools/derive_sites.py --campaign SUB04 --write    # persiste (target_frame)
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import registry  # noqa: E402
import site_derivation  # noqa: E402

DL_DIR = ROOT / ".pdb_cache"


def _iter_derived_sites(campaign_id):
    # registry.list_targets() devolve os OBJETOS (dicts), nao os ids.
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
                    help="persiste no registry (default: so relata)")
    ap.add_argument("--crystal-frame", action="store_true",
                    help="grava centro no frame do cristal (exige receptor movido)")
    ap.add_argument("--rmsd-max", type=float, default=3.0)
    args = ap.parse_args()

    DL_DIR.mkdir(exist_ok=True)
    rows, ok, blocked = [], 0, 0
    for tid, sname, site, cif in _iter_derived_sites(args.campaign):
        rec = {"target": tid, "site": sname,
               "reference_pdb": site.get("reference_pdb")}
        if not cif or not Path(cif).exists():
            rec["status"] = f"BLOQUEADO: .cif ausente ({cif})"
            blocked += 1
            rows.append(rec)
            continue
        try:
            d = site_derivation.derive_site_center(
                cif, site["reference_pdb"], site.get("reference_ligand_resname", ""),
                DL_DIR, rmsd_max_A=args.rmsd_max)
            rec.update(status="OK", rmsd=d["superposition_rmsd_A"],
                       n_ca=d["n_ca_matched"],
                       center_target=d["box_center_A_target_frame"],
                       center_crystal=d["box_center_A_crystal_frame"],
                       box_size=d["box_size_A"])
            ok += 1
            if args.write:
                key = "box_center_A_crystal_frame" if args.crystal_frame else "box_center_A_target_frame"
                path = registry.REGISTRY_DIR / "targets" / f"{tid}.json"
                obj = json.loads(path.read_text(encoding="utf-8"))
                s = obj["sites"][sname]
                s["box_center_A"] = d[key]
                s["box_size_A"] = d["box_size_A"]
                s["superposition_rmsd_A"] = d["superposition_rmsd_A"]
                s["n_ca_matched"] = d["n_ca_matched"]
                s["box_center_frame"] = "crystal" if args.crystal_frame else "target"
                path.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + "\n",
                                encoding="utf-8")
        except Exception as exc:  # noqa: BLE001
            rec["status"] = f"BLOQUEADO: {type(exc).__name__}: {exc}"
            blocked += 1
        rows.append(rec)

    print(f"\n=== derive_sites {args.campaign} "
          f"({'WRITE' if args.write else 'DRY-RUN'}) ===")
    for r in rows:
        if r.get("status") == "OK":
            print(f"[OK]  {r['target']:8} {r['site']:22} ref={r['reference_pdb']:5} "
                  f"RMSD={r['rmsd']:.2f} n_ca={r['n_ca']:3} "
                  f"center_target={r['center_target']} size={r['box_size']}")
        else:
            print(f"[!!]  {r['target']:8} {r['site']:22} ref={r.get('reference_pdb')} "
                  f"-> {r['status']}")
    print(f"\nOK: {ok}   BLOQUEADOS: {blocked}   TOTAL: {len(rows)}")
    if args.write:
        print("PERSISTIDO no registry. Rode 'python registry.py' e um redock gate.")
    else:
        print("DRY-RUN: nada foi escrito. Revise antes de --write + redock gate.")


if __name__ == "__main__":
    main()
