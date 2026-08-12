# -*- coding: utf-8 -*-
"""tools/validate_oligomer_sites.py -- valida (sem docking, sem GPU) os sitios de
receptores OLIGOMERICOS cujo ligante ja esta EMBUTIDO na pose nativa (ex.: pentamero
de F-actina com jasplaquinolida 9UE). Diferente de validate_derived_sites.py, aqui NAO
ha superposicao: o ligante ja esta no frame do receptor e o sitio e INTER-PROTOMERICO
(multi-cadeia). Usa site_derivation.validate_embedded_ligand, que escolhe a copia mais
imersa (central) do ligante.

Dupla funcao:
  * RELATORIO honesto do bolso (imersao da pose nativa) -- o contraste com o monomero
    (imersao fraca) e o resultado 'estado oligomerico como determinante'.
  * PARITY: reconfere que o box_center_A persistido bate com o recomputado (gate de
    regressao -- se o modelo mudar em disco sem o registry ser atualizado, quebra alto).

Uso:
    python tools/validate_oligomer_sites.py                 # relatorio + parity
    python tools/validate_oligomer_sites.py --campaign SUB04
    python tools/validate_oligomer_sites.py --tol 0.05      # tolerancia da parity (A)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import registry            # noqa: E402
import site_derivation     # noqa: E402


def _oligomer_targets(campaign_id=None):
    ids = registry.list_targets(campaign_id) if campaign_id else \
        [registry.get_target(t) for t in registry._ids("targets")]
    for t in ids:
        struct = t.get("structure") or {}
        if (struct.get("n_protomers_required") or 1) > 1:
            yield t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--campaign", default=None)
    ap.add_argument("--tol", type=float, default=0.05, help="tolerancia da parity do centro (A)")
    args = ap.parse_args()

    rows, mismatches = [], []
    print(f"{'ALVO':22}{'SITIO':16}{'copia':7}{'imersao':>8}{'contato':>8}  veredito")
    for t in _oligomer_targets(args.campaign):
        tid = t.get("id")
        local = (t.get("structure") or {}).get("local_file")
        if not local:
            continue  # source=pdb (ex.: tubulina 1SA0) -- valida no docking, nao aqui
        cif = registry.resolve_local_file(local)
        if not cif or not Path(cif).exists():
            continue
        for sid, site in (t.get("sites") or {}).items():
            resname = site.get("reference_ligand_resname") or (
                (site.get("ligand_resnames") or [None])[0])
            if not resname:
                continue
            try:
                r = site_derivation.validate_embedded_ligand(cif, resname)
            except Exception as exc:  # noqa: BLE001
                print(f"{tid:22}{sid:16}  ERRO: {type(exc).__name__}: {str(exc)[:50]}")
                mismatches.append((tid, sid, "erro"))
                continue
            copy = f"{r['copy_chain']}{r['copy_resseq']}"
            print(f"{tid:22}{sid:16}{copy:7}{r['immersion_mean_A']:>7}A{r['frac_contact']:>8}  {r['verdict']}")
            rows.append((tid, sid, r))
            # parity contra o persistido
            stored = site.get("box_center_A")
            if stored is not None:
                d = max(abs(a - b) for a, b in zip(stored, r["box_center_A"]))
                if d > args.tol:
                    mismatches.append((tid, sid, f"centro difere {d:.3f}A > tol {args.tol}A"))

    n_ok = sum(1 for _, _, r in rows if r["verdict"] == "center_ok")
    print(f"\ncenter_ok: {n_ok}   total sitios oligomericos: {len(rows)}")
    if mismatches:
        print("\n[FALHA] parity/erros:")
        for tid, sid, msg in mismatches:
            print(f"  - {tid} / {sid}: {msg}")
        sys.exit(1)
    print("[ok] parity: box_center_A persistido bate com o recomputado.")


if __name__ == "__main__":
    main()
