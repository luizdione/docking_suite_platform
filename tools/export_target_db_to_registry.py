# -*- coding: utf-8 -*-
"""
Converte docking_v2.TARGET_DB nos arquivos declarativos de registry/targets/.

Gerado, nunca transcrito a mao: transcricao de 10 alvos x 12 campos e onde
nascem os erros silenciosos -- foi assim que o fallback de alvos do dashboard
acabou com tres PDBs errados em relacao ao TARGET_DB.

Uso:
    python tools/export_target_db_to_registry.py --dry-run
    python tools/export_target_db_to_registry.py
"""
import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

OUT_DIR = ROOT / "registry" / "targets"

# O TARGET_DB nao registra o nome do sitio; derivamos um rotulo estavel a partir
# da estrategia, para que o par possa referenciar o sitio por nome.
SITE_NAME_BY_STRATEGY = {
    "ligand_center": "orthosteric",
    "metal_center": "metal_active_site",
    "residue_center": "catalytic_site",
}


def build_target(tid, info):
    strategy = info.get("grid_strategy", "ligand_center")
    site_name = SITE_NAME_BY_STRATEGY.get(strategy, "primary")

    site = {
        "definition": strategy,
        "box_size_A": list(info.get("grid_size", (22, 22, 22))),
        "box_policy": "fixed",
    }
    if info.get("ligand_resnames"):
        site["ligand_resnames"] = list(info["ligand_resnames"])
    if info.get("grid_fallback_residues"):
        site["fallback_residues"] = list(info["grid_fallback_residues"])
    if info.get("pdb"):
        # A numeracao de residuo do sitio so faz sentido em relacao a uma
        # estrutura; sem isso, transferir sitio entre ortologos vira chute.
        site["residue_numbering_reference"] = info["pdb"]

    return {
        "schema_version": "2.0",
        "id": tid,
        "aliases": [info["pdb"]] if info.get("pdb") else [],
        "name": info.get("target"),
        "organism": info.get("organism"),
        "category": info.get("category"),
        "ortholog_of": info.get("human_ortholog"),
        "structure": {
            "source": "pdb",
            "pdb_id": info.get("pdb"),
            "assembly": None,
            "n_protomers_required": 1,
            "local_file": None,
            "confidence": None,
        },
        "chemistry": {
            "keep_cofactors": list(info.get("keep_cofactors") or []),
            "essential_cofactors": list(info.get("keep_cofactors") or []),
            "cofactor_source": "crystal" if info.get("keep_cofactors") else "none",
            "nonstandard_residues": [],
            "protonation_pH": 7.4,
        },
        "sites": {site_name: site},
        "default_site": site_name,
        "notes": info.get("notes", ""),
        "literature": info.get("literature", ""),
        "provenance": {
            "generated_from": "docking_v2.TARGET_DB",
            "generator": "tools/export_target_db_to_registry.py",
        },
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--force", action="store_true", help="sobrescreve alvos ja existentes")
    args = ap.parse_args()

    import docking_v2

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    written, skipped = 0, 0

    for tid, info in docking_v2.TARGET_DB.items():
        payload = build_target(tid, info)
        path = OUT_DIR / f"{tid.lower()}.json"

        if path.exists() and not args.force:
            print(f"[skip] {path.name} ja existe (use --force)")
            skipped += 1
            continue
        if args.dry_run:
            print(f"[dry] {path.name}")
            print(json.dumps(payload, indent=2, ensure_ascii=False)[:400] + " ...")
            continue

        with open(path, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2, ensure_ascii=False, sort_keys=False)
            fh.write("\n")
        print(f"[ok] {path.name}")
        written += 1

    print(f"\n{written} escrito(s), {skipped} pulado(s)")


if __name__ == "__main__":
    main()
