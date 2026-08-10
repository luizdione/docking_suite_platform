# -*- coding: utf-8 -*-
"""
Migra o conteudo cientifico de target_selector.py antes de aposenta-lo.

target_selector.py era a SEGUNDA fonte de verdade de alvos: dicionario indexado
por inteiro, com grid_center FIXO e divergente do resto da plataforma. O codigo
dele nao e aproveitavel -- mas o conteudo em portugues (justificativa de escolha
de cada alvo e referencias) e ciencia real e nao pode ser perdido.

O que esta rotina faz:
  1. injeta `notes_pt` e `references` nos alvos ja existentes em registry/targets/;
  2. arquiva o dicionario inteiro em docs/historico/, inclusive os alvos que nunca
     entraram no painel (TcNDH2, LmTR) e os grid_center fixos, marcados como NAO
     confiaveis;
  3. reporta as divergencias encontradas entre as duas fontes.

O que NAO migra, de proposito: os `grid_center` fixos. Para TcTR o valor fixo era
(15,0, -2,0, -18,0) enquanto o calculo sobre o co-cristal da (43,57, 7,66, 0,09).
Sao numeros aproximados que nunca foram usados no docking real.

Uso:
    python tools/migrate_target_selector.py --dry-run
    python tools/migrate_target_selector.py
"""
import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

TARGETS_DIR = ROOT / "registry" / "targets"
ARCHIVE = ROOT / "docs" / "historico" / "target_selector_TARGET_DATABASE.json"


def clean(text):
    """Normaliza a prosa multi-linha do dicionario original."""
    return " ".join(str(text or "").split())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    import registry
    import target_selector

    legacy = target_selector.TARGET_DATABASE
    by_id = {v["id"]: v for v in legacy.values()}

    updated, missing, divergences = [], [], []

    for entry_id, legacy_entry in by_id.items():
        path = TARGETS_DIR / f"{entry_id.lower()}.json"
        if not path.exists():
            missing.append(entry_id)
            continue

        with open(path, encoding="utf-8") as fh:
            target = json.load(fh)

        # divergencia de PDB entre as duas fontes: registrar, nunca sobrescrever
        legacy_pdb = legacy_entry.get("pdb")
        current_pdb = (target.get("structure") or {}).get("pdb_id")
        if legacy_pdb and current_pdb and legacy_pdb != current_pdb:
            divergences.append(
                f"{entry_id}: target_selector dizia {legacy_pdb}, registry diz {current_pdb} "
                "(mantido o do registry, que e o que o docking usou)"
            )

        target["notes_pt"] = clean(legacy_entry.get("justification"))
        refs = legacy_entry.get("references")
        if refs:
            target["references"] = refs if isinstance(refs, list) else [clean(refs)]
        if legacy_entry.get("resolution"):
            target.setdefault("structure", {})["resolution_A"] = legacy_entry["resolution"]
        if legacy_entry.get("target_name_en") and not target.get("name_en"):
            target["name_en"] = clean(legacy_entry["target_name_en"])

        if not args.dry_run:
            with open(path, "w", encoding="utf-8") as fh:
                json.dump(target, fh, indent=2, ensure_ascii=False)
                fh.write("\n")
        updated.append(entry_id)

    archive = {
        "_meta": {
            "source": "target_selector.py TARGET_DATABASE",
            "archived_by": "tools/migrate_target_selector.py",
            "warning": (
                "Arquivo de memoria. Os grid_center aqui sao FIXOS e nao confiaveis: "
                "foram substituidos pelo calculo dinamico de _calculate_grid_center() "
                "sobre o co-cristal. Nao reintroduzir."
            ),
            "targets_never_in_panel": ["TcNDH2", "LmTR"],
        },
        "targets": {str(k): v for k, v in legacy.items()},
    }
    if not args.dry_run:
        ARCHIVE.parent.mkdir(parents=True, exist_ok=True)
        with open(ARCHIVE, "w", encoding="utf-8") as fh:
            json.dump(archive, fh, indent=2, ensure_ascii=False, default=str)
            fh.write("\n")

    print(f"alvos atualizados com notes_pt: {len(updated)} -> {updated}")
    print(f"alvos do target_selector sem correspondente no registry: {missing}")
    if divergences:
        print("\ndivergencias entre as duas fontes (registradas, nao aplicadas):")
        for d in divergences:
            print(f"  - {d}")
    print(f"\narquivo historico: {ARCHIVE.relative_to(ROOT)}")
    if args.dry_run:
        print("\n(dry-run: nada foi escrito)")

    registry.clear_cache()
    errors = registry.validate_all()
    print(f"\nvalidacao do registry apos migracao: {len(errors)} erro(s)")
    for e in errors[:10]:
        print(f"  - {e}")


if __name__ == "__main__":
    main()
