# -*- coding: utf-8 -*-
"""
Leitura e escrita do resultado agregado de docking, com uma convencao so.

A plataforma tinha tres nomes/raizes para o mesmo arquivo, e o sintoma era
"rodei o docking e o dashboard nao mostra nada": o pipeline escrevia num lugar
e a UI lia de outro.

Regra de transicao: **le de todos os lugares, escreve nos dois nomes**. A escrita
dupla cai na Fase 9, quando nenhum leitor antigo restar.

    read_results(campaign, ligand)          -> dict | None
    write_results(campaign, ligand, payload) -> Path canonico
"""
import json
import shutil
from pathlib import Path

import paths

ENCODING = "utf-8"


def read_results(campaign_id=None, ligand=None, run_id=None, strict=False):
    """Le o resultado agregado do primeiro caminho que existir.

    Devolve o dict com a chave extra `_source_path`, para que a UI possa dizer
    de onde leu -- diagnostico que faltava quando as tres convencoes divergiam.
    """
    campaign_id = campaign_id or paths.default_campaign()
    if ligand is None:
        ligand = _default_ligand(campaign_id)

    for path in paths.candidate_results_paths(campaign_id, ligand, run_id):
        if not path.exists():
            continue
        try:
            with open(path, encoding=ENCODING) as fh:
                data = json.load(fh)
        except UnicodeDecodeError:
            # arquivos antigos foram escritos sem encoding explicito no Windows
            with open(path, encoding="cp1252", errors="replace") as fh:
                data = json.load(fh)
        except json.JSONDecodeError as exc:
            if strict:
                raise
            print(f"[AVISO] resultado ilegivel em {path}: {exc}")
            continue
        if isinstance(data, dict):
            data.setdefault("_source_path", str(path))
        return data

    if strict:
        tried = "\n  ".join(str(p) for p in paths.candidate_results_paths(campaign_id, ligand, run_id))
        raise FileNotFoundError(
            f"nenhum resultado para campanha={campaign_id} ligante={ligand}. Procurei em:\n  {tried}"
        )
    return None


def write_results(campaign_id, ligand, payload, run_id=None, write_legacy=True):
    """Grava o resultado no caminho canonico (e no nome legado, na transicao).

    Devolve o caminho canonico.
    """
    campaign_id = campaign_id or paths.default_campaign()
    canonical = paths.results_path(campaign_id, ligand, run_id=run_id, canonical=True)
    canonical.parent.mkdir(parents=True, exist_ok=True)

    if isinstance(payload, dict):
        payload = dict(payload)
        payload.setdefault("campaign", campaign_id)
        payload.setdefault("ligand", ligand)
        if run_id:
            payload.setdefault("run_id", run_id)

    with open(canonical, "w", encoding=ENCODING) as fh:
        json.dump(payload, fh, indent=2, ensure_ascii=False, default=str)
        fh.write("\n")

    if write_legacy:
        legacy = paths.results_path(campaign_id, ligand, run_id=run_id, canonical=False)
        try:
            legacy.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(canonical, legacy)
        except OSError as exc:
            print(f"[AVISO] nao consegui gravar o nome legado {legacy}: {exc}")

    return canonical


def _default_ligand(campaign_id):
    """Ligante default da campanha -- substitui o '2-SCN-NQ' soldado no codigo."""
    try:
        import registry

        campaign = registry.get_campaign(campaign_id)
        declared = campaign.get("default_ligand")
        if declared:
            return declared
        ligands = campaign.get("ligands") or []
        if ligands:
            return ligands[0]
    except Exception:
        pass
    return "2-SCN-NQ"


def describe(campaign_id=None, ligand=None):
    """Diagnostico: onde o resultado esta e onde passaria a ser escrito."""
    campaign_id = campaign_id or paths.default_campaign()
    ligand = ligand or _default_ligand(campaign_id)
    lines = [f"campanha={campaign_id} ligante={ligand}", "candidatos de leitura:"]
    for path in paths.candidate_results_paths(campaign_id, ligand):
        lines.append(f"  {'ACHADO' if path.exists() else '      '}  {path}")
    lines.append(f"escrita canonica : {paths.results_path(campaign_id, ligand)}")
    lines.append(f"escrita legada   : {paths.results_path(campaign_id, ligand, canonical=False)}")
    return "\n".join(lines)


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description="Diagnostico do resultado agregado.")
    ap.add_argument("--campaign")
    ap.add_argument("--ligand")
    args = ap.parse_args()
    print(describe(args.campaign, args.ligand))

    data = read_results(args.campaign, args.ligand)
    if data is None:
        print("\nnenhum resultado encontrado")
    else:
        keys = [k for k in data if not k.startswith("_")]
        print(f"\nlido de: {data.get('_source_path')}")
        print(f"chaves : {keys[:12]}")
