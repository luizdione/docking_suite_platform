# -*- coding: utf-8 -*-
"""
Resolucao de caminhos por campanha.

Antes disso conviviam tres convencoes para "onde esta o resultado do docking":

    dashboard_v2  -> SCRIPT_DIR/docking_project/results/docking_results.json
                     (ignorando NQ_DATA_DIR, o proprio env var do dashboard)
    /data_status  -> BASE_DIR/results/<composto>/results.json
    analysis.py   -> docking_project/results/results.json

Tres nomes, tres raizes, o mesmo dado. Este modulo estabelece uma convencao e
resolve as antigas em cascata, para que os dados do SUB01 continuem sendo lidos
exatamente de onde ja estao -- sem mover um byte.

Layout canonico:

    <data_root>/<CAMPANHA>/
        proteins/                       receptores (deixou de ser global)
        ligand/<ligante>/
        runs/<run_id>/                  unidade de execucao, com proveniencia
        results/<ligante>/results.json  convencao UNICA
        figures/<ligante>/<sessao>/
        logs/<ligante>/<sessao>/
        manuscript_texts/<ligante>/<sessao>/

Precedencia da raiz de dados de uma campanha:
    1. campaign.paths.data_root declarado no registry
    2. $DOCKING_DATA_DIR/<CAMPANHA>, se existir
    3. $NQ_DATA_DIR             (layout legado -- so para a campanha default)
    4. <repo>/docking_project   (fallback historico)

Variaveis de ambiente: DOCKING_DATA_DIR (nova), NQ_DATA_DIR (legada).
"""
import os
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parent

# Subdiretorios da arvore de uma campanha. Nomes iguais aos ja usados por
# compound_db.get_session_dirs(), de proposito: a arvore legada continua valida.
SUBDIRS = (
    "proteins",
    "ligand",
    "runs",
    "results",
    "figures",
    "logs",
    "manuscript_texts",
)

# Nome canonico do arquivo de resultado agregado de docking.
RESULTS_FILENAME = "results.json"

# Nome legado, escrito por docking_v2. Continua sendo lido; deixa de ser escrito
# na Fase 9.
LEGACY_RESULTS_FILENAME = "docking_results.json"


def _env_path(name):
    value = os.environ.get(name, "").strip()
    return Path(value) if value else None


def default_campaign():
    """Campanha corrente: env DOCKING_CAMPAIGN, senao o default do programa."""
    env = os.environ.get("DOCKING_CAMPAIGN", "").strip()
    if env:
        return env
    try:
        import registry

        return registry.load_program().get("default_campaign", "SUB01")
    except Exception:
        return "SUB01"


def _campaign_declared_root(campaign_id):
    """data_root declarado no registry da campanha, se houver."""
    try:
        import registry

        campaign = registry.get_campaign(campaign_id)
    except Exception:
        return None
    declared = (campaign.get("paths") or {}).get("data_root")
    return Path(declared) if declared else None


def data_root(campaign_id=None, warn_legacy=True):
    """Raiz de dados de uma campanha, resolvida em cascata.

    Nao cria diretorio: resolver caminho e efeito colateral demais para uma
    funcao que o dashboard chama em toda requisicao.
    """
    campaign_id = campaign_id or default_campaign()

    declared = _campaign_declared_root(campaign_id)
    if declared:
        return declared

    new_root = _env_path("DOCKING_DATA_DIR")
    if new_root:
        candidate = new_root / campaign_id
        if candidate.exists():
            return candidate
        # raiz nova declarada mas ainda sem a arvore da campanha: usamos assim
        # mesmo, para que a primeira escrita ja nasca no layout novo.
        if campaign_id != default_campaign():
            return candidate

    legacy = _env_path("NQ_DATA_DIR")
    if legacy:
        if campaign_id == default_campaign():
            # SUB01 vive na arvore antiga (E:/data/docking). Ler dali e o certo.
            return legacy
        candidate = legacy / campaign_id
        if candidate.exists():
            return candidate
        if warn_legacy:
            warnings.warn(
                f"campanha '{campaign_id}' sem raiz de dados propria; declare "
                "paths.data_root no registry dela ou defina DOCKING_DATA_DIR. "
                f"Usando {candidate} (pode nao existir).",
                RuntimeWarning,
                stacklevel=2,
            )
        return candidate

    fallback = ROOT / "docking_project"
    if campaign_id == default_campaign():
        return fallback
    return fallback / campaign_id


def campaign_dirs(campaign_id=None, create=False):
    """Mapa subdiretorio -> caminho para uma campanha."""
    base = data_root(campaign_id)
    dirs = {name: base / name for name in SUBDIRS}
    dirs["root"] = base
    if create:
        for name, path in dirs.items():
            if name != "root":
                path.mkdir(parents=True, exist_ok=True)
    return dirs


def run_dir(campaign_id, run_id, create=False):
    """Diretorio de uma execucao. E onde vive a proveniencia (run.json)."""
    path = data_root(campaign_id) / "runs" / str(run_id)
    if create:
        path.mkdir(parents=True, exist_ok=True)
    return path


def ligand_dirname(ligand):
    """Nome de diretorio de um ligante, no case historico.

    O id do registry e minusculo por estabilidade entre sistemas de arquivos
    ('2-scn-nq'), mas os diretorios de dados ja existentes preservam o case
    original ('2-SCN-NQ', '2-Br-NQ'). Traduzimos id -> name antes de tocar o
    disco: sem isso, o layout novo criaria uma arvore paralela ao lado da
    antiga -- invisivel no Windows, que e case-insensitive, e duplicada no
    Linux.
    """
    name = str(ligand)
    try:
        import registry

        entry = registry.get_ligand(name)
        name = entry.get("name") or name
    except Exception:
        pass
    try:
        import compound_db

        return compound_db._safe_name(name)
    except Exception:
        return name


# alias interno, mantido por legibilidade nas chamadas
_safe_ligand = ligand_dirname


def results_path(campaign_id, ligand, run_id=None, canonical=True):
    """Caminho do resultado agregado.

    Com run_id, o resultado da execucao especifica; sem, o cumulativo do ligante
    (que e como o dashboard sempre leu).
    """
    name = RESULTS_FILENAME if canonical else LEGACY_RESULTS_FILENAME
    if run_id:
        return run_dir(campaign_id, run_id) / name
    return data_root(campaign_id) / "results" / _safe_ligand(ligand) / name


def candidate_results_paths(campaign_id, ligand, run_id=None):
    """Todos os lugares onde o resultado pode estar, em ordem de preferencia.

    Existe para que a leitura funcione tanto no layout novo quanto nos tres
    layouts historicos, sem exigir migracao previa dos dados.
    """
    base = data_root(campaign_id)
    safe = _safe_ligand(ligand)
    candidates = []

    if run_id:
        rd = run_dir(campaign_id, run_id)
        candidates += [rd / RESULTS_FILENAME, rd / LEGACY_RESULTS_FILENAME]

    candidates += [
        base / "results" / safe / RESULTS_FILENAME,
        base / "results" / safe / LEGACY_RESULTS_FILENAME,
        base / "results" / RESULTS_FILENAME,
        base / "results" / LEGACY_RESULTS_FILENAME,
        # convencao que docking_v2/dashboard usavam, relativa ao repo e nao ao
        # diretorio de dados -- a origem do bug de "rodei mas nao aparece".
        ROOT / "docking_project" / "results" / safe / RESULTS_FILENAME,
        ROOT / "docking_project" / "results" / LEGACY_RESULTS_FILENAME,
        ROOT / "docking_project" / "results" / RESULTS_FILENAME,
    ]

    seen, ordered = set(), []
    for path in candidates:
        key = str(path).lower()
        if key not in seen:
            seen.add(key)
            ordered.append(path)
    return ordered


def latest_results_path(campaign_id, ligand, run_id=None):
    """Primeiro candidato que existe em disco, ou None."""
    for path in candidate_results_paths(campaign_id, ligand, run_id):
        if path.exists():
            return path
    return None


def session_dirs(campaign_id, ligand, session_tag, create=True):
    """Delegado a compound_db.get_session_dirs() com a raiz da campanha.

    A logica de sessao por timestamp ja estava certa; o que faltava era o
    prefixo de campanha. Reimplementar criaria uma quarta convencao.
    """
    import compound_db

    return compound_db.get_session_dirs(data_root(campaign_id), session_tag, ligand)


def resolve(template, campaign_id=None):
    """Expande ${DATA}, ${REG}, ${CODE} e variaveis de ambiente num caminho.

    Permite que o registry declare caminhos portateis em vez de absolutos --
    os caminhos absolutos espalhados pelos scripts sao o que hoje impede rodar
    a mesma campanha em outra maquina.
    """
    if template is None:
        return None
    campaign_id = campaign_id or default_campaign()
    text = str(template)

    try:
        import registry as _registry

        reg_dir = str(_registry.REGISTRY_DIR)
    except Exception:
        reg_dir = str(ROOT / "registry")

    code_root = str(ROOT)
    try:
        import registry as _registry

        declared = (_registry.get_campaign(campaign_id).get("paths") or {}).get("code_root")
        if declared:
            code_root = str(declared)
    except Exception:
        pass

    for key, value in (
        ("${DATA}", str(data_root(campaign_id))),
        ("${REG}", reg_dir),
        ("${CODE}", code_root),
        ("${CAMPAIGN}", campaign_id),
    ):
        text = text.replace(key, value)

    return Path(os.path.expandvars(text)).expanduser()


def describe(campaign_id=None):
    """Diagnostico legivel de onde a campanha esta lendo e escrevendo."""
    campaign_id = campaign_id or default_campaign()
    base = data_root(campaign_id)
    lines = [
        f"campanha        : {campaign_id}",
        f"data_root       : {base}  {'(existe)' if base.exists() else '(NAO existe)'}",
        f"DOCKING_DATA_DIR: {os.environ.get('DOCKING_DATA_DIR') or '(nao definido)'}",
        f"NQ_DATA_DIR     : {os.environ.get('NQ_DATA_DIR') or '(nao definido)'}",
    ]
    for name in SUBDIRS:
        path = base / name
        lines.append(f"  {name:18s} {'ok ' if path.exists() else '-- '} {path}")
    return "\n".join(lines)


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description="Diagnostico de caminhos por campanha.")
    ap.add_argument("--campaign")
    ap.add_argument("--ligand", default=None, help="mostra onde o resultado seria procurado")
    args = ap.parse_args()

    campaign = args.campaign or default_campaign()
    print(describe(campaign))

    if args.ligand:
        print(f"\ncandidatos de resultado para '{args.ligand}':")
        for path in candidate_results_paths(campaign, args.ligand):
            print(f"  {'ACHADO' if path.exists() else '      '}  {path}")
