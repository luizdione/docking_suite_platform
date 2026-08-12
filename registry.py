# -*- coding: utf-8 -*-
"""
Registry declarativo da plataforma (guarda-chuva multi-alvo).

Fonte unica de verdade sobre PROGRAMA / CAMPANHA / ALVO / LIGANTE / PAR /
PROTOCOLO. Codigo le o registry; codigo nao contem dicionario de dominio.

Hierarquia:

    ALVO e LIGANTE sao globais e reutilizaveis entre campanhas.
    O PAR (alvo x ligante, escopado por campanha) e onde vive a peculiaridade
    -- estado oligomerico exigido pelo sitio, Kd experimental com o estado em
    que foi medido, qual estimador e confiavel, vazamento de treino do preditor.

    Foi a ausencia desse conceito que permitiu, por exemplo, md_simulation.py
    reusar cofactor_FAD_TcTR.xml para alvos que nao sao flavoenzimas.

Layout em disco (JSON; ver docs/REGISTRY.md):

    registry/
      program.json
      campaigns/<ID>.json
      targets/<id>.json
      ligands/<id>.json
      pairs/<CAMPANHA>__<alvo>__<ligante>.json
      protocols/<id>.json
      _schema/*.schema.json

Variavel de ambiente: DOCKING_REGISTRY sobrepoe o diretorio.
"""
import hashlib
import json
import os
import threading
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REGISTRY_DIR = Path(os.environ.get("DOCKING_REGISTRY", ROOT / "registry"))

# I/O sempre em UTF-8 explicito: os campos `notes` carregam prosa cientifica
# acentuada em portugues e o default do Windows (cp1252) corrompe silenciosamente.
ENCODING = "utf-8"

_LOCK = threading.RLock()
_CACHE = {}          # path -> (mtime_ns, payload)
_ALIAS_CACHE = None  # alias -> id canonico de alvo


class RegistryError(Exception):
    """Erro de leitura, validacao ou resolucao no registry."""


# ---------------------------------------------------------------- leitura base

def _read(path):
    """Le um JSON do registry com cache invalidado por mtime."""
    path = Path(path)
    if not path.exists():
        raise RegistryError(f"arquivo nao encontrado no registry: {path}")
    with _LOCK:
        try:
            mtime = path.stat().st_mtime_ns
        except OSError as exc:
            raise RegistryError(f"nao foi possivel ler {path}: {exc}") from exc
        cached = _CACHE.get(path)
        if cached and cached[0] == mtime:
            return cached[1]
        try:
            with open(path, encoding=ENCODING) as fh:
                payload = json.load(fh)
        except json.JSONDecodeError as exc:
            raise RegistryError(f"JSON invalido em {path}: {exc}") from exc
        except UnicodeDecodeError as exc:
            raise RegistryError(
                f"{path} nao esta em UTF-8: {exc}. Todo arquivo do registry deve ser UTF-8."
            ) from exc
        _CACHE[path] = (mtime, payload)
        return payload


def clear_cache():
    """Descarta o cache. Util em testes e apos escrita programatica."""
    global _ALIAS_CACHE
    with _LOCK:
        _CACHE.clear()
        _ALIAS_CACHE = None


def _dir(kind):
    return REGISTRY_DIR / kind


def _ids(kind):
    d = _dir(kind)
    if not d.exists():
        return []
    # Ignora sidecars com extensao composta (ex.: SUB01.paper.json -- dados
    # congelados do manuscrito consumidos por export_web.build_viewer, NAO uma
    # campanha). p.stem de "SUB01.paper.json" seria "SUB01.paper".
    return sorted(
        p.stem for p in d.glob("*.json")
        if not p.stem.startswith("_") and not p.name.endswith(".paper.json")
    )


# ---------------------------------------------------------------- programa

def load_program():
    """Metadados do guarda-chuva (nome, responsavel, ORCID, defaults quimicos)."""
    return _read(REGISTRY_DIR / "program.json")


def program_defaults():
    """Defaults quimicos do programa, sobrepostos por target.chemistry.*.

    Herdeiros de ESSENTIAL_COFACTORS / ESSENTIAL_METALS / REMOVE_RESIDUES, que
    viviam em docking_v2.py enviesados para flavoenzimas.
    """
    return load_program().get("chemistry_defaults", {})


# ---------------------------------------------------------------- campanhas

def list_campaigns():
    return _ids("campaigns")


def get_campaign(campaign_id):
    try:
        return _read(_dir("campaigns") / f"{campaign_id}.json")
    except RegistryError:
        raise RegistryError(
            f"campanha '{campaign_id}' nao existe. Disponiveis: {list_campaigns()}"
        ) from None


# ---------------------------------------------------------------- alvos

def _alias_map():
    """alias (case-insensitive) -> id canonico."""
    global _ALIAS_CACHE
    with _LOCK:
        if _ALIAS_CACHE is not None:
            return _ALIAS_CACHE
        mapping = {}
        for tid in _ids("targets"):
            data = _read(_dir("targets") / f"{tid}.json")
            canonical = data.get("id", tid)
            mapping[canonical.lower()] = canonical
            mapping[tid.lower()] = canonical
            for alias in data.get("aliases", []):
                mapping.setdefault(str(alias).lower(), canonical)
        _ALIAS_CACHE = mapping
        return mapping


def _target_path(target_id):
    """Resolve id/alias -> caminho do arquivo, tolerando maiusculas."""
    direct = _dir("targets") / f"{target_id}.json"
    if direct.exists():
        return direct
    lowered = _dir("targets") / f"{str(target_id).lower()}.json"
    if lowered.exists():
        return lowered
    canonical = _alias_map().get(str(target_id).lower())
    if canonical:
        for cand in (canonical, canonical.lower()):
            p = _dir("targets") / f"{cand}.json"
            if p.exists():
                return p
    raise RegistryError(
        f"alvo '{target_id}' nao existe no registry. Disponiveis: {_ids('targets')}"
    )


def get_target(target_id):
    """Alvo por id ou alias. Nunca devolve None: erra alto."""
    return _read(_target_path(target_id))


def resolve_local_file(raw):
    """Resolve structure.local_file para um caminho utilizavel em qualquer CWD.

    Aceita ABSOLUTO (usado como esta -- back-compat com maquinas que apontam p/
    fora do repo) OU RELATIVO a raiz do repo (ex.: 'docking_project/SUB04/models/
    X.cif'), que e portavel: resolve local e dentro do container (/app) igual.
    Nao levanta se o arquivo nao existir -- quem consome/valida decide."""
    if not raw:
        return raw
    p = Path(raw)
    if p.is_absolute():
        return str(p)
    return str((ROOT / raw).resolve())


def list_targets(campaign_id=None):
    """Alvos do registry, ou apenas os da campanha (na ordem declarada nela)."""
    if campaign_id is None:
        return [get_target(t) for t in _ids("targets")]
    campaign = get_campaign(campaign_id)
    return [get_target(t) for t in campaign.get("targets", [])]


def resolve_site(target_id, site_id=None):
    """Sitio de ligacao de um alvo, com o id resolvido dentro do dict.

    Sem site_id, devolve o `default_site`. Se o alvo tiver um unico sitio, ele
    e o default implicito.
    """
    target = get_target(target_id)
    sites = target.get("sites") or {}
    if not sites:
        raise RegistryError(f"alvo '{target_id}' nao declara nenhum sitio")

    if site_id is None:
        site_id = target.get("default_site")
        if site_id is None:
            if len(sites) == 1:
                site_id = next(iter(sites))
            else:
                raise RegistryError(
                    f"alvo '{target_id}' tem {len(sites)} sitios e nenhum default_site; "
                    f"escolha um: {sorted(sites)}"
                )
    if site_id not in sites:
        raise RegistryError(
            f"sitio '{site_id}' nao existe em '{target_id}'. Disponiveis: {sorted(sites)}"
        )
    site = dict(sites[site_id])
    site["id"] = site_id
    site["target"] = target.get("id", target_id)
    return site


# ---------------------------------------------------------------- ligantes

def get_ligand(ligand_id):
    path = _dir("ligands") / f"{ligand_id}.json"
    if not path.exists():
        path = _dir("ligands") / f"{_safe_id(ligand_id)}.json"
    if not path.exists():
        raise RegistryError(
            f"ligante '{ligand_id}' nao existe no registry. Disponiveis: {_ids('ligands')}"
        )
    return _read(path)


def list_ligands(campaign_id=None):
    if campaign_id is None:
        return [_read(_dir("ligands") / f"{lid}.json") for lid in _ids("ligands")]
    campaign = get_campaign(campaign_id)
    return [get_ligand(l) for l in campaign.get("ligands", [])]


def _safe_id(value):
    """Normaliza um nome para id de arquivo do registry.

    Difere de compound_db._safe_name() de proposito: ali o nome preserva
    maiusculas porque vira diretorio de dados ja existente em disco ('2-SCN-NQ');
    aqui o id e minusculo e com hifen para que o nome do arquivo seja estavel
    entre sistemas de arquivos case-insensitive. get_ligand() aceita as duas formas.
    """
    out = []
    for ch in str(value).strip().lower():
        out.append(ch if (ch.isalnum() or ch in "-_") else "-")
    return "".join(out).strip("-") or "unnamed"


# ---------------------------------------------------------------- pares

def pair_id(campaign_id, target_id, ligand_id):
    return f"{campaign_id}__{_safe_id(target_id)}__{_safe_id(ligand_id)}"


def get_pair(campaign_id, target_id, ligand_id, strict=True):
    """Par alvo-ligante de uma campanha.

    Com strict=False devolve None quando o par nao foi declarado -- util para
    percorrer combinacoes onde nem todo par faz sentido cientifico (ex.: um
    ligante que so tem sitio no filamento nao se pareia com um alvo monomerico).
    """
    path = _dir("pairs") / f"{pair_id(campaign_id, target_id, ligand_id)}.json"
    if not path.exists():
        if strict:
            raise RegistryError(
                f"par nao declarado: campanha={campaign_id} alvo={target_id} "
                f"ligante={ligand_id}. Crie {path.name} em registry/pairs/."
            )
        return None
    return _read(path)


def list_pairs(campaign_id=None):
    pairs = []
    for pid in _ids("pairs"):
        data = _read(_dir("pairs") / f"{pid}.json")
        if campaign_id is None or data.get("campaign") == campaign_id:
            pairs.append(data)
    return pairs


# ---------------------------------------------------------------- protocolos

def get_protocol(protocol_id):
    path = _dir("protocols") / f"{protocol_id}.json"
    if not path.exists():
        raise RegistryError(
            f"protocolo '{protocol_id}' nao existe. Disponiveis: {_ids('protocols')}"
        )
    return _read(path)


def list_protocols(engine=None):
    out = []
    for pid in _ids("protocols"):
        data = _read(_dir("protocols") / f"{pid}.json")
        if engine is None or data.get("engine") == engine:
            out.append(data)
    return out


def resolve_protocol(campaign_id, kind, target_id=None, ligand_id=None):
    """Protocolo efetivo para uma etapa, com o par tendo precedencia.

    kind: "docking" | "structure" | "md" | ...
    Ordem: par -> campanha. Sem isso, o protocolo de MD do SUB01 (ff14SB/HMR/300K)
    vazaria para uma campanha que exige CHARMM36m/2fs/310K.
    """
    if target_id and ligand_id:
        pair = get_pair(campaign_id, target_id, ligand_id, strict=False)
        if pair:
            pid = (pair.get("protocols") or {}).get(kind)
            if pid:
                return get_protocol(pid)
    campaign = get_campaign(campaign_id)
    pid = (campaign.get("protocols") or {}).get(kind)
    if not pid:
        raise RegistryError(
            f"campanha '{campaign_id}' nao declara protocolo de '{kind}' "
            f"(nem o par {target_id}/{ligand_id})"
        )
    return get_protocol(pid)


# ---------------------------------------------------------------- ponte legada

# Chaves que docking_v2.TARGET_DB expunha. A ponte tem que devolver exatamente
# estas, com os mesmos tipos -- inclusive grid_size como TUPLA (uma lista
# desempacota igual e a regressao passaria despercebida em runtime).
_LEGACY_FIELDS = (
    "pdb",
    "organism",
    "target",
    "category",
    "keep_cofactors",
    "ligand_resnames",
    "grid_strategy",
    "grid_fallback_residues",
    "grid_size",
    "human_ortholog",
    "notes",
    "literature",
)


def as_legacy_target_db(campaign_id="SUB01", site_id=None):
    """Projeta o registry na forma do antigo docking_v2.TARGET_DB.

    Ponte da Fase 2: docking_v2 passa a consumir o registry sem que nenhuma
    outra linha dele mude. O campo extra `id` e incluido para permitir eliminar
    o lookup reverso por PDB em clean_protein(), que estoura com alvo sem PDB
    (modelo AlphaFold) e fica ambiguo com PDB repetido.
    """
    legacy = {}
    for target in list_targets(campaign_id):
        tid = target["id"]
        site = resolve_site(tid, site_id) if target.get("sites") else {}
        chem = target.get("chemistry") or {}
        struct = target.get("structure") or {}

        entry = {
            "id": tid,
            "pdb": struct.get("pdb_id"),
            "organism": target.get("organism"),
            "target": target.get("name"),
            "category": target.get("category"),
            "keep_cofactors": list(chem.get("keep_cofactors") or []),
            "ligand_resnames": list(site.get("ligand_resnames") or []),
            "grid_strategy": site.get("definition"),
            "grid_fallback_residues": list(site.get("fallback_residues") or []),
            "grid_size": tuple(site.get("box_size_A") or ()),
            "human_ortholog": target.get("ortholog_of"),
            "notes": target.get("notes", ""),
            "literature": target.get("literature", ""),
        }
        legacy[tid] = {k: entry[k] for k in _LEGACY_FIELDS if k in entry} | {"id": tid}
    return legacy


# ---------------------------------------------------------------- validacao

def validate_all():
    """Valida o registry inteiro. Devolve lista de erros ([] = tudo certo).

    Validador proprio e deliberado: jsonschema nao esta em requirements.txt e o
    ambiente Windows/conda ja e fragil. Se jsonschema estiver disponivel, os
    arquivos de _schema/ sao usados como checagem adicional.
    """
    errors = []

    if not REGISTRY_DIR.exists():
        return [f"diretorio do registry nao existe: {REGISTRY_DIR}"]

    clear_cache()

    # --- programa
    try:
        load_program()
    except RegistryError as exc:
        errors.append(f"program.json: {exc}")

    known_targets, known_ligands, known_protocols = set(), set(), set()

    # --- alvos
    for tid in _ids("targets"):
        try:
            t = _read(_dir("targets") / f"{tid}.json")
        except RegistryError as exc:
            errors.append(str(exc))
            continue
        ctx = f"targets/{tid}.json"
        canonical = t.get("id")
        if not canonical:
            errors.append(f"{ctx}: campo obrigatorio 'id' ausente")
            continue
        known_targets.add(canonical)

        for field in ("name", "organism", "category"):
            if not t.get(field):
                errors.append(f"{ctx}: campo obrigatorio '{field}' ausente")

        struct = t.get("structure") or {}
        source = struct.get("source")
        if source not in ("pdb", "alphafold", "local"):
            errors.append(f"{ctx}: structure.source deve ser pdb|alphafold|local (got {source!r})")
        if source == "pdb" and not struct.get("pdb_id"):
            errors.append(f"{ctx}: structure.source='pdb' exige structure.pdb_id")
        if source in ("alphafold", "local") and not struct.get("local_file"):
            errors.append(f"{ctx}: structure.source={source!r} exige structure.local_file")
        elif source in ("alphafold", "local") and struct.get("local_file"):
            # local_file aponta pra fora do repo (dado ingerido) -- se o projeto-fonte
            # for reorganizado, isto quebra em silencio sem esta checagem (ja aconteceu
            # uma vez: Rho_OCRL/.../structures/ -> colaboracoes/tcactina/.../models/).
            if not Path(resolve_local_file(struct["local_file"])).exists():
                errors.append(
                    f"{ctx}: structure.local_file nao existe em disco: "
                    f"{struct['local_file']!r} (resolvido: "
                    f"{resolve_local_file(struct['local_file'])!r})"
                )
        if source == "alphafold" and not (struct.get("confidence") or {}).get("global"):
            errors.append(
                f"{ctx}: modelo AlphaFold sem confidence.global -- pLDDT e obrigatorio "
                "para modelo predito (o revisor vai perguntar)"
            )

        # Estado oligomerico do RECEPTOR (ex.: pentamero de F-actina). Ate agora
        # structure.assembly / n_protomers_required existiam nos JSON mas NUNCA eram
        # validados nem lidos (campos mortos). Um receptor oligomerico que nao declara
        # quantos protomeros tem esconde do runner/relatorio a natureza do modelo.
        n_prot = struct.get("n_protomers_required")
        if n_prot is not None:
            if not isinstance(n_prot, int) or isinstance(n_prot, bool) or n_prot < 1:
                errors.append(
                    f"{ctx}: structure.n_protomers_required deve ser int >= 1 (got {n_prot!r})")
            elif n_prot > 1 and not struct.get("assembly"):
                errors.append(
                    f"{ctx}: n_protomers_required={n_prot} exige structure.assembly descrito "
                    "(ex.: 'F-actin filament (5-mer)') -- nao deixar o estado oligomerico implicito")
        if struct.get("assembly") is not None and not isinstance(struct.get("assembly"), str):
            errors.append(
                f"{ctx}: structure.assembly deve ser string|null "
                f"(got {type(struct.get('assembly')).__name__})")

        chem = t.get("chemistry") or {}
        cof_source = chem.get("cofactor_source")
        if cof_source not in (None, "crystal", "transferred", "none"):
            errors.append(f"{ctx}: chemistry.cofactor_source invalido: {cof_source!r}")
        if cof_source == "transferred":
            tr = chem.get("cofactor_transfer") or {}
            for field in ("method", "donor_pdb", "align_scope"):
                if not tr.get(field):
                    errors.append(f"{ctx}: cofactor_source='transferred' exige cofactor_transfer.{field}")
            if tr.get("align_scope") not in (None, "global", "site_pocket"):
                errors.append(f"{ctx}: cofactor_transfer.align_scope invalido: {tr.get('align_scope')!r}")

        sites = t.get("sites") or {}
        if not sites:
            errors.append(f"{ctx}: nenhum sitio declarado")
        for sid, site in sites.items():
            sctx = f"{ctx} sitio '{sid}'"
            definition = site.get("definition")
            valid = (
                "ligand_center",
                "metal_center",
                "residue_center",
                "derived_from_reference_complex",
                "explicit_center",
            )
            if definition not in valid:
                errors.append(f"{sctx}: definition invalido {definition!r}; esperado um de {valid}")
            if definition == "ligand_center" and not site.get("ligand_resnames"):
                errors.append(f"{sctx}: definition='ligand_center' exige ligand_resnames")
            if definition == "residue_center" and not site.get("fallback_residues"):
                errors.append(f"{sctx}: definition='residue_center' exige fallback_residues")
            if definition == "derived_from_reference_complex":
                for field in ("reference_pdb", "reference_ligand_resname"):
                    if not site.get(field):
                        errors.append(f"{sctx}: definition='derived_from_reference_complex' exige {field}")
            if definition == "explicit_center" and not site.get("box_center_A"):
                errors.append(f"{sctx}: definition='explicit_center' exige box_center_A")

            box = site.get("box_size_A")
            if not (isinstance(box, (list, tuple)) and len(box) == 3):
                errors.append(f"{sctx}: box_size_A deve ser [x,y,z] (got {box!r})")
            center = site.get("box_center_A")
            if center is not None and not (isinstance(center, (list, tuple)) and len(center) == 3):
                errors.append(f"{sctx}: box_center_A deve ser [x,y,z] (got {center!r})")

        default_site = t.get("default_site")
        if default_site and default_site not in sites:
            errors.append(f"{ctx}: default_site '{default_site}' nao existe entre {sorted(sites)}")
        if not default_site and len(sites) > 1:
            errors.append(f"{ctx}: {len(sites)} sitios exigem default_site explicito")

    # --- ligantes
    for lid in _ids("ligands"):
        try:
            l = _read(_dir("ligands") / f"{lid}.json")
        except RegistryError as exc:
            errors.append(str(exc))
            continue
        ctx = f"ligands/{lid}.json"
        if not l.get("id"):
            errors.append(f"{ctx}: campo obrigatorio 'id' ausente")
            continue
        known_ligands.add(l["id"])
        if not l.get("smiles"):
            errors.append(f"{ctx}: campo obrigatorio 'smiles' ausente")

        topo = l.get("topology") or {}
        prep = l.get("prep") or {}
        if topo.get("is_macrocycle"):
            pdbqt = prep.get("pdbqt") or {}
            if pdbqt.get("rigid_macrocycles") is not False:
                errors.append(
                    f"{ctx}: macrociclo com prep.pdbqt.rigid_macrocycles != false -- "
                    "o anel ficaria rigido e a pose seria errada por construcao"
                )
            if not topo.get("largest_ring_size"):
                errors.append(f"{ctx}: macrociclo sem topology.largest_ring_size")
        prot = prep.get("protonation") or {}
        if prot and prot.get("required") and not prot.get("post_check"):
            errors.append(
                f"{ctx}: prep.protonation.required sem post_check -- a protonacao "
                "falha em silencio e o docking entrega pose errada sem aviso"
            )

    # --- protocolos
    for pid in _ids("protocols"):
        try:
            p = _read(_dir("protocols") / f"{pid}.json")
        except RegistryError as exc:
            errors.append(str(exc))
            continue
        ctx = f"protocols/{pid}.json"
        if not p.get("id"):
            errors.append(f"{ctx}: campo obrigatorio 'id' ausente")
            continue
        known_protocols.add(p["id"])
        if not p.get("engine"):
            errors.append(f"{ctx}: campo obrigatorio 'engine' ausente")

    # --- campanhas
    for cid in _ids("campaigns"):
        try:
            c = _read(_dir("campaigns") / f"{cid}.json")
        except RegistryError as exc:
            errors.append(str(exc))
            continue
        ctx = f"campaigns/{cid}.json"
        if not c.get("id"):
            errors.append(f"{ctx}: campo obrigatorio 'id' ausente")
        for tid in c.get("targets", []):
            if tid not in known_targets:
                errors.append(f"{ctx}: alvo '{tid}' nao existe em registry/targets/")
        for lid in c.get("ligands", []):
            if lid not in known_ligands:
                errors.append(f"{ctx}: ligante '{lid}' nao existe em registry/ligands/")
        for kind, pid in (c.get("protocols") or {}).items():
            if pid not in known_protocols:
                errors.append(f"{ctx}: protocolo '{pid}' (etapa {kind}) nao existe")
        default_ligand = c.get("default_ligand")
        if default_ligand and default_ligand not in c.get("ligands", []):
            errors.append(f"{ctx}: default_ligand '{default_ligand}' fora da lista de ligantes")

    # --- pares (onde vive a peculiaridade; validacao mais severa)
    for pid in _ids("pairs"):
        try:
            p = _read(_dir("pairs") / f"{pid}.json")
        except RegistryError as exc:
            errors.append(str(exc))
            continue
        ctx = f"pairs/{pid}.json"
        for field in ("campaign", "target", "ligand"):
            if not p.get(field):
                errors.append(f"{ctx}: campo obrigatorio '{field}' ausente")
        if p.get("target") and p["target"] not in known_targets:
            errors.append(f"{ctx}: alvo '{p['target']}' nao existe")
        if p.get("ligand") and p["ligand"] not in known_ligands:
            errors.append(f"{ctx}: ligante '{p['ligand']}' nao existe")

        asm = p.get("assembly_state") or {}
        required = asm.get("required")
        if required not in ("monomer", "dimer", "trimer", "filament", "complex"):
            errors.append(
                f"{ctx}: assembly_state.required invalido ({required!r}). "
                "O estado oligomerico e propriedade do PAR, nao do alvo."
            )
        if required == "filament" and not asm.get("min_protomers"):
            errors.append(f"{ctx}: assembly_state='filament' exige min_protomers")

        # Kd sem o estado em que foi medido e erro de categoria: citocalasina D
        # e 4,6 uM em G-actina e ~2 nM na extremidade barbada.
        for i, kd in enumerate(p.get("experimental_kd") or []):
            if not kd.get("state"):
                errors.append(
                    f"{ctx}: experimental_kd[{i}] sem 'state' -- comparar afinidade "
                    "sem o estado do alvo e erro de categoria"
                )
            for field in ("value", "unit"):
                if kd.get(field) in (None, ""):
                    errors.append(f"{ctx}: experimental_kd[{i}] sem '{field}'")

        val = p.get("validation") or {}
        if val.get("rmsd_function") == "GetBestRMS":
            errors.append(
                f"{ctx}: rmsd_function='GetBestRMS' alinha antes de medir e mascara "
                "pose errada (7,86 A pareceu 0,40 A). Use CalcRMS."
            )
        leak = val.get("training_leakage") or {}
        if leak and leak.get("held_out") is None:
            errors.append(f"{ctx}: training_leakage sem 'held_out' explicito")

        for kind, prot_id in (p.get("protocols") or {}).items():
            if prot_id not in known_protocols:
                errors.append(f"{ctx}: protocolo '{prot_id}' (etapa {kind}) nao existe")

        # Vocabulario de CONTROLE (matriz de cross-docking de especificidade). Um par
        # de controle negativo (droga sem sitio no alvo, ex.: colchicina x actina) TEM
        # de dizer isso EXPLICITAMENTE -- senao "afinidade fraca" fica indistinguivel de
        # "controle que funcionou". expected_binding coerente com o tipo (fail-loud).
        ctrl = p.get("control")
        if ctrl is not None:
            if not isinstance(ctrl, dict):
                errors.append(f"{ctx}: control deve ser objeto (got {type(ctrl).__name__})")
            else:
                ctype = ctrl.get("type")
                if ctype not in ("positive", "negative", None):
                    errors.append(
                        f"{ctx}: control.type invalido ({ctype!r}); esperado positive|negative|null")
                if ctype in ("positive", "negative"):
                    exp = ctrl.get("expected_binding")
                    if not isinstance(exp, bool):
                        errors.append(f"{ctx}: control.type={ctype!r} exige expected_binding booleano")
                    elif (ctype == "positive") != exp:
                        errors.append(
                            f"{ctx}: control incoerente -- type={ctype!r} mas expected_binding={exp} "
                            "(positivo=>true, negativo=>false)")
                    if not ctrl.get("rationale"):
                        errors.append(
                            f"{ctx}: control.type={ctype!r} exige rationale (por que liga/nao liga)")

    errors.extend(_validate_with_jsonschema())
    return errors


def _validate_with_jsonschema():
    """Checagem extra via _schema/*.schema.json, se jsonschema estiver instalado."""
    try:
        import jsonschema
    except ImportError:
        return []

    schema_dir = _dir("_schema")
    if not schema_dir.exists():
        return []

    kinds = {
        "target": "targets",
        "ligand": "ligands",
        "pair": "pairs",
        "campaign": "campaigns",
        "protocol": "protocols",
    }
    errors = []
    for kind, folder in kinds.items():
        schema_path = schema_dir / f"{kind}.schema.json"
        if not schema_path.exists():
            continue
        schema = _read(schema_path)
        for item_id in _ids(folder):
            data = _read(_dir(folder) / f"{item_id}.json")
            try:
                jsonschema.validate(data, schema)
            except jsonschema.ValidationError as exc:
                loc = "/".join(str(p) for p in exc.absolute_path) or "(raiz)"
                errors.append(f"{folder}/{item_id}.json [{loc}]: {exc.message}")
    return errors


# ---------------------------------------------------------------- proveniencia

# Ordem canonica dos tipos no snapshot. Fixa de proposito: o hash tem que ser
# estavel entre execucoes e maquinas, e qualquer reordenacao mudaria o sha256
# sem o conteudo ter mudado.
_SNAPSHOT_KINDS = ("campaigns", "targets", "ligands", "pairs", "protocols")


def snapshot():
    """Copia integral e determinista do registry, para embutir em proveniencia.

    O runner grava isto dentro de runs/<run_id>/run.json. Sem uma copia
    congelada, editar o registry amanha reescreve retroativamente o significado
    de um resultado de hoje -- e a generalizacao guarda-chuva destruiria a
    reprodutibilidade retroativa em vez de preserva-la. Le cru do disco (via
    _read), sem projecoes, para ser fiel ao que estava valendo no momento.
    """
    snap = {
        "registry_dir": str(REGISTRY_DIR),
        "program": _read(REGISTRY_DIR / "program.json"),
    }
    for kind in _SNAPSHOT_KINDS:
        snap[kind] = {
            item_id: _read(_dir(kind) / f"{item_id}.json")
            for item_id in _ids(kind)
        }
    return snap


def _canonical_json(obj):
    """JSON canonico (chaves ordenadas, sem espacos, UTF-8): a base do hash."""
    return json.dumps(obj, sort_keys=True, ensure_ascii=False,
                      separators=(",", ":"))


def snapshot_sha256(snap=None):
    """sha256 hex do snapshot canonico. Identifica a versao do registry usada.

    Estavel sob reformatacao (indentacao, ordem de chaves) porque hasheia a
    forma canonica, nao o texto bruto dos arquivos: dois registries com o mesmo
    conteudo logico produzem o mesmo hash.
    """
    if snap is None:
        snap = snapshot()
    return hashlib.sha256(_canonical_json(snap).encode("utf-8")).hexdigest()


def main():
    import argparse

    ap = argparse.ArgumentParser(description="Inspeciona e valida o registry.")
    ap.add_argument("--validate", action="store_true", help="valida tudo e sai")
    ap.add_argument("--campaign", help="lista os alvos de uma campanha")
    ap.add_argument("--sha256", action="store_true",
                    help="imprime o sha256 do snapshot do registry e sai")
    args = ap.parse_args()

    if args.sha256:
        print(snapshot_sha256())
        return 0

    print(f"registry: {REGISTRY_DIR}")

    if args.campaign:
        for t in list_targets(args.campaign):
            site = resolve_site(t["id"])
            print(f"  {t['id']:12s} {str((t.get('structure') or {}).get('pdb_id')):8s} "
                  f"{site['definition']:32s} {t.get('organism')}")
        return 0

    print(f"  campanhas : {list_campaigns()}")
    print(f"  alvos     : {len(_ids('targets'))}")
    print(f"  ligantes  : {len(_ids('ligands'))}")
    print(f"  pares     : {len(_ids('pairs'))}")
    print(f"  protocolos: {len(_ids('protocols'))}")

    errors = validate_all()
    if errors:
        print(f"\n{len(errors)} erro(s) de validacao:")
        for e in errors:
            print(f"  - {e}")
        return 1
    print("\n[ok] registry valido")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
