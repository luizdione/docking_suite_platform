# -*- coding: utf-8 -*-
"""
haddock_analysis.py -- parsing e scoring da SAIDA do HADDOCK3 (clustfcc.tsv,
capri_clt.tsv), sem depender do binario haddock3/CNS.

Motivacao: HADDOCK score bruto (CNS, unidade arbitraria) NAO e comparavel
entre sistemas de tamanhos diferentes (ponto 1 do revisor). Este modulo da o
Z-score de cada cluster contra a propria populacao de decoys do sistema
(`zscore_by_cluster`) e um gate pass/fail documentado (`decide`) que usa Z e,
quando disponivel, a banda DockQ da comunidade CAPRI -- em vez de comparar
score CNS cru entre projetos/alvos diferentes.

Consumido por fla_fla e tcactina (nao duplicar o parser la -- ver
`pip install -e C:\\github\\docking_suite_platform` / pyproject.toml deste
repo). Este modulo so LE arquivos .tsv que o haddock3 ja escreveu -- nao roda
o binario, nao depende de CNS/WSL.

Contrato de arquivos (confirmado contra dado real, ver tests/ Tarefa 4):
    - clustfcc.tsv (escrito pelo modulo `clustfcc`, ex.:
      `run.../5_clustfcc/clustfcc.tsv`): score CNS INDIVIDUAL de cada modelo
      avaliado -- rank, model_name, score, cluster_id ("-" se nao
      clusterizado). E a "populacao de decoys" usada no Z-score.
    - capri_clt.tsv (escrito pelo `caprieval` em modo cluster-based, ex.:
      `run.../analysis/<N>_<modulo>_analysis/capri_clt.tsv`): UMA linha por
      cluster, com score medio (dos top-N do cluster, N = clt_threshold do
      cfg) e metricas de qualidade estrutural (irmsd/lrmsd/dockq/...).
"""
from __future__ import annotations

import csv
import math
from pathlib import Path

# Colunas de capri_clt.tsv que sao inteiras por natureza (rank/id/contagem),
# nao metrica continua.
_CAPRI_INT_COLS = {"cluster_rank", "cluster_id", "n", "caprieval_rank"}
# Colunas de texto cru (nao numerico) -- "yes"/"no"/"-" no formato real.
_CAPRI_STR_COLS = {"under_eval"}


def _clean(raw) -> str | None:
    """'-' ou vazio -> None (convencao do HADDOCK3 pra celula ausente/NA)."""
    val = (raw or "").strip()
    return None if val in ("", "-") else val


def parse_clustfcc_tsv(path) -> list[dict]:
    """Le clustfcc.tsv: score CNS individual de CADA modelo avaliado (a
    populacao de decoys, pre-agregacao por cluster).

    Colunas reais do arquivo (header na 1a linha, TAB-separado):
        rank  model_name  score  cluster_id

    `cluster_id` e a string "-" quando o modelo nao caiu em nenhum cluster
    (fica None no dict de saida) -- esses modelos SAO incluidos na populacao
    usada por `zscore_by_cluster` (a populacao de decoys inclui os
    nao-clusterizados, so nao entram na media de nenhum cluster especifico).

    Retorna uma lista com uma linha por MODELO (nao por cluster), na ordem do
    arquivo:
        [{"rank": int, "model_name": str, "score": float,
          "cluster_id": int | None}, ...]
    """
    path = Path(path)
    rows: list[dict] = []
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for raw in reader:
            model_name = _clean(raw.get("model_name"))
            if model_name is None:
                continue  # linha em branco final (comum no arquivo real)
            cid = _clean(raw.get("cluster_id"))
            rows.append({
                "rank": int(raw["rank"]),
                "model_name": model_name,
                "score": float(raw["score"]),
                "cluster_id": int(cid) if cid is not None else None,
            })
    return rows


def parse_capri_clt(path) -> list[dict]:
    """Le capri_clt.tsv (`analysis/<N>_<modulo>_analysis/capri_clt.tsv`) -- UM
    dict por CLUSTER, com as metricas do modulo `caprieval` (cluster-based).

    O arquivo real comeca com um bloco de comentarios (linhas iniciando com
    '#', com o `sortby_key`/`clt_threshold` usados) antes da linha de colunas
    de verdade; comentarios e linhas em branco sao ignorados automaticamente.

    Colunas (schema real confirmado contra dado do fla_fla, ver Tarefa 4):
        cluster_rank  cluster_id  n  under_eval  score  score_std  irmsd
        irmsd_std  fnat  fnat_std  lrmsd  lrmsd_std  dockq  dockq_std  ilrmsd
        ilrmsd_std  rmsd  rmsd_std  air  air_std  bsa  bsa_std  desolv
        desolv_std  elec  elec_std  total  total_std  vdw  vdw_std
        caprieval_rank

    IMPORTANTE -- auto-referencia vs. referencia real: se o .cfg de origem NAO
    tinha um bloco [caprieval] final com `reference_fname` explicito (ver
    `haddock3.build_cfg(caprieval_reference=...)`), os campos
    irmsd/lrmsd/dockq/rmsd aqui sao AUTO-REFERENCIADOS -- CNS compara cada
    modelo contra o proprio melhor modelo do HADDOCK, por isso o cluster
    caprieval_rank=1 sempre sai com irmsd=0.0/dockq=1.0 por construcao. Isso
    NAO e uma medida de qualidade contra a realidade (ex.: o modelo do Boltz)
    -- so vira isso quando `reference_fname` aponta pra uma estrutura de
    referencia de verdade. `decide()` documenta a mesma ressalva.

    Retorna uma lista de dicts (um por cluster, na ordem do arquivo) com todos
    os campos numericos convertidos: int para cluster_rank/cluster_id/n/
    caprieval_rank, float para as demais colunas numericas, None onde o
    arquivo tem "-", e string crua para `under_eval`.
    """
    path = Path(path)
    with path.open(encoding="utf-8") as fh:
        data_lines = [ln for ln in fh
                      if ln.strip() and not ln.lstrip().startswith("#")]
    if not data_lines:
        raise ValueError(f"{path}: nenhuma linha de dados (so comentarios/vazio)")

    reader = csv.DictReader(data_lines, delimiter="\t")
    clusters: list[dict] = []
    for raw in reader:
        row: dict = {}
        for key, val in raw.items():
            key = key.strip()
            cleaned = _clean(val)
            if key in _CAPRI_STR_COLS:
                row[key] = cleaned
            elif cleaned is None:
                row[key] = None
            elif key in _CAPRI_INT_COLS:
                row[key] = int(float(cleaned))
            else:
                row[key] = float(cleaned)
        clusters.append(row)
    return clusters


def zscore_by_cluster(clustfcc_tsv_path) -> dict:
    """Z-score de CADA cluster vs. a populacao COMPLETA de modelos avaliados
    (le SO clustfcc.tsv -- nao junta com capri_clt.tsv).

    Populacao = campo `score` de TODAS as linhas de clustfcc.tsv, clusterizados
    ou nao -- e a populacao de decoys real que o HADDOCK avaliou naquele
    sistema, nao so os que sobreviveram ao clustering. Media e desvio-padrao
    POPULACIONAIS (ddof=0, divide por N): os N modelos avaliados SAO a
    populacao de interesse aqui (nao uma amostra de uma populacao maior).

        Z(cluster) = (media_score_do_cluster - media_populacao) / desvio_populacao

    `media_score_do_cluster` usa TODOS os modelos daquele cluster_id em
    clustfcc.tsv (nao so o top-N que aparece em capri_clt.tsv -- os dois
    arquivos descrevem populacoes diferentes de proposito: clustfcc.tsv = todo
    modelo avaliado, capri_clt.tsv = media dos top-N/cluster (clt_threshold do
    cfg) usados no relatorio final).

    Mais negativo = cluster melhor (mais favoravel) que a media dos decoys
    daquele MESMO sistema -- e o que substitui comparar score CNS bruto entre
    sistemas de tamanhos diferentes (pedido explicito do revisor).

    Retorna {cluster_id (int): z (float)} para os clusters presentes (exclui
    "-"/nao-clusterizado, que nao tem cluster_id). Se o desvio-padrao da
    populacao for 0 (todos os scores identicos -- caso degenerado, nao
    esperado com dados reais), devolve 0.0 para todos os clusters em vez de
    levantar ZeroDivisionError.
    """
    rows = parse_clustfcc_tsv(clustfcc_tsv_path)
    scores = [r["score"] for r in rows]
    n = len(scores)
    if n == 0:
        raise ValueError(f"{clustfcc_tsv_path}: nenhum modelo em clustfcc.tsv")

    mean = sum(scores) / n
    variance = sum((s - mean) ** 2 for s in scores) / n
    std = math.sqrt(variance)

    by_cluster: dict[int, list[float]] = {}
    for r in rows:
        cid = r["cluster_id"]
        if cid is None:
            continue
        by_cluster.setdefault(cid, []).append(r["score"])

    if std == 0:
        return {cid: 0.0 for cid in by_cluster}
    return {
        cid: (sum(vals) / len(vals) - mean) / std
        for cid, vals in by_cluster.items()
    }


def decide(cluster: dict, z_max: float = -1.0, dockq_min: float = 0.23) -> dict:
    """Gate pass/fail de UM cluster HADDOCK -- substitui comparar score CNS
    bruto entre sistemas de tamanhos diferentes (ponto 1 do revisor).

    `cluster` e um dict simples que PRECISA ter a chave "z" ja calculada --
    esta funcao nao le arquivo nenhum. Padrao de uso: junte a mao o resultado
    de `zscore_by_cluster(clustfcc_tsv)` (dict cluster_id -> z) com cada dict
    de `parse_capri_clt(capri_clt_tsv)` pelo campo `cluster_id`, por exemplo:

        zs = zscore_by_cluster(clustfcc_tsv)
        for c in parse_capri_clt(capri_clt_tsv):
            c["z"] = zs.get(c["cluster_id"])
            print(decide(c))

    Chave "dockq" e OPCIONAL dentro de `cluster`:
      - None/ausente: DockQ indisponivel (capri_clt.tsv veio de um cfg SEM
        [caprieval] reference_fname explicito) -- decide-se so pelo Z-score.
      - presente: SO deve vir de um capri_clt.tsv gerado COM reference_fname
        apontando pra uma estrutura de referencia REAL (ex.: o modelo do
        Boltz). NUNCA passe aqui o dockq auto-referenciado (cfg sem
        reference_fname) -- ele e trivialmente ~1.0 pro cluster
        caprieval_rank=1 por construcao (comparado contra si mesmo, ver
        docstring de parse_capri_clt), nao mede qualidade nenhuma.

    Banda DockQ CAPRI (Basu & Wallner 2016, PLoS ONE -- convencao da
    comunidade CAPRI): dockq >= 0.23 aceitavel, >= 0.49 medio, >= 0.80 alto.
    `dockq_min` default (0.23) usa o piso "aceitavel". Ver
    `haddock_config.HADDOCK_THRESH` para os defaults documentados do projeto.

    Regra de decisao: passa por Z-score (z <= z_max) E, SE dockq estiver
    disponivel, tambem passa por DockQ (dockq >= dockq_min) -- os dois
    criterios tem que concordar quando ambos existem; nenhum dos dois sozinho
    basta para reprovar o outro em silencio (o "reason" sempre explicita qual
    criterio decidiu).

    Retorna {"passes": bool, "reason": str, "z": float | None,
             "dockq": float | None}.
    """
    z = cluster.get("z")
    dockq = cluster.get("dockq")

    if z is None:
        return {
            "passes": False,
            "reason": ("z ausente no dict do cluster -- junte o resultado de "
                       "zscore_by_cluster() ao dict antes de chamar decide()"),
            "z": None,
            "dockq": dockq,
        }

    z_ok = z <= z_max
    if dockq is None:
        passes = z_ok
        reason = (f"z={z:.3f} {'<=' if z_ok else '>'} z_max={z_max} vs. "
                  f"populacao de decoys (dockq indisponivel -- capri_clt.tsv "
                  f"sem reference_fname vs. referencia real no cfg de origem)")
    else:
        dockq_ok = dockq >= dockq_min
        passes = z_ok and dockq_ok
        if z_ok and dockq_ok:
            reason = (f"passa: z={z:.3f} <= z_max={z_max} e dockq={dockq:.3f} "
                      f">= dockq_min={dockq_min} (banda CAPRI aceitavel)")
        elif not z_ok and not dockq_ok:
            reason = (f"reprova em ambos: z={z:.3f} > z_max={z_max} e "
                      f"dockq={dockq:.3f} < dockq_min={dockq_min}")
        elif not z_ok:
            reason = (f"reprova por z: z={z:.3f} > z_max={z_max} (apesar de "
                      f"dockq={dockq:.3f} >= dockq_min={dockq_min})")
        else:
            reason = (f"reprova por dockq: dockq={dockq:.3f} < "
                      f"dockq_min={dockq_min} (abaixo da banda CAPRI "
                      f"aceitavel, apesar de z={z:.3f} <= z_max={z_max})")

    return {"passes": passes, "reason": reason, "z": z, "dockq": dockq}
