# -*- coding: utf-8 -*-
"""
Gate de haddock_analysis.py / haddock_config.py.

Cobre:
  - parse_clustfcc_tsv / parse_capri_clt contra fixtures SINTETICAS (formato
    exato confirmado contra dado real, ver abaixo);
  - zscore_by_cluster (populacao completa vs. media/cluster, ddof=0);
  - decide() nos 4 quadrantes (z ok/nao-ok x dockq ausente/ok/nao-ok);
  - VALIDACAO contra dado REAL de outro projeto (fla_fla, run1-ecto):
    confirma que os 3 clusters batem com os numeros ja publicados no
    relatorio (score/score_std do capri_clt.tsv) -- pulado (nao falha) se a
    pasta de dados nao existir nesta maquina (dado vive fora deste repo, em
    outro projeto do usuario).

Uso:
    pytest tests/test_haddock_analysis.py -v
    python tests/test_haddock_analysis.py
"""
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import pytest  # noqa: E402

from haddock_analysis import (  # noqa: E402
    decide,
    parse_capri_clt,
    parse_clustfcc_tsv,
    zscore_by_cluster,
)
from haddock_config import HADDOCK_THRESH  # noqa: E402

# ---------------------------------------------------------------- fixtures

CLUSTFCC_TSV = """rank\tmodel_name\tscore\tcluster_id
1\temref_6.pdb\t-169.76\t-
2\temref_16.pdb\t-160.16\t-
3\temref_2.pdb\t-137.35\t-
4\temref_19.pdb\t-134.87\t3
5\temref_12.pdb\t-129.95\t2
6\temref_17.pdb\t-129.45\t-
7\temref_8.pdb\t-129.12\t1
8\temref_27.pdb\t-122.49\t-
9\temref_30.pdb\t-119.33\t3
10\temref_34.pdb\t-115.22\t1
"""

CAPRI_CLT_TSV = """########################################
# `caprieval` cluster-based analysis
#
# > sortby_key=score
# > sort_ascending=True
# > clt_threshold=4
#
########################################
cluster_rank\tcluster_id\tn\tunder_eval\tscore\tscore_std\tirmsd\tirmsd_std\tfnat\tfnat_std\tlrmsd\tlrmsd_std\tdockq\tdockq_std\tilrmsd\tilrmsd_std\trmsd\trmsd_std\tair\tair_std\tbsa\tbsa_std\tdesolv\tdesolv_std\telec\telec_std\ttotal\ttotal_std\tvdw\tvdw_std\tcaprieval_rank
1\t1\t7\t-\t-114.575\t9.055\t24.186\t0.168\t0.055\t0.005\t87.845\t0.559\t0.023\t0.002\t54.043\t0.401\t38.765\t0.704\t57.903\t22.577\t2384.900\t83.688\t-7.915\t2.229\t-197.164\t39.490\t-212.278\t60.375\t-73.018\t4.138\t1
2\t3\t4\t-\t-112.688\t15.538\t1.269\t0.746\t0.672\t0.196\t4.115\t2.634\t0.678\t0.195\t2.932\t1.796\t1.646\t1.001\t46.174\t33.121\t2471.185\t137.991\t-8.168\t1.803\t-176.686\t20.985\t-204.314\t48.123\t-73.801\t10.952\t2
3\t2\t5\t-\t-105.703\t16.522\t20.697\t0.160\t0.020\t0.005\t79.267\t0.503\t0.012\t0.002\t45.150\t0.374\t30.473\t0.177\t52.173\t12.708\t2105.762\t135.927\t-5.551\t3.501\t-270.537\t49.903\t-269.627\t63.074\t-51.262\t2.964\t3
"""


@pytest.fixture
def clustfcc_path(tmp_path):
    p = tmp_path / "clustfcc.tsv"
    p.write_text(CLUSTFCC_TSV, encoding="utf-8")
    return p


@pytest.fixture
def capri_clt_path(tmp_path):
    p = tmp_path / "capri_clt.tsv"
    p.write_text(CAPRI_CLT_TSV, encoding="utf-8")
    return p


# ---------------------------------------------------------- parse_clustfcc_tsv

def test_parse_clustfcc_tsv_basic(clustfcc_path):
    rows = parse_clustfcc_tsv(clustfcc_path)
    assert len(rows) == 10
    first = rows[0]
    assert first == {"rank": 1, "model_name": "emref_6.pdb", "score": -169.76,
                      "cluster_id": None}
    clustered = [r for r in rows if r["cluster_id"] is not None]
    assert len(clustered) == 5
    by_name = {r["model_name"]: r for r in rows}
    assert by_name["emref_8.pdb"]["cluster_id"] == 1
    assert by_name["emref_8.pdb"]["score"] == -129.12


def test_parse_clustfcc_tsv_trailing_blank_line_ignored(tmp_path):
    p = tmp_path / "clustfcc.tsv"
    p.write_text(CLUSTFCC_TSV + "\n\n", encoding="utf-8")
    rows = parse_clustfcc_tsv(p)
    assert len(rows) == 10


# ------------------------------------------------------------ parse_capri_clt

def test_parse_capri_clt_skips_comment_header(capri_clt_path):
    clusters = parse_capri_clt(capri_clt_path)
    assert len(clusters) == 3


def test_parse_capri_clt_field_types_and_values(capri_clt_path):
    clusters = parse_capri_clt(capri_clt_path)
    c1 = clusters[0]
    assert c1["cluster_rank"] == 1
    assert c1["cluster_id"] == 1
    assert c1["n"] == 7
    assert c1["under_eval"] is None  # "-" no arquivo real
    assert c1["score"] == pytest.approx(-114.575)
    assert c1["score_std"] == pytest.approx(9.055)
    assert c1["dockq"] == pytest.approx(0.023)
    assert c1["caprieval_rank"] == 1
    assert isinstance(c1["cluster_rank"], int)
    assert isinstance(c1["score"], float)


# ------------------------------------------------------------ zscore_by_cluster

def test_zscore_by_cluster_matches_manual_population_stats(clustfcc_path):
    rows = parse_clustfcc_tsv(clustfcc_path)
    scores = [r["score"] for r in rows]
    mean = sum(scores) / len(scores)
    pstd = math.sqrt(sum((s - mean) ** 2 for s in scores) / len(scores))

    zs = zscore_by_cluster(clustfcc_path)
    assert set(zs) == {1, 2, 3}

    cluster1_scores = [r["score"] for r in rows if r["cluster_id"] == 1]
    expected_z1 = (sum(cluster1_scores) / len(cluster1_scores) - mean) / pstd
    assert zs[1] == pytest.approx(expected_z1)


def test_zscore_by_cluster_ranks_clusters_by_mean_score(clustfcc_path):
    # cluster 2 (so -129.95) < cluster 3 (media -127.10) < cluster 1 (media
    # -122.17) em score -- mais negativo = melhor, entao o Z-score de cluster 2
    # deve ser o mais negativo dos tres, e o de cluster 1 o menos negativo.
    zs = zscore_by_cluster(clustfcc_path)
    assert zs[2] < zs[3] < zs[1]


def test_zscore_by_cluster_zero_std_no_crash(tmp_path):
    p = tmp_path / "degenerate.tsv"
    p.write_text(
        "rank\tmodel_name\tscore\tcluster_id\n"
        "1\ta.pdb\t-100.0\t1\n"
        "2\tb.pdb\t-100.0\t1\n"
        "3\tc.pdb\t-100.0\t2\n",
        encoding="utf-8",
    )
    zs = zscore_by_cluster(p)
    assert zs == {1: 0.0, 2: 0.0}


def test_zscore_by_cluster_empty_raises(tmp_path):
    p = tmp_path / "empty.tsv"
    p.write_text("rank\tmodel_name\tscore\tcluster_id\n", encoding="utf-8")
    with pytest.raises(ValueError):
        zscore_by_cluster(p)


# ------------------------------------------------------------------- decide()

def test_decide_missing_z_fails_closed():
    result = decide({"dockq": 0.9})
    assert result["passes"] is False
    assert result["z"] is None
    assert "z ausente" in result["reason"]


def test_decide_z_ok_no_dockq_passes():
    result = decide({"z": -1.5})
    assert result["passes"] is True
    assert result["dockq"] is None


def test_decide_z_not_ok_no_dockq_fails():
    result = decide({"z": -0.2})
    assert result["passes"] is False


def test_decide_z_ok_dockq_ok_passes():
    result = decide({"z": -2.0, "dockq": 0.5})
    assert result["passes"] is True


def test_decide_z_ok_dockq_not_ok_fails():
    result = decide({"z": -2.0, "dockq": 0.1})
    assert result["passes"] is False
    assert "dockq" in result["reason"]


def test_decide_z_not_ok_dockq_ok_fails():
    result = decide({"z": 0.5, "dockq": 0.9})
    assert result["passes"] is False


def test_decide_uses_thresh_defaults_matching_haddock_config():
    # decide() default kwargs devem bater com HADDOCK_THRESH (mesma fonte de
    # verdade, so nao importado automaticamente para deixar a chamada explicita
    # a criterio de quem consome).
    result_default = decide({"z": -1.0, "dockq": 0.23})
    result_explicit = decide({"z": -1.0, "dockq": 0.23},
                              z_max=HADDOCK_THRESH["z_max"],
                              dockq_min=HADDOCK_THRESH["dockq_min"])
    assert result_default["passes"] == result_explicit["passes"]


def test_decide_custom_thresholds():
    result = decide({"z": -0.5}, z_max=-0.3)
    assert result["passes"] is True


# ------------------------------------------------------------------ end-to-end

def test_end_to_end_join_capri_clt_with_zscore(clustfcc_path, capri_clt_path):
    zs = zscore_by_cluster(clustfcc_path)
    clusters = parse_capri_clt(capri_clt_path)
    for c in clusters:
        c["z"] = zs.get(c["cluster_id"])
    decisions = {c["cluster_id"]: decide(c) for c in clusters}
    # capri_clt_path so tem dockq auto-referenciado (sem reference_fname no
    # cfg de origem) -- este teste NAO passa dockq pra decide() de proposito
    # (simula o uso correto: dockq so entra quando vem de referencia real).
    for c in clusters:
        assert decisions[c["cluster_id"]]["z"] == pytest.approx(zs[c["cluster_id"]])


# ------------------------------------------------------- validacao dado real

FLA_FLA_HADDOCK_DIR = Path(
    r"C:\Users\luizd\OneDrive\Data_Science\Projetos\colaboracoes\fla_fla"
    r"\results\haddock_ecto\run1-ecto"
)
REAL_CLUSTFCC = FLA_FLA_HADDOCK_DIR / "5_clustfcc" / "clustfcc.tsv"
REAL_CAPRI_CLT = (FLA_FLA_HADDOCK_DIR / "analysis" / "5_clustfcc_analysis"
                  / "capri_clt.tsv")

_real_data_available = REAL_CLUSTFCC.exists() and REAL_CAPRI_CLT.exists()


@pytest.mark.skipif(not _real_data_available,
                    reason="dado real do fla_fla (run1-ecto) nao encontrado "
                           "nesta maquina -- e SOMENTE LEITURA de outro "
                           "projeto do usuario, nao versionado neste repo")
def test_real_fla_fla_run1_ecto_matches_published_report():
    """Tarefa 4: valida parse_capri_clt/zscore_by_cluster contra o dado REAL
    ja em disco do fla_fla (run1-ecto). Confirma que os 3 clusters batem com
    os numeros ja publicados no relatorio de fla_fla (PROJECT_STATE.md /
    haddock_summary.json): cluster 1 (n=7) -114.575+-9.055 [cluster_id=1],
    cluster 2 (n=4) -112.688+-15.538 [cluster_id=3, rank 2], cluster 3 (n=5)
    -105.703+-16.522 [cluster_id=2, rank 3]. So LEITURA -- nada e escrito na
    pasta do fla_fla.
    """
    clusters = parse_capri_clt(REAL_CAPRI_CLT)
    assert len(clusters) == 3
    by_rank = {c["cluster_rank"]: c for c in clusters}

    c1 = by_rank[1]
    assert c1["cluster_id"] == 1
    assert c1["n"] == 7
    assert c1["score"] == pytest.approx(-114.575, abs=0.01)
    assert c1["score_std"] == pytest.approx(9.055, abs=0.01)

    c2 = by_rank[2]
    assert c2["cluster_id"] == 3
    assert c2["n"] == 4
    assert c2["score"] == pytest.approx(-112.688, abs=0.01)
    assert c2["score_std"] == pytest.approx(15.538, abs=0.01)

    c3 = by_rank[3]
    assert c3["cluster_id"] == 2
    assert c3["n"] == 5
    assert c3["score"] == pytest.approx(-105.703, abs=0.01)
    assert c3["score_std"] == pytest.approx(16.522, abs=0.01)

    # zscore_by_cluster le SO clustfcc.tsv (populacao completa de modelos
    # avaliados, clusterizados ou nao) -- valores conferidos manualmente
    # contra os 40 modelos reais do run1-ecto (ver PROJECT_STATE.md /
    # validacao desta tarefa).
    zs = zscore_by_cluster(REAL_CLUSTFCC)
    assert set(zs) == {1, 2, 3}
    assert zs[1] == pytest.approx(-0.2938, abs=0.001)
    assert zs[2] == pytest.approx(-0.1245, abs=0.001)
    assert zs[3] == pytest.approx(-0.5178, abs=0.001)

    # decide() com os defaults do projeto (z_max=-1.0): nenhum dos 3 clusters
    # passa no corte provisorio -- consistente com a leitura "WEAK" ja no
    # relatorio do fla_fla (nenhum cluster e >=1 desvio-padrao melhor que a
    # media dos proprios decoys).
    for cid, z in zs.items():
        result = decide({"z": z})
        assert result["passes"] is False, (cid, z, result)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
