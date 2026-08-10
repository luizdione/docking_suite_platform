# -*- coding: utf-8 -*-
"""
haddock_config.py -- limiares de decisao para triagem de resultados HADDOCK3.

Constantes PROVISORIAS/AJUSTAVEIS consumidas por `haddock_analysis.decide()`.
Mesmo espirito de `ppi/config.py:SCORE` no projeto irmao ppi_platform: nao sao
verdade absoluta, sao um ponto de partida documentado. Qualquer projeto que
consome este pacote (fla_fla, tcactina, ...) pode sobrescrever chamando
`decide(cluster, z_max=..., dockq_min=...)` com seus proprios valores em vez
destes defaults.
"""
from __future__ import annotations

HADDOCK_THRESH = {
    # Z-score do cluster vs. a populacao de decoys avaliados no MESMO sistema
    # (ver haddock_analysis.zscore_by_cluster). Mais negativo = complexo
    # melhor (mais favoravel, HADDOCK score) que a media dos decoys daquele
    # sistema. -1.0 (1 desvio-padrao populacional abaixo da media) e um corte
    # PROVISORIO -- ao contrario do DockQ (banda CAPRI publicada), nao ha uma
    # convencao da comunidade equivalente para "quanto de Z-score e bom o
    # bastante" em docking guiado por AIR; ajuste por sistema/projeto.
    "z_max": -1.0,
    # DockQ (Basu & Wallner 2016, "DockQ: A Quality Measure for Protein-Protein
    # Docking Models", PLoS ONE) so e um sinal valido quando o cluster vem de
    # um capri_clt.tsv gerado com [caprieval] reference_fname apontando pra
    # uma estrutura de referencia REAL (ex.: o modelo do Boltz) -- NUNCA o
    # dockq auto-referenciado do proprio HADDOCK (ver docstring de
    # haddock_analysis.parse_capri_clt / decide). Banda CAPRI (convencao da
    # comunidade): dockq >= 0.23 aceitavel | >= 0.49 medio | >= 0.80 alto.
    # dockq_min=0.23 usa o piso "aceitavel".
    "dockq_min": 0.23,
}
