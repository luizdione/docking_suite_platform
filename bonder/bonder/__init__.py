#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""bonder - plataforma de analise de CoIP-MS / AP-MS para identificacao de
parceiros de interacao no proteoma de camundongo (Mus musculus, taxid 10090).

Filosofia: NAO reinventa a roda. A plataforma ORQUESTRA ferramentas consagradas
e revisadas por pares da comunidade de proteomica (SAINTexpress, MiST, MSstats,
artMS, limma/DEP, DIA-NN, MS2Rescore, STRING/BioGRID) dentro de conteineres
reprodutiveis, expondo tudo num dashboard Flask acessivel via navegador local.

Camadas (fluxo de dados unidirecional io -> db -> scoring/stats -> gpu -> pipeline -> web):
  bonder.io       parsers dos motores de busca -> um MSDataset unificado
  bonder.db       bancos murinos (UniProt/STRING/BioGRID/cRAP/CRAPome) + iscas
  bonder.scoring  formatacao e execucao de SAINTexpress / MiST (especificidade)
  bonder.stats    wrappers R (MSstats/artMS/limma/DEP) para o quantitativo
  bonder.annotate anotacao (gene/GO/via) e montagem da rede de interacao
  bonder.gpu      DIA-NN / MS2Rescore / AlphaFold / ESM2 com fallback CPU
  bonder.pipeline orquestrador + relatorio + manifesto de reprodutibilidade
  bonder.web      dashboard Flask (localhost) + fila de jobs

O nucleo cientifico (io/scoring-formatacao/annotate) e testavel em CPU, sem GPU
nem R, para rodar em CI e em qualquer maquina.
"""

__version__ = "0.1.0"
__all__ = ["__version__"]
