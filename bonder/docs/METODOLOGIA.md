# Metodologia cientifica e citacoes

Este documento descreve o metodo do `bonder` e as referencias das ferramentas
orquestradas. Todo escore que sustenta uma conclusao vem de uma ferramenta
**publicada e revisada por pares** — nada e reimplementado.

## 1. Desenho

- **Amostra:** *Mus musculus* (taxid 10090), celulas in vitro.
- **Iscas (exogenas):** NSs do ortobunyavirus *Icoaraci* (viral) e Galectina-3
  humana (LGALS3, P17931). Parceiros = proteinas murinas.
- **Controles:** IgG/mock/linhagem parental, pareados por isca.
- **Replicas:** recomenda-se **>= 3 biologicas** por grupo. Com n=2-4, a decisao
  prioriza especificidade (SAINT/MiST) e presenca/ausencia vs controle, com o
  quantitativo (DEP/limma) como suporte.

## 2. Ingesta e QC

- Parsers por motor (MaxQuant/DIA-NN/FragPipe/generico) -> matriz unica
  (`MSDataset`), **chaveada pelo acesso UniProt lider**, com **ausente = NaN**
  (0 do MaxQuant/FragPipe convertido para NaN em DDA; MaxLFQ do DIA-NN **nao** e
  re-normalizado).
- Contaminantes/decoys sinalizados por coluna booleana **e** por prefixo
  (`CON__`/`REV__`/`contam_`); iscas nunca sao tratadas como contaminante.
- QC: valores validos por amostra/condicao, correlacao entre replicas, log2.

## 3. Especificidade isca-presa (nucleo do CoIP)

- **SAINTexpress** (Teo et al., *J Proteomics* 2014; Choi et al., *Nat Methods*
  2011) — probabilidade SAINT + **BFDR** por prey vs controles. Padrao de fato
  em AP-MS. Iscas exogenas sao excluidas do conjunto de preys murinos.
- **MiST** (Verschueren et al., 2015; Jager et al., *Nature* 2011) — combina
  reprodutibilidade, abundancia e especificidade; recomendado para a **isca
  viral NSs** (paradigma host-pathogen, ex.: Gordon et al., *Nature* 2020).

## 4. Quantitativo diferencial

- **DEP** (Zhang et al., *Nat Protoc* 2018) — enriquecimento diferencial em
  nivel de proteina: filtragem, normalizacao **vsn**, imputacao MNAR
  (**MinProb**) e teste via **limma** (Ritchie et al., *NAR* 2015; teste-t
  moderado por Bayes empirico). Ideal para matrizes de proteinGroups.
- **MSstats** (Choi et al., *Bioinformatics* 2014) — quando ha dados em nivel de
  feature (evidence/report), usando os conversores oficiais.

## 5. Background / contaminantes

- **CRAPome** (Mellacheruvu et al., *Nat Methods* 2013) — frequencia de background;
  aqui priorizamos a frequencia observada nos **proprios controles** do usuario.

## 6. Anotacao e rede

- **UniProt** (proteoma murino), **STRING v12.0** (Szklarczyk et al.), **BioGRID**
  (Oughtred et al.), GO via **g:Profiler** (Kolberg et al.). Rede exportada em
  GraphML/Cytoscape (Shannon et al., *Genome Res* 2003).

## 7. GPU (opcional, RTX 3060)

- **DIA-NN** (Demichev et al., *Nat Methods* 2020) — busca DIA acelerada.
- **MS2Rescore** (Declercq et al.) + **AlphaPeptDeep** (Zeng et al., *Nat
  Commun* 2022) — rescoring por deep learning; +IDs a FDR fixo.
- **AlphaFold2-Multimer** (Evans et al.) / **ColabFold** (Mirdita et al., *Nat
  Methods* 2022) — validacao estrutural isca-presa dos melhores candidatos.
- **ESM2** (Lin et al., *Science* 2023) — priorizacao **exploratoria/suplementar**;
  nunca a base do escore de confianca.

## 8. Ranqueamento final

Interatores ordenados por **BFDR (SAINT)**, depois **FDR (DEP/limma)**; os
descritivos (fold change de medias, contagem de observados) entram apenas como
QC/desempate. Rotulo de confianca: alta (BFDR <= 0.05), media (<= 0.10), baixa.

## 9. Reporte para publicacao

Depositar dados brutos no **ProteomeXchange/PRIDE**; reportar versoes das
ferramentas (do `manifest.json`), parametros de busca, criterios de FDR e o
desenho de replicas/controles. Ver `docs/REPRODUTIBILIDADE.md`.
