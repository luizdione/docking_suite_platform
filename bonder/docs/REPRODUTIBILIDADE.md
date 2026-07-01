# Reprodutibilidade

Principios (mesma filosofia da suite docking: nada silencioso, tudo rastreavel).

## 1. Bancos fixados por versao + checksum

`config/databases.yaml` fixa versao, URL e sha256 de cada fonte (UniProt release,
STRING v12.0, BioGRID build, cRAP). O downloader:

- grava no cache content-addressed (`data/`, ou `BONDER_DATA_ROOT`);
- mantem `data/registry.json` com url, sha256, tamanho e timestamp;
- **falha** se o checksum divergir (nao re-baixa silenciosamente);
- suporta **modo offline** (`BONDER_OFFLINE=1`): usa so o cache.

Na primeira `bonder db sync`, preencha os `sha256_*` em `databases.yaml` com os
hashes observados e faca commit — travando as versoes para sempre.

## 2. Manifesto por corrida

Cada `runs/<run_id>/manifest.json` registra:

- versao do `bonder`, Python e plataforma;
- ferramentas externas disponiveis (SAINTexpress, Rscript, DIA-NN, ColabFold);
- decisao de device (GPU/CPU) e o motivo;
- proveniencia da ingesta (motor, aquisicao, tipo de quant, arquivos-fonte);
- status de cada passo (ok/pulado/erro) — cobertura honesta.

## 3. Sem instalacao silenciosa

`bonder/_deps.py`: dependencias faltantes levantam erro claro; auto-pip so com
`BONDER_AUTO_PIP=1`. Conteineres fixam as versoes (conda env + wheels).

## 4. Determinismo

- Ausente = NaN (nunca 0); MaxLFQ do DIA-NN nunca e re-normalizado.
- Cache de embeddings ESM2 (se usado) chaveado por (sha256 da sequencia, modelo).
- Sementes fixadas onde aplicavel.

## 5. Deposito de dados (exigido por revisores)

- Deposite os **.raw** e resultados de busca no **ProteomeXchange via PRIDE**
  (https://www.ebi.ac.uk/pride/).
- Inclua: parametros de busca, FASTA de busca (`data/search_db.fasta`), criterios
  de FDR, desenho de replicas/controles e as **versoes** (do `manifest.json`).
- Para AP-MS, reporte SAINT/MiST BFDR e limiares usados para chamar interatores.

## 6. Checklist de submissao

- [ ] `.raw` + saidas de busca no PRIDE (identificador ProteomeXchange)
- [ ] `manifest.json` de cada contraste anexado como material suplementar
- [ ] versoes de MaxQuant/FragPipe/DIA-NN, SAINTexpress, DEP/limma/MSstats
- [ ] tabela de interatores com SAINT-score/BFDR e (quando houver) log2FC/FDR
- [ ] rede de interacao (GraphML) e, se aplicavel, validacoes AF2-Multimer
