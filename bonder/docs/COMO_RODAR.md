# Como rodar o bonder (passo a passo)

## Pre-requisitos

- Docker Desktop (Windows 11 + WSL2) ou Docker Engine (Linux).
- Para GPU: NVIDIA Container Toolkit + driver compativel com CUDA 12.4.
- Alternativa sem Docker: miniforge3 com o env `datascience`.

## 1. Preparar iscas e bancos

```bash
# sequencias das iscas (ver config/sequences/README.md)
curl -L "https://rest.uniprot.org/uniprotkb/P17931.fasta" -o config/sequences/LGALS3_HUMAN.fasta
# (adicione config/sequences/NSs_Icoaraci.fasta do seu isolado)

# baixar bancos murinos fixados por versao (uma vez)
python -m bonder.cli db sync
# validar iscas + montar o FASTA de busca combinado
python -m bonder.cli baits check
```

## 2. (Opcional) Gerar a saida de busca a partir dos .raw

Se voce so tem os arquivos brutos:

- **DDA:** rode MaxQuant ou FragPipe usando `data/search_db.fasta` -> obtenha
  `proteinGroups.txt` (MaxQuant) ou `combined_protein.tsv` (FragPipe).
- **DIA:** rode DIA-NN (GPU) -> obtenha `report.pg_matrix.tsv` (+ `report.tsv`).
  No conteiner GPU:
  ```bash
  docker compose -f docker/docker-compose.yml -f docker/docker-compose.gpu.yml \
    run --rm gpu python -m bonder.gpu.tools   # wrappers DIA-NN/rescoring
  ```

Coloque a saida em `./inputs/`.

## 3. Subir a plataforma

```bash
# CPU (sempre funciona)
docker compose -f docker/docker-compose.yml up --build
# GPU (RTX 3060)
docker compose -f docker/docker-compose.yml -f docker/docker-compose.gpu.yml up --build
```

Abra **http://localhost:8080**.

## 4. Rodar uma analise (dashboard)

1. Informe o caminho da saida de busca, ex.: `/inputs/proteinGroups.txt`.
2. Escolha o contraste (ex.: `NSs_vs_CTRL` ou `LGALS3_vs_CTRL`).
3. Clique **Rodar analise**. Acompanhe o status.
4. Veja **Resumo/QC**, **Interatores** (tabela ranqueada por BFDR/FDR),
   **Volcano** (se DEP rodou) e **Rede** (cytoscape).

## 4b. Rodar via CLI (headless / CI)

```bash
python -m bonder.cli run --input inputs/proteinGroups.txt --contrast NSs_vs_CTRL
# saidas em runs/<run_id>/: interactors.tsv, network.graphml, volcano.png,
#                           results.json, manifest.json
```

## 5. Saidas

- `interactors.tsv` — tabela ranqueada (prey, gene, SAINT/BFDR, log2FC/FDR,
  background, confianca).
- `network.graphml` / `network.json` — rede de interacao (Cytoscape/cytoscape.js).
- `volcano.png` — volcano do DEP/limma.
- `manifest.json` — reprodutibilidade (versoes/ferramentas/device/passos).

## Aba Manutencao (auditoria + autoteste)

No dashboard, a aba **Manutencao** faz uma auditoria de sanidade:

- dependencias Python (nucleo x opcionais) e suas versoes;
- ferramentas externas no PATH (SAINTexpress, Rscript, DIA-NN, ColabFold, nvidia-smi);
- pacotes R gold-standard (DEP/limma/MSstats) quando o Rscript existe;
- GPU (device escolhido + VRAM livre) e config (iscas/contrastes/avisos);
- bancos de referencia em cache + `search_db.fasta`; caminhos com escrita; corridas.

O botao **Rodar autoteste** executa o pipeline em dados sinteticos (num diretorio
temporario) e mostra o status de cada etapa e verificacao — a checagem de
sanidade dos scripts e processos, sem depender de GPU nem de R. Tambem por API:

```bash
curl http://localhost:8080/api/health/full
curl -X POST http://localhost:8080/api/health/selftest
```

## Dicas de GPU (compartilhar com o Ollama)

O `gpu_gate` espera VRAM livre e ausencia de processo do LLM (Ollama/llama) antes
de carregar etapas de GPU, no mesmo espirito da suite docking. Ajuste limiares em
`config/gpu.yaml`.
