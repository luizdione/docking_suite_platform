# bonder

**Plataforma local (navegador + Docker) para analise de CoIP-MS / AP-MS e
identificacao de parceiros de interacao no proteoma de camundongo**
(*Mus musculus*, taxid 10090).

`bonder` **nao reinventa a roda**: ele **orquestra ferramentas consagradas e
revisadas por pares** da comunidade de proteomica dentro de conteineres
reprodutiveis, com um dashboard Flask acessivel em `http://localhost:8080`.
Resultado: rigor e reprodutibilidade compativeis com publicacao, rodando **no seu
servidor**, sem software fechado/pago.

> Projeto-piloto: duas iscas **exogenas** ao camundongo, cada uma em uma Co-IP
> separada contra amostras murinas —
> **NSs** do ortobunyavirus *Icoaraci* (isca viral) e **Galectina-3** humana
> (LGALS3, UniProt P17931). Os parceiros buscados sao proteinas de camundongo.

---

## Por que open source tem rigor para publicar?

A proteomica de interacao publicada em alto impacto roda majoritariamente em
software **aberto**. O que confere publicabilidade nao e a marca, e sim: usar
metodos **validados e citaveis**, controle de **FDR** adequado, desenho com
replicas + controles, e **deposito dos dados brutos** (PRIDE/ProteomeXchange).
O interactoma host-pathogen do SARS-CoV-2 (Gordon et al., *Nature* 2020) e o
BioPlex foram feitos com estas mesmas ferramentas.

## Ferramentas orquestradas (todas com paper)

| Etapa | Ferramenta (aberta) |
|---|---|
| Busca/quantificacao (raw -> proteina) | **MaxQuant** (DDA) · **FragPipe/MSFragger** · **DIA-NN** (DIA, GPU) |
| **Especificidade isca-presa** | **SAINTexpress** (todas as iscas) · **MiST** (host-pathogen, isca viral) |
| Quantitativo diferencial | **DEP / limma / MSstats** (R/Bioconductor) |
| Background / contaminantes | **CRAPome** interno (controles) + cRAP |
| Anotacao / redes | **UniProt, STRING, BioGRID**, GO (g:Profiler), Cytoscape |
| GPU (RTX 3060) | **DIA-NN**, **MS2Rescore + AlphaPeptDeep** (rescoring), **AlphaFold2-Multimer/ColabFold** (validacao estrutural) |

> ESM2 (embeddings) entra apenas como **priorizacao exploratoria suplementar**,
> nunca como base do escore de confianca que vai para a figura principal.

## Arquitetura

```
 navegador (localhost:8080)
        |
   [web] Flask ---- fila (Redis/RQ) ----> [worker]
                                             |
   io -> db -> scoring(SAINT/MiST) -> stats(DEP/limma R) -> annotate -> pipeline -> relatorio
                                             |
                        [gpu] (opcional) DIA-NN / rescoring / AlphaFold / ESM2
```

- **Fluxo unidirecional** `io -> db -> scoring/stats -> gpu -> pipeline -> web`.
- O **nucleo cientifico e testavel em CPU**, sem GPU e sem R (ver `tests/`).
- Etapas cujas ferramentas externas faltam sao **puladas com registro no
  manifesto** — o pipeline nunca finge cobertura que nao teve.

## Como rodar (rapido)

```bash
# 1) CPU (funciona em qualquer maquina):
docker compose -f docker/docker-compose.yml up --build
#    -> abra http://localhost:8080

# 2) Com GPU (RTX 3060; requer NVIDIA Container Toolkit):
docker compose -f docker/docker-compose.yml -f docker/docker-compose.gpu.yml up --build
```

Coloque suas saidas de busca em `./inputs`. No dashboard: informe o caminho
(`/inputs/proteinGroups.txt`), escolha o contraste (isca vs controle) e rode.
O **motor e o modo de aquisicao (DDA/DIA) sao detectados automaticamente**.

Sem Docker (conda `datascience`):

```bash
conda env update -n datascience -f environment.yml
pip install -e .
python -m bonder.cli list-contrasts
python -m bonder.web.app          # dashboard em 127.0.0.1:8080
```

## Configuracao (declarativa)

- `config/baits.yaml` — as duas iscas (viral NSs, humana LGALS3), com ortologo
  murino `Lgals3` (P16110) mantido **distinto** da isca humana.
- `config/samples.yaml` — amostras, replicas e **contrastes** isca-vs-controle.
- `config/databases.yaml` — bancos murinos **fixados por versao** + checksums.
- `config/gpu.yaml` — device, gpu_gate, e toggles das etapas de GPU.
- `config/sequences/` — FASTA das iscas (veja o README de la).

## CLI

```bash
bonder db sync            # baixa/verifica UniProt/STRING/contaminantes (murinos)
bonder baits check        # valida iscas e monta o search_db.fasta
bonder qc --input F       # resumo de QC de uma saida de busca
bonder run --input F --contrast NSs_vs_CTRL
```

## Reprodutibilidade

Cada corrida grava `manifest.json` (versoes, ferramentas disponiveis, device,
proveniencia, passos) e `results.json`. Bancos sao **fixados por versao** e
verificados por sha256. Ver `docs/REPRODUTIBILIDADE.md` e o guia de deposito no
PRIDE.

## Testes

```bash
python -m pytest -q     # nucleo CPU (IO, SAINT-format, pipeline e2e, web jobs)
```

## Requisitos alvo

Windows 11 · RTX 3060 12GB · 32GB RAM · miniforge3 (env `datascience`) ·
Python 3.11 · Docker Desktop + WSL2 + NVIDIA Container Toolkit (para GPU).

## Licenca

MIT (ver `LICENSE`). As ferramentas orquestradas mantem suas proprias licencas;
cite-as conforme `docs/METODOLOGIA.md`.
