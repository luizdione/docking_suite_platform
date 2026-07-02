# Bonder - Plataforma de Analise de Interactoma (CoIP-MS/AP-MS)

## Projeto
Plataforma Dockerizada para analise de espectrometria de massas por co-imunoprecipitacao
(CoIP-MS/AP-MS), com foco no proteoma murino (Mus musculus, taxid 10090). Identifica
parceiros de interacao proteica usando ferramentas gold-standard da comunidade (SAINTexpress,
DEP/limma, MSstats, MiST). Nada e reimplementado - as estatisticas rodam via Rscript e
binarios compilados.

## Repositorios
- **Desenvolvimento**: `/home/user/docking_suite` (branch: `claude/bonder-proteome-platform-k8m6pa`)
- **Standalone**: `https://github.com/luizdione/docking_suite` (push direto)
- **Sync**: `sync.ps1` copia de docking_suite -> bonder-src local do usuario

## Estrutura do Codigo (`bonder/`)

### Core Python (`bonder/bonder/`)
- `io/` — parsers (MaxQuant, DIA-NN, FragPipe, generico), auto-detect motor, contaminantes
- `scoring/` — saint.py (SAINTexpress-int), mist.py (MiST host-pathogen)
- `stats/` — prepare.py + runner.py (DEP/limma via R), qc.py (somente QC, sem inferencia)
- `db/` — downloader.py (cache reproduzivel), uniprot.py, string_biogrid.py, crapome.py, baits.py
- `gpu/` — detect.py, gpu_gate.py (cooperativo c/ Ollama), tools.py, esm2.py (suplementar)
- `pipeline/` — orchestrator.py (coordenador unico), report.py (montagem + volcano + JSON)
- `web/` — app.py (Flask), jobs.py (RQ+Redis ou thread), health.py (selftest), docs.py
- `annotate/` — anotacao + rede (networkx -> cytoscape json)
- `cli.py` — CLI: db sync, baits check, qc, list-contrasts, selftest, run
- `config.py` — carregamento YAML com overrides de ambiente

### R Scripts (`bonder/R/`)
- `run_dep_limma.R` — DEP: make_se → filter → VSN → impute(MinProb) → test_diff → results
- `run_msstats.R` — MSstats: dataProcess → groupComparison

### Configs (`bonder/config/`)
- `baits.yaml` — 2 iscas exogenas: NSs_ICOARACI (viral) e LGALS3_HUMAN (humana P17931)
- `samples.yaml` — mapa de amostras + contrastes (NSs: 3 IP x 3 CTRL; LGALS3: 7 condicoes multi-infeccao)
- `databases.yaml` — versoes pinadas (UniProt, STRING, BioGRID, cRAP)
- `gpu.yaml` — configuracao GPU, gpu_gate, ESM2

### Docker (`bonder/docker/`)
- `Dockerfile.analysis` — mambaorg/micromamba + R/Bioconductor + SAINTexpress (compilado de SourceForge + GSL)
- `Dockerfile.gpu` — nvidia/cuda:12.4.1 + torch cu124
- `docker-compose.yml` — web + worker + redis
- `docker-compose.gpu.yml` — overlay GPU (RTX 3060 12GB)

### Testes (`bonder/tests/`)
21 testes: io parsers (9), scoring format (3), pipeline e2e (1), web jobs (1), health (3), docs (4)

## Iscas do Projeto Piloto
1. **NSs_ICOARACI** — proteina nao-estrutural do orthobunyavirus Icoaraci (viral, sem ortologo murino)
   - Scoring: SAINT + MiST (paradigma host-pathogen)
   - FASTA: `config/sequences/NSs_Icoaraci.fasta` (sequencia real ja inserida)
2. **LGALS3_HUMAN** — Galectina-3 humana (P17931), ortologo murino Lgals3 (P16110) mantido DISTINTO
   - Scoring: SAINT (IP vs IgG) + DEP/limma diferencial (infectado vs controle)
   - FASTA: `config/sequences/LGALS3_HUMAN.fasta` (baixar de UniProt: P17931)

## Desenho Experimental LGALS3 (7 condicoes)
- `IP_LGALS3_CTRL` — controle nao infectado (3 replicas)
- `IP_LGALS3_LA6h/LA18h/LA24h` — L. amazonensis 6h/18h/24h (3 replicas cada)
- `IP_LGALS3_TC6h/TC18h/TC24h` — T. cruzi 6h/18h/24h (3 replicas cada)
- `CTRL_IgG_LGALS3` — IgG controle negativo (3 replicas, celulas nao infectadas)
- Contrastes SAINT: cada condicao IP vs IgG
- Contrastes diferenciais (`type: differential`): cada infeccao vs CTRL — so DEP/limma, sem SAINT

## Convencoes Importantes
- MSDataset: accession-keyed, NaN para missing (NUNCA zero), DIA-NN MaxLFQ nunca re-normalizado
- Iscas exogenas NUNCA entram como preys no scoring
- P17931 (humana) != P16110 (murina) — NAO colapsar
- SAINTexpress compilado de tar.gz (SourceForge) + GSL; build e best-effort e nao-fatal
- Pyarrow com fallback pickle (pyarrow quebra em Windows/conda)
- Dashboard: Flask em localhost:8080, 6 abas (Resumo/QC, Interatores, Volcano, Rede, Manutencao, Ajuda)

## Como Rodar
```bash
# Docker (principal)
docker compose -f bonder/docker/docker-compose.yml up
# http://localhost:8080

# Demo (validacao)
docker compose -f bonder/docker/docker-compose.yml run --rm worker \
  bonder run --input /examples/demo/proteinGroups_demo.txt --contrast NSs_vs_CTRL

# Selftest
bonder selftest

# Listar contrastes disponiveis
bonder list-contrasts
```

## Scripts do Usuario (Windows)
- `run.ps1` / `run-gpu.ps1` — sobe o Docker (CPU / GPU)
- `stop.ps1` — para os containers
- `iniciar_bonder.ps1` — abre http://localhost:8080 no navegador
- `criar_atalho.ps1` — cria atalho na area de trabalho
- `sync.ps1` — sincroniza upstream para repo standalone
