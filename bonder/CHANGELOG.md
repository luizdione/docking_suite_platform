# Changelog - bonder

Todas as mudancas relevantes da plataforma. Formato baseado em
[Keep a Changelog](https://keepachangelog.com/pt-BR/1.1.0/); versionamento
[SemVer](https://semver.org/lang/pt-BR/).

## [Nao lancado]

### Adicionado
- Aba **Ajuda e Tutoriais** no dashboard: serve os documentos `.md` (tutoriais,
  metodologia, reprodutibilidade, versionamento) e a licenca, renderizados no
  navegador. Endpoints `GET /api/docs` e `GET /api/docs/<id>`.

## [0.1.0] - 2026-07-02

Primeira versao do subprojeto `bonder` (plataforma CoIP-MS / AP-MS para
identificacao de parceiros no proteoma murino).

### Adicionado
- **Ingesta (`bonder.io`)**: parsers MaxQuant/DIA-NN/FragPipe/generico para um
  `MSDataset` unificado; deteccao automatica de motor e aquisicao (DDA/DIA);
  ausente = NaN (nunca 0); chave por acesso UniProt lider; iscas exogenas
  marcadas e mantidas distintas do ortologo murino.
- **Especificidade (`bonder.scoring`)**: formatacao de entrada + wrappers de
  **SAINTexpress** (todas as iscas) e **MiST** (isca viral NSs, host-pathogen).
- **Quantitativo (`bonder.stats`)**: wrappers R para **DEP/limma/MSstats**
  (Bioconductor) via Rscript; QC em Python (sem valor inferencial). MaxLFQ do
  DIA-NN nunca e re-normalizado.
- **Bancos (`bonder.db`)**: UniProt/STRING/BioGRID/cRAP murinos com cache
  reproduzivel (checksum, modo offline) + montagem do `search_db.fasta`.
- **GPU (`bonder.gpu`)**: DIA-NN / MS2Rescore / AlphaPeptDeep / AlphaFold e ESM2
  (suplementar), com fallback CPU e `gpu_gate` cooperativo (compartilha a placa
  com o Ollama).
- **Pipeline + Web**: orquestrador, relatorio (tabela/volcano/rede), dashboard
  Flask + fila RQ/Redis (fallback thread).
- **Docker**: imagem de analise (Python + R + SAINT) e imagem GPU (CUDA 12.4 +
  torch); `docker compose` CPU por padrao + overlay de GPU (RTX 3060).
- **Aba Manutencao**: auditoria (deps, ferramentas, R, GPU, config, bancos,
  caminhos, corridas) + autoteste do pipeline em dados sinteticos.
- Documentacao PT-BR (metodologia com citacoes, como rodar, reprodutibilidade)
  e testes do nucleo (CPU, sem GPU/R).

### Corrigido
- `saint` nao existe como pacote conda: removido do env; SAINTexpress agora e
  compilado do codigo-fonte (SourceForge) no Dockerfile, best-effort e nao-fatal.
- Mapeamento de peptideos/q-value do DIA-NN pela `Protein.Group` original.
- Fallback de gene no SAINT gravava o literal `nan` (NaN e truthy) -> cai para o acesso.
- Dashboard: coluna log2fc vazia sem R -> usa `log2fc_descriptive`; o polling
  agora para em status terminal (done/error).
- CI: dependencia `requests` faltando no job de testes.
