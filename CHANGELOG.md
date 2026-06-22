# CHANGELOG — Docking / MD TcTR + 2-SCN-NQ

Formato baseado em [Keep a Changelog](https://keepachangelog.com/pt-BR/1.1.0/) e
[Semantic Versioning](https://semver.org/lang/pt-BR/). Datas em horário de
Brasília (BRT). Pendências e contexto vivo continuam em `docs/BRIEFING_*.md`.

---

## [2.4.1] — 2026-06-22

Fase 2 da curadoria: correção dos bugs de backend e de figuras mapeados no *Roadmap* da
2.4.0, cada um validado em runtime (Flask `test_client`, sem rede/GPU). Um bug adicional de
serialização foi descoberto durante o *smoke test* e corrigido na mesma rodada.

### Fixed
- **Paper AI — autenticação Anthropic** (`/generate_paper_text`): agora envia os headers
  `x-api-key` e `anthropic-version` e lê `ANTHROPIC_API_KEY` do ambiente, com *guard* e
  mensagem clara quando ausente (antes: 401 cego capturado pelo `except` → "API error").
- **`reload(None)` de `md_module`** (`/run_all_figures`, `/generate_all_sections`): o módulo
  aposentado deixa de cair em `reload() argument must be a module`; retorna `[SKIP]`/`[deprecated]`.
- **`results.json` por sessão** (`/generate_paper_text`, `/static_text`): resolvido via
  `compound_db.get_session_dirs` (caminho por composto) em vez de um *path* fixo inexistente
  que sempre forçava o texto ao *fallback* hardcoded; quando não há dados, registra
  `app.logger.warning` em vez de engolir a falha em silêncio.
- **Fallback de `matplotlib`** (`admet_module`, `cyclodex_module`): o `except ImportError`
  passa a reimportar todos os símbolos do `try` (`FancyBboxPatch`, `GridSpec`, `mpatches`),
  evitando `NameError` numa primeira execução sem matplotlib (CI/máquina limpa).
- **Serialização JSON da ciclodextrina** (`cyclodex_module.predict_cd_inclusion`): com RDKit,
  os critérios de viabilidade e as dimensões eram `numpy.bool_`/`numpy.float64`, quebrando
  `run_cyclodextrin` com *"Object of type bool is not JSON serializable"*; convertidos para
  tipos nativos (`bool`/`float`) na origem, deixando o `dict` 100% JSON-safe.

### Changed
- `VERSION`: 2.4.0 → 2.4.1.

---

## [2.4.0] — 2026-06-21

Auditoria e curadoria abrangente da plataforma (multiagente, 86 achados) + registro do
dataset científico **n=3** (três redutases) e do **controle inibidor TcTR a 50 ns**, já
propagados ao manuscrito (ACS Infect. Dis.). Preparação para publicação no GitHub
(repositório privado `docking_suite`, MIT).

### Added
- `LICENSE` (MIT, 2026) e `.gitignore` consolidado (exclui `PAPER/`, venv, PyMOL,
  trajetórias `.dcd` de GB, logs, segredos e o manuscrito).
- `paper_tools/` (fora do repo) — scaffolding do manuscrito separado da plataforma
  (`edit_tex_*`, `edit_bib_refs`, `finalize_md_n3`, `rerun_*`, `n3_resultados.json`…) + MANIFESTO.
- `docking_project/notion_config.example.json` — modelo de config (sem token real).

### Changed
- `VERSION`: 2.3.0 → 2.4.0.
- `requirements.txt`: + `pdbfixer` (usado por `md_simulation.py`).
- `README.md`: seção MD/"Achado atual" reescrita para o estado real (subversivo n=3 ~50 ns
  como hipótese central; inibidor 50 ns como controle); removida a árvore `fig13_3d_kit/` (inexistente).

### Fixed
- `toxicity_module._to_float`: `bare except` → `except (ValueError, TypeError)`.
- `compound_db`: variável de exceção morta no fallback RDKit/obabel.

### Housekeeping
- Diagnósticos mortos (`_fad_diagnose`, `diag_fad_bonds`, `diag_warmup_nan`) → `OBSOLETOS/2026-06-21_audit/`.
- Veredito de segurança da auditoria: **zero segredos hardcoded** (notion/Anthropic leem de env/arquivo).

### Roadmap (da auditoria — próximas releases, requerem teste de runtime)
- Bugs do app: auth Anthropic (401), `reload(None)` de md_module, `results.json` no path da sessão, fallback matplotlib.
- Melhorias: bridge `/build_paper`, botão EN+PT, Glossário na Ajuda, `_deps.py` unificado, logging nos `except` silenciosos.

---

## [2.3.0] — 2026-06-19

Camada **guarda-chuva de subprojetos** adicionada à plataforma docking, alinhando-a ao
padrão maduro da plataforma irmã epigenoma (R/Shiny). Run autônomo da tarefa agendada
`desenvolvimento-docking-suite`. Detalhes completos em
`RELATORIO_DESENVOLVIMENTO_PLATAFORMAS.md`.

### Added
- `subprojeto/subprojeto_config.json` — registry guarda-chuva (NOME GERAL
  `NQ-Antiparasitarios`; `SUB01` = título atual do paper; schema irmão da epigenoma).
- `subprojeto/SUB01/RESUMO_SUB01.md` — abstract-âncora + tabela de resultados n=3 com
  legenda minuciosa + tabela de scripts (input/output/logs).
- `subprojeto/SUB01/status.json` — etapas (done/pending/blocked) + métricas-chave n=3.
- `subprojeto/pipeline_svg.py` — gerador de SVG de pipeline por subprojeto (espelha
  `epigenoma/scripts/python/pipeline_svg.py`), lê `status.json` e pinta por status.
- `docking_project/figures/pipeline_sub01.svg` (+ cópia em `subprojeto/SUB01/`) — gerado.
- Memória estruturada em `.claude/projects/.../memory/` (índice + 5 fichas).
- `RELATORIO_DESENVOLVIMENTO_PLATAFORMAS.md` — mapeamento das 10 demandas × 2 plataformas.
- **UI Subprojetos no `dashboard_v2.py`** (2ª passagem do mesmo run): rotas Flask aditivas
  `GET /subprojects` (página-resumo), `GET /api/subprojects`, `GET /api/subproject/<ID>` e
  `GET /api/subproject/<ID>/pipeline.svg` — consomem o layer autoritativo `subprojeto/`
  (config + `status.json` + `RESUMO_<ID>.md`) e servem o SVG de `subprojeto/pipeline_svg.py`
  (regenerando sob demanda). Botão "Subprojetos" adicionado à barra de abas. Verificado
  em servidor Flask vivo (HTTP 200 em todas as rotas; abas existentes intactas).
- `docs/SUBPROJETOS_GUARDA_CHUVA.md` — documentação do padrão guarda-chuva (surge na aba Ajuda).

### Changed
- `VERSION`: 2.1.0 → 2.3.0 (corrige defasagem frente ao CHANGELOG, que já estava em 2.2.0).

### Validated (smoke-test, demanda 10)
- `py_compile` em 10 arquivos-núcleo da docking → todos OK.
- Import de runtime de `dashboard_v2` → OK, app Flask instanciado (rotas sem erro); deps
  científicas (rdkit, openmm, mdtraj, scipy...) presentes em Python 3.13.9 + Flask 3.1.3.
- `subprojeto/pipeline_svg.py` executado → SVG XML bem-formado.
- Epigenoma: checagem estática (arquivos-núcleo presentes, 3 SVGs). Runtime do Shiny
  **pendente** — Rscript ausente nesta máquina.

### Housekeeping
- Movidos para `OBSOLETOS/2026-06-19_auto/` (com `MANIFESTO.md`; sem referências em código):
  `md_simulation.py.bak_pre_mode_resume_2026-05-29`,
  `md_simulation.py.bak_pre_perf_2026-05-30`, `_io_test.txt`,
  e a pasta corrompida `docking_project/figures/md_clean/_corrupt_bk`.

### Pendente (requer sessão interativa / R)
- ~~UI: aba "Subprojetos" em `dashboard_v2.py`~~ — **FEITO** nesta versão (rotas Flask + página).
- Bridge `/build_paper` + calibração ML na aba "Produção de Texto" da docking.
- Enriquecer Glossário (verbetes TR/substrato subversivo) — doc geral já adicionado em `docs/`.
- Instalar R para smoke-test de runtime do Shiny.

---

## [2.2.0] — 2026-05-28

Sessão de produção bem-sucedida em RTX 3060 (5 ns smoke) + dois bugs críticos
na função `analyze()` corrigidos. Detalhes completos em
`docs/BRIEFING_2026-05-28_md_module.md` (§U1–U8).

### Added
- `run_md_safe.ps1` (criado em 05-26): wrapper PowerShell de rodada segura
  com gate `validate_ff.py` + power-cap (`nvidia-smi -pl 120`) + serialização
  do RNA-seq antes da MD. Usado pela primeira vez com sucesso em 05-28.
- `…/md_clean/md_prod_v5_smoke.{dcd,csv,chk}`: trajetória 5 ns @ 10 ps/frame
  (500 frames, 1,10 GB), CSV de StateData e checkpoint retomável.
- `…/md_clean/minimized_prod_v5_smoke.pdb`: topologia minimizada (15 MB) para
  análise pós-rodada.
- `CHANGELOG.md` (este arquivo).

### Changed
- `md_simulation.py` — função `analyze()`:
  - **Bug A (PBC):** RMSD do ligante "explodia" para 175 Å quando o ligante
    atravessava a borda da caixa cúbica (lado 124,15 Å). Aplicado
    *minimum-image convention* frame-a-frame só no ligante (mantém
    proteína/águas intactas). Validado: `max |shift| = 122,58 Å` (≈ 1
    caixa). RMSD pós-fix passou de 175,30 Å → 5,84 Å.
  - **Bug B (eixo X 10× menor):** `traj.time` do mdtraj a partir do DCD do
    OpenMM 8.2 + plugin CUDA 110 vinha com fator 10× errado (provável
    `DELTA` em AKMA no header). `analyze()` agora aceita
    `frame_interval_ps: float = 10.0` e computa o tempo explicitamente:
    `time_ns = (np.arange(1, n_frames+1) * frame_interval_ps) / 1000.0`.
    Pós-fix: `t_total_ns` passou de 0,499 → 5,0.
- `md_simulation.py` — docstring de cabeçalho atualizado (referenciava
  `md_clean.py … ~290 linhas` desde a renomeação de 05-25; agora descreve o
  pipeline real).

### Fixed
- Causa do corte de produção identificada e neutralizada (**Hipótese A
  confirmada como única, B descartada**): power-cut + Kernel-Power 41 +
  RNA-seq concorrente. Com `nvidia-smi -pl 120` + RNA-seq parado, a MD
  passou o boundary "produção iniciando" e rodou 5 ns sem hiccup
  (38,2 ns/dia, RTX 3060 @ 120 W).
- Item carry-over `DLL load failed 126` (CUDA) confirmado **resolvido** em
  runtime: `cudart64_110.dll` + plugins ok, plataforma SELECIONADA: CUDA
  (single precision).

### Validated
- Gate `validate_ff.py` rodou pela primeira vez em rodada real (passo `[2]`
  do `run_md_safe.ps1`):
  - **FAD:** 86 átomos, 91 bonds, 248 torsions, **soma de cargas = 0,000**,
    `|q|max = 0,48`. `RESULTADO: PASS`. (Nota: literatura sugere FAD
    oxidado a pH 7,4 ≈ −2; o XML canônico atual está neutro — decisão
    científica a discutir.)
  - **Ligante 2-SCN-NQ:** 20 átomos, 21 bonds, 57 torsions, soma de
    cargas = −0,000, `|q|max = 0,29`. `RESULTADO: PASS`.

### Housekeeping
- Backups movidos para `obsoletos/`:
  - `md_clean.py.bak_2026-05-21`
  - `md_clean.py.truncated_2026-05-23.bak`
- `*.broken_bak` de `docking_project/figures/smoke_test_fad/` movidos para
  `obsoletos/smoke_test_fad_broken_bak/`:
  - `cofactor_FAD_TcTR.{pdb,xml}.broken_bak`
  - `lig_TcTR{.pdb,_smirnoff.xml}.broken_bak`
  - `cofactor_FAD_TcTR.pdb.bak_2026-05-21`

### Observado, ainda pendente
- **Ligante NÃO converge em 5 ns** — RMSD sobe monotonicamente 0 → 5,69 Å
  (corpo rígido, `rmsd_lig_intra = 0,34 Å`). Sinal de instabilidade da pose
  Vina; não conclusivo. **Requer ≥ 20 ns** para discriminar
  unbinding / pose alternativa / flutuação transitória.
- Proteína: ✅ converge em ~2,18 Å, platô a partir de ~2,5 ns.
- `git init` ainda não feito (recomendação do briefing — não aplicado nesta
  sessão automatizada para preservar autonomia do usuário sobre `.gitignore`
  e commit inicial).

---

## [2.1.0] — 2026-05-25

Marco da renomeação `md_module.py → md_clean.py → md_simulation.py` e
regeneração canônica de force fields. Conteúdo reconstruído a partir do
`docs/BRIEFING_2026-05-25_md_module.md` (não houve CHANGELOG na época).

### Added
- `regen_fad_xml.py`, `regen_ligand_xml.py`: scripts para regenerar XML do
  FAD e do ligante a partir de SMILES + PDB canônicos.
- `validate_ff.py`: gate de validação rápida de force-field (sai 0=PASS,
  1=FAIL, 2=XML ausente) — base do passo `[2]` do `run_md_safe.ps1`.
- FAD canônico v3 (`…/smoke_test_fad/cofactor_FAD_TcTR.{xml,pdb}`,
  70.954 / 8.913 B) e ligante regenerado
  (`lig_TcTR_smirnoff.xml`, 15.353 B).
- `md_clean.py` (961 B): shim que delega para `md_simulation.main()`.

### Changed
- `md_module.py` → `md_clean.py` → `md_simulation.py` (~48 KB, 1.016 linhas).

### Fixed
- NaN no warm-up: FAD v3 com geometria limpa + ligante regenerado
  eliminaram bonds anômalas que disparavam NaN no `[2/5] Warm-up`.
- Produções curtas (`md_prod_v3 / v4`, 0,019 ns) confirmaram **0 NaN** e
  protocolo são (RMSD ligante 1,68 / 2,88 Å, RMSF médio 0,73 / 0,74 Å).

---

## Versões anteriores

Anteriormente à 2.1.0, as alterações estão registradas exclusivamente em
`docs/BRIEFING_*.md` (séries diárias 05-08, 05-20, 05-21, 05-23, 05-23
UPDATE, 05-25, 05-26 e 05-28). A partir da 2.2.0, novas releases entram
aqui em paralelo aos briefings.
