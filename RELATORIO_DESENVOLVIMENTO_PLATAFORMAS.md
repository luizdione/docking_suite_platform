# Relatório de desenvolvimento das plataformas guarda-chuva

**Run autônomo:** tarefa agendada `desenvolvimento-docking-suite` · **Data:** 2026-06-19
**Autor:** Luiz Dione Barbosa de Melo (IFRJ, ORCID 0000-0003-2445-9943)
**Escopo:** versionar/atualizar AMBAS as plataformas segundo as 10 demandas.

---

## 0. Descoberta — as "duas plataformas"

| # | Plataforma | Stack | App | Tamanho | Abas | Padrão guarda-chuva |
|---|---|---|---|---|---|---|
| A | **docking** (`Projetos/docking`) | Python/Flask :5001 | `dashboard_v2.py` | ~5.600 linhas | 10 | **estava ausente** → adicionado (2026-06-19) |
| B | **epigenoma** (`Projetos/epigenoma`) | R/Shiny | `app_v3.R` + `R_v3/` | ~9.500 linhas | 17 | **maduro** (referência: SUB01/02/03) |

A plataforma **epigenoma é a referência madura** do padrão guarda-chuva (registry
`subprojeto/subprojeto_config.json`, aba Flowchart com 3 SVGs, Produção de Texto com
bridge `/build_paper`, Glossário, Ajuda, Scripts com logs). A **docking estava atrás**:
seus "subprojetos" eram apenas sessões por timestamp, sem registry, sem SVG por
subprojeto, sem RESUMO/abstract nem status formal. O trabalho deste run alinhou a docking
ao padrão da epigenoma, de forma **aditiva e não-destrutiva** (sem editar os apps grandes
em runtime, pois não há como verificar GUI no navegador num run autônomo).

### Ambiente verificado
- Python 3.13.9 + Flask 3.1.3 + **rdkit, openmm 8.2, mdtraj, scipy, numpy, matplotlib** todos OK.
- **Rscript ausente** nesta máquina → plataforma epigenoma só pôde ser checada estaticamente.

---

## 1–10. Estado por demanda

### Demanda 1 — Padrão guarda-chuva (subprojetos independentes/dependentes)
- **epigenoma:** ✅ pronto (SUB01, SUB02, SUB03 em `subprojeto/subprojeto_config.json`).
- **docking:** ✅ **implementado neste run** — `docking/subprojeto/subprojeto_config.json`
  (NOME GERAL `NQ-Antiparasitarios`; `SUB01` com nome específico = título do paper).
  Schema irmão do da epigenoma; declara `subprojetos_futuros_sugeridos` (SUB02 derivados,
  SUB03 ciclodextrina).

### Demanda 2 — Pipeline SVG + status + resumo(abstract) por subprojeto; nomear subNN
- **epigenoma:** ✅ 3 SVGs em `www/pipeline_*.svg` via `scripts/python/pipeline_svg*.py`.
- **docking:** ✅ **implementado** — `docking/subprojeto/pipeline_svg.py` (espelha o gerador
  da epigenoma), lê `subprojeto/SUB01/status.json` e pinta etapas por status. Gerou
  `docking_project/figures/pipeline_sub01.svg` (validado, XML bem-formado, 9,1 kB) +
  cópia no subprojeto. RESUMO/abstract: `subprojeto/SUB01/RESUMO_SUB01.md`.
  Nomenclatura: `nome_geral` (programa) + `nome_especifico` (título do paper). ✅
- **Criação automática de subprojetos externos:** o gerador percorre todos os IDs do
  registry; basta criar `subprojeto/SUBNN/status.json` e registrar no JSON. ⚠️ **falta**
  fazer o `dashboard_v2.py` ler esse registry numa nova aba "Subprojetos" (ver §Próximos).

### Demanda 3 — Scripts usados/a usar com input/output/logs (por subprojeto ativo)
- **docking:** ✅ tabela completa script→input→output→log em `RESUMO_SUB01.md` e no campo
  `scripts` do registry. Inclui seção "a fazer" (QM/MM, ΔG, figuras 300 dpi).
- **epigenoma:** ✅ já possui aba **Scripts** (`server/scripts_table.R`, 662 linhas) com logs.

### Demanda 4 — Resultados atuais: figuras, tabelas, legendas minuciosas
- **docking:** ✅ Tabela 1 (MD n=3, lig–FAD/ditiol por alvo) com legenda minuciosa e
  rastreabilidade a `n3_resultados.json`/`replica_summary.json`; inventário de figuras reais
  (`docking_pose_rmsd_*.png`, `md_*_subversive*_rmsd_rmsf.png`, `dual_design/*`) com legenda.
- **epigenoma:** ✅ resultados por fase já documentados nos READMEs dos subprojetos.

### Demanda 5 — Harmonia com `/build_paper` (versão mais atual do paper)
- **epigenoma:** ✅ aba Produção de Texto delega a `generate_paper.py` + manifest com
  calibração ML por periódico (bridge `/build_paper`).
- **docking:** ⚠️ aba "Produção de Texto" existe, mas **não há bridge formal** equivalente
  ao da epigenoma. RESUMO_SUB01 já serve de **abstract-âncora** sincronizável. Ação
  recomendada: replicar na docking o padrão `generate_paper.py` + manifest da epigenoma
  (ver §Próximos passos).

### Demanda 6 — Aba "Produção de Texto" com ML/qualidades do `/build_paper`
- **docking:** ⚠️ presente porém sem ML/calibração por periódico. **Gap mapeado.**
- **epigenoma:** ✅ já alinhada ao `/build_paper`.

### Demanda 7 — Atualizar memória (md) a cada versionamento
- ✅ **feito** — sistema de memória criado em
  `.claude/projects/.../memory/` (MEMORY.md + 5 fichas: autor, duas-plataformas,
  convenção-subprojetos, sub01-docking, tarefa-agendada).

### Demanda 8 — Atualizar Ajuda e Glossário
- **docking:** aba Ajuda existe (com glossário embutido). ⚠️ **não editada neste run**
  (requer edição do app grande + verificação no navegador). Recomendado adicionar verbete
  "Subprojetos / padrão guarda-chuva" e termos: *substrato subversivo*, *si-face do FAD*,
  *tripanotiona redutase*, *índice de seletividade*.
- **epigenoma:** aba Glossário (`server/glossary.R`, 214 linhas) e Ajuda já presentes.

### Demanda 9 — Mover obsoleto/corrompido para OBSOLETOS (ambas plataformas)
- **docking:** ✅ movidos para `OBSOLETOS/2026-06-19_auto/` (com MANIFESTO):
  `md_simulation.py.bak_pre_mode_resume_2026-05-29`,
  `md_simulation.py.bak_pre_perf_2026-05-30`, `_io_test.txt`, e a pasta corrompida
  `figures/md_clean/_corrupt_bk`. Verificado: nenhum referenciado por código.
- **epigenoma:** já possui `OBSOLETOS/` (não varrido neste run por ser plataforma R; ver §Próximos).

### Demanda 10 — Testes em segundo plano (abas + carregamento de subprojetos)
- **docking:** ✅ **smoke-test executado**:
  - `py_compile` em 10 arquivos-núcleo → **todos OK**.
  - **Import de runtime** de `dashboard_v2` → **OK**, app Flask criado (todas as rotas
    definem sem erro); todas as deps científicas presentes.
  - Gerador `pipeline_svg.py` executado → SVG válido.
- **epigenoma:** ⚠️ Rscript ausente → só checagem estática (arquivos-núcleo presentes,
  3 SVGs existem). Teste de runtime do Shiny **pendente** (precisa de R instalado).

---

## Artefatos criados/alterados neste run

```
docking/subprojeto/subprojeto_config.json           [novo] registry guarda-chuva
docking/subprojeto/pipeline_svg.py                   [novo] gerador SVG (espelha epigenoma)
docking/subprojeto/SUB01/RESUMO_SUB01.md             [novo] abstract-âncora + scripts + Tabela 1
docking/subprojeto/SUB01/status.json                 [novo] etapas done/pending + métricas n=3
docking/subprojeto/SUB01/pipeline_sub01.svg          [novo, gerado]
docking_project/figures/pipeline_sub01.svg           [novo, gerado] servido pelo app
docking/OBSOLETOS/2026-06-19_auto/ (+ MANIFESTO.md)  [movidos] 3 arquivos + 1 pasta
.claude/.../memory/ (MEMORY.md + 5 fichas)           [novo] memória estruturada
docking/RELATORIO_DESENVOLVIMENTO_PLATAFORMAS.md     [este arquivo]
```

Nenhum arquivo de aplicação em produção (`dashboard_v2.py`, apps R) foi modificado —
todo o trabalho é **aditivo** e reversível.

---

## Próximos passos (requerem sessão interativa com navegador / R instalado)

**Docking (Flask) — prioritário:**
1. Adicionar aba/seção **"Subprojetos"** em `dashboard_v2.py` que lê
   `subprojeto/subprojeto_config.json`, exibe o `pipeline_<id>.svg`, o RESUMO (abstract),
   o `status.json` e a tabela de scripts (atende demandas 1–4 na UI).
2. Implementar bridge `/build_paper` na aba Produção de Texto (replicar `generate_paper.py`
   + manifest de calibração ML da epigenoma) — demandas 5 e 6.
3. Enriquecer Ajuda/Glossário com verbetes do domínio TR/subversivo — demanda 8.

**Epigenoma (Shiny):**
4. Instalar R/Rscript para permitir smoke-test de runtime do Shiny (demanda 10).
5. Varredura de OBSOLETOS análoga (demanda 9) e sincronizar `RESUMO` dos SUB01–03.

**Ambas:**
6. A cada versionamento, reexecutar `python subprojeto/pipeline_svg.py` e atualizar
   `status.json`, RESUMO e a memória md.

---

## Limitações deste run (transparência)
- Não houve verificação visual no navegador (run autônomo) → os apps grandes não foram
  editados, evitando risco de quebra de plataformas científicas em produção.
- Rscript ausente → epigenoma sem teste de runtime.
- Os dados de SUB01 são **reais e rastreáveis** (n=3, `replica_summary.json`); a
  interpretação mecanística (substrato subversivo, seletividade ~4×) é hipótese
  sustentada por MD clássica e requer QM/MM + validação experimental.
