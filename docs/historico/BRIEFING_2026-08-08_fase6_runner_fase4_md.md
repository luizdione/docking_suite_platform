# Contexto da Sessão — Fase 6 (runner/proveniência) + início da Fase 4 (correções md_simulation)

- **Data/hora:** 2026-08-08 00:22 -03:00
- **Branch:** feat/guarda-chuva-multialvo
- **Repositório:** docking_suite_platform (C:\github\docking_suite_platform)
- **Foco desta sessão:** Concluir a Fase 6 (runner + proveniência) e iniciar a Fase 4 (correções não-negociáveis de md_simulation.py)

## 1. Objetivo / tarefa atual
Generalizar a plataforma de docking/MD para **guarda-chuva multi-alvo**, sem destruir o SUB01 nem apagar a peculiaridade de cada par proteína-ligante. Plano de 10 fases aprovado em `C:\Users\luizd\.claude\plans\analisando-a-plataforma-http-localhost-5-serene-steele.md`. Esta sessão executou a **Fase 6 inteira** e as **correções #1/#3 da Fase 4**.

## 2. Estado / progresso
- **Fase 6 COMPLETA e verificada.** Novo `runner.py` (RunSpec/plan/execute/resume/status) grava `runs/<run_id>/run.json` auto-contido (~67 KB) com `registry_snapshot` embutido + `registry_sha256`, `toolchain` (versão/hash de vina/obabel/boltz/python + GPU) e `tasks[]` (cmdline, hashes in/out, wall time, gates de 1ª classe). Gate reprovado = `gate_failed`, nunca exceção. Executor **injetável** (fronteira exata com a Fase 4). Merge não-destrutivo de `subprojeto/<ID>/status.json` (preserva as 11 etapas curadas do SUB01).
- Adicionado ao `registry.py`: `snapshot()`, `snapshot_sha256()` e CLI `--sha256`.
- **Fase 4 PARCIAL (#1 e #3 feitas).** Extraída `resolve_target_paths(args, base)` em `md_simulation.py`: alvo sem MD preparada e XML de force field ausente agora **ABORTAM** (sys.exit 2/3) em vez de cair para TcTR / reusar XML em silêncio; `--target` sem `choices=` fixo; reuso legítimo só via `--fad-xml`/`--lig-xml` explícitos.
- **Gates todos verdes:** `test_runner_provenance` 6/6; `test_md_target_safety` 5/5; `test_registry_parity` 5/5; `registry.py` válido; `audit_hardcode --check` = **439** (baseline 538, delta −99).
- Descoberta de ambiente: **sem vina.exe nesta máquina** (toolchain do run.json: `vina=False`) → docking real não roda aqui; verificação por executor injetável + `--dry-run`.

## 3. Arquivos criados / modificados
(no repo `C:\github\docking_suite_platform`, branch `feat/guarda-chuva-multialvo` — **NADA commitado**)
- `runner.py` — **NOVO**. Orquestração + proveniência (Fase 6).
- `registry.py` — +`snapshot()`/`snapshot_sha256()` + CLI `--sha256` (arquivo untracked, da Fase 1).
- `md_simulation.py` — **MODIFICADO**. `resolve_target_paths()` + correções #1/#3 (fallback de alvo e reuso de XML agora abortam).
- `tests/test_runner_provenance.py` — **NOVO**. 6 testes (executor fake, sem Vina).
- `tests/test_md_target_safety.py` — **NOVO**. 5 testes (módulo openmm fake, sem GPU).
- Memória PMM (OneDrive): `.project-memory/PROJECT_STATE.md` e `TASKS.md` atualizados (T-0011 DONE; T-0012 em andamento).

## 4. Decisões tomadas (e justificativa)
- **snapshot() no `registry.py`** (não no runner) — "o registry sabe se serializar"; aditivo, não entra no WATCHED do auditor nem quebra o teste de paridade.
- **Executor injetável no runner** — permite testar a proveniência sem binário de docking e é o ponto exato onde a Fase 4/engines entra sem reescrever o runner.
- **Merge não-destrutivo do `status.json`** — o SUB01 curou 11 etapas à mão (manuscrito, figuras); o runner atualiza só a etapa executada e nunca inventa subprojeto/etapa.
- **Fase 4: fazer #1 e #3, adiar #2** — #1/#3 matam degradação silenciosa e NÃO mudam ciência; #2 (TARGET_MAP→registry) tem armadilha científica (ver seção 5) e exige mudança de schema + verificação dedicada.

## 5. Pendências / próximos passos
- [ ] **Fase 4 #2**: TARGET_MAP (`md_simulation.py:45-49`) → `registry.get_target()`. **ARMADILHA**: `n_fads` é curado (TcTR=2, LiTR=2, **hGR=1**) e NÃO deriva de `keep_cofactors × min_protomers` (hGR é dímero, `min_protomers=2` no par, mas usa 1 FAD no arquivo 3DJJ_clean.pdb). Migrar o **VALOR** curado para o registry (campo novo, ex. `structure.n_cofactor_copies`), **NUNCA re-derivar** por fórmula.
- [ ] **Fase 4**: pacote `engines/` (base, vina, boltz2, openmm_md, haddock3) + `toolchain.py`, plugado no executor injetável do runner.
- [ ] Fases 5.2-5.4 (UI multi-campanha), 7 (SUB02 actina/Rho_OCRL), 8-9 (`export_web`, deleções).
- [ ] Decidir se commita o branch `feat/guarda-chuva-multialvo` (nada commitado ainda).
- [ ] **Recomendado: reiniciar a sessão** antes de emendar o restante da Fase 4 (regra de custo de contexto do autor: ~200k tokens).

## 6. Como retomar
- Gates a rodar sempre (do repo `C:\github\docking_suite_platform`, python do env `datascience`):
  `python tests/test_registry_parity.py`, `python tests/test_runner_provenance.py`, `python tests/test_md_target_safety.py`, `python registry.py`, `python tools/audit_hardcode.py --check`.
- Runner: `python runner.py plan --campaign SUB01` | `python runner.py run --campaign SUB01 --dry-run`.
- Ao reabrir (após `/compact`, `/clear` ou `/branch`), **releia** `docs/historico/CONTEXTO_SESSAO.md` e a memória PMM (`.project-memory/`), e continue pela seção 5.

## 7. Instruções de compactação (o que preservar)
- Preserve: objetivo (seção 1), decisões (seção 4) e pendências (seção 5), com destaque para a **ARMADILHA do n_fads**.
- Pode descartar: saídas de ferramentas antigas e trechos longos já resumidos aqui.
