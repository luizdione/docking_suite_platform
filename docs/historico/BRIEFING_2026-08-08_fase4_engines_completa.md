# Briefing — Fase 6 + Fase 4 completas; arquitetura/UX da plataforma

- **Data/hora:** 2026-08-08 01:13 -03:00
- **Branch:** feat/guarda-chuva-multialvo
- **Repositório:** docking_suite_platform (C:\github\docking_suite_platform)
- **Snapshot de:** sessão que concluiu Fase 6 (runner/proveniência) e Fase 4 inteira (correções de md_simulation + pacote engines/), e mapeou a arquitetura de UX.

## 1. Objetivo / tarefa atual
Generalizar a plataforma de docking/MD para **guarda-chuva multi-alvo**, sem destruir o SUB01. Plano de 10 fases em `C:\Users\luizd\.claude\plans\analisando-a-plataforma-http-localhost-5-serene-steele.md`. Fases 6 e 4 fechadas; próxima é a Fase 5 (UI multi-campanha), que absorve os requisitos de UX abaixo.

## 2. Estado / progresso
- **Fase 6 COMPLETA**: `runner.py` (RunSpec/plan/execute/resume/status) + `runs/<id>/run.json` (registry_snapshot+sha256, toolchain, tasks[] com gates); `registry.snapshot()/snapshot_sha256()`+CLI `--sha256`; merge não-destrutivo do `status.json`. `test_runner_provenance` 6/6.
- **Fase 4 COMPLETA**: (a) 3 correções de `md_simulation.py` (#1 fallback de alvo, #2 `TARGET_MAP`→registry via `structure.cofactor_copies` com hGR=1 curado, #3 reuso de XML) — `test_md_target_safety` 9/9; (b) pacote **`engines/`** (base/toolchain/vina/boltz2/openmm_md/haddock3) plugado no runner (`_engine_executor`) — `test_engines` 10/10.
- **Gates todos verdes**; `audit_hardcode --check` = **439** (delta −99).
- **Arquitetura atual (verificada no código):** SIM, é **um único localhost** — `dashboard_v2.py` (Flask, http://localhost:5001) com ~40 rotas faz TUDO local: docking, MD, ADMET, toxicidade, drug design, ciclodextrina, figuras, biblioteca, texto do paper. UI = `static/app.html`. Viewer Hostinger é SEPARADO, estático, read-only (só paper).
- **Visualização 3D no navegador: NÃO existe ainda** (sem 3Dmol/py3Dmol/molstar/NGL no app.html). Requisito novo → T-0019.
- **Pipeline E2E:** `runner.py` é o motor E2E, mas AINDA NÃO conectado ao dashboard (rotas hoje chamam docking_v2/md_simulation direto). Falta plugar.
- **CLI:** `runner.py` tem `plan/run/status/resume`; `docking_v2`/`md_simulation` têm CLIs próprios.
- **Toolchain nesta máquina:** vina.exe no OneDrive `Projetos\docking\`; boltz no env `foldgpu`; obabel (OB 3.1.0) no PATH; **haddock3 ausente**; GPU ok.

## 3. Arquivos criados / modificados (NADA commitado)
- `runner.py` (NOVO+integração engines); `engines/` (NOVO: __init__/base/toolchain/vina/boltz2/openmm_md/haddock3).
- `registry.py` (+snapshot/sha256/CLI); `registry/targets/{tctr,litr,hgr}.json` (+cofactor_copies 2/2/1).
- `md_simulation.py` (resolve_target_paths/_target_md_info/_md_ready_targets; TARGET_MAP removido).
- `tests/` (test_runner_provenance 6, test_md_target_safety 9, test_engines 10); `claude.md` (+compactação).
- Memória PMM: `.project-memory/PROJECT_STATE.md` e `TASKS.md` (T-0011/T-0012 DONE; T-0013 c/ UX; T-0019 novo).

## 4. Decisões tomadas (e justificativa)
- Executor injetável no runner = fronteira dos engines; testes de proveniência intactos.
- `n_fads` migrado como valor curado (hGR=1), não derivado (fórmula daria 2 e mudaria a MD).
- engines: imports pesados tardios; dry_run planeja; toolchain ausente → erro claro.
- HADDOCK3 gera .cfg+ambig.tbl de protocolo+sites (lógica pura, testável).
- Requisitos de UX (py3Dmol, botão E2E, amigabilidade) REGISTRADOS (T-0013/T-0019), não implementados nesta sessão.

## 5. Pendências / próximos passos
- [ ] Fase 5.2-5.4 (T-0013): `/api/config?campaign=`; JS data-driven uma aba por commit; seletor de campanha.
- [ ] UX/3D (T-0019): 3Dmol.js/py3Dmol no dash (receptor+pose); botão "rodar pipeline" via `runner.execute` com progresso do `run.json`.
- [ ] Fase 7 (SUB02 actina): 7.1 ingestão read-only; 7.2 execução nativa.
- [ ] Fase 8-9: `verify_paper_vs_pipeline.py`; `export_web.py` → `SUB01.paper.json` frozen; deleções.
- [ ] Decidir se commita o branch (nada commitado).

## 6. Como retomar
- Gates (repo `C:\github\docking_suite_platform`, python env `datascience`): test_registry_parity, test_runner_provenance, test_md_target_safety, test_engines, `registry.py`, `audit_hardcode --check`.
- Dashboard: `python dashboard_v2.py` → http://localhost:5001. Runner: `python runner.py plan --campaign SUB01`.
- Ao reabrir, releia `docs/historico/CONTEXTO_SESSAO.md` e `.project-memory/`, continue pela seção 5.

## 7. Instruções de compactação (o que preservar)
- Preserve: objetivo (1), arquitetura (2), decisões (4), pendências (5).
- Pode descartar: saídas de ferramentas antigas e trechos longos já resumidos aqui.
