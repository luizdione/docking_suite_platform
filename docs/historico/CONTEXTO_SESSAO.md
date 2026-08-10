# Contexto da Sessão — Fase 6 + Fase 4 completas; arquitetura/UX da plataforma

- **Data/hora:** 2026-08-08 01:13 -03:00
- **Branch:** feat/guarda-chuva-multialvo
- **Repositório:** docking_suite_platform (C:\github\docking_suite_platform)
- **Foco desta sessão:** Concluir Fase 6 (runner/proveniência) e Fase 4 inteira (correções de md_simulation + pacote engines/); mapear a arquitetura de UX (localhost único, py3Dmol, pipeline E2E, CLI).

## 1. Objetivo / tarefa atual
Generalizar a plataforma de docking/MD para **guarda-chuva multi-alvo**, sem destruir o SUB01. Plano de 10 fases em `C:\Users\luizd\.claude\plans\analisando-a-plataforma-http-localhost-5-serene-steele.md`. Esta sessão fechou as Fases 6 e 4; próxima é a Fase 5 (UI multi-campanha) — que absorve os requisitos de UX abaixo.

## 2. Estado / progresso
- **Fase 6 COMPLETA**: `runner.py` (RunSpec/plan/execute/resume/status) + `runs/<id>/run.json` (registry_snapshot+sha256, toolchain, tasks[] com gates de 1a classe); `registry.snapshot()/snapshot_sha256()`+CLI `--sha256`; merge não-destrutivo do `status.json`. `test_runner_provenance` 6/6.
- **Fase 4 COMPLETA**: (a) 3 correções de `md_simulation.py` (#1 fallback de alvo, #2 `TARGET_MAP`→registry via `structure.cofactor_copies` com hGR=1 curado, #3 reuso de XML) — `test_md_target_safety` 9/9; (b) pacote **`engines/`** (base/toolchain/vina/boltz2/openmm_md/haddock3) plugado no runner (`_engine_executor`) — `test_engines` 10/10.
- **Gates todos verdes**; `audit_hardcode --check` = **439** (delta −99).
- **Arquitetura atual (verificada no código):** SIM, é **um único localhost** — `dashboard_v2.py` (Flask, http://localhost:5001) com ~40 rotas faz TUDO local: docking, MD, ADMET, toxicidade, drug design, ciclodextrina, figuras, biblioteca, geração de texto do paper. UI = `static/app.html` (extraída na Fase 5.1). O viewer no Hostinger é SEPARADO, estático e read-only (só o paper).
- **Visualização 3D no navegador: NÃO existe ainda** (grep confirmou: sem 3Dmol/py3Dmol/molstar/NGL no app.html). Requisito novo → T-0019.
- **Pipeline E2E:** o `runner.py` é o motor E2E, mas AINDA NÃO está conectado ao dashboard (as rotas hoje chamam docking_v2/md_simulation direto). Falta plugar o runner no dash.
- **CLI amigável:** `runner.py` já tem `plan/run/status/resume`; `docking_v2`/`md_simulation` têm CLIs próprios.
- **Toolchain nesta máquina:** vina.exe no OneDrive `Projetos\docking\`; boltz no env `foldgpu`; obabel (OB 3.1.0) no PATH; **haddock3 ausente**; GPU ok.

## 3. Arquivos criados / modificados (branch feat/guarda-chuva-multialvo, NADA commitado)
- `runner.py` — NOVO (Fase 6) + integração com engines (`_engine_executor`, `toolchain_snapshot` delega a `engines.toolchain`).
- `engines/` — NOVO: `__init__.py`, `base.py`, `toolchain.py`, `vina.py`, `boltz2.py`, `openmm_md.py`, `haddock3.py`.
- `registry.py` — +`snapshot()`/`snapshot_sha256()`+CLI `--sha256`.
- `registry/targets/{tctr,litr,hgr}.json` — +`structure.cofactor_copies` (n_fads curado: 2/2/1).
- `md_simulation.py` — `resolve_target_paths()`/`_target_md_info()`/`_md_ready_targets()`; `TARGET_MAP` removido; `import registry`.
- `tests/` — NOVO: `test_runner_provenance.py` (6), `test_md_target_safety.py` (9), `test_engines.py` (10).
- `claude.md` — +seção "Instruções de Compactação".
- Memória PMM (OneDrive): `.project-memory/PROJECT_STATE.md` e `TASKS.md` atualizados (T-0011/T-0012 DONE; T-0013 c/ requisitos UX; T-0019 novo).

## 4. Decisões tomadas (e justificativa)
- **Executor injetável no runner** = fronteira exata onde os engines entram; testes de proveniência (fake executor) intactos após a integração.
- **`n_fads` migrado como valor curado** (não derivado): hGR é dímero mas usa 1 FAD no 3DJJ; fórmula daria 2 e mudaria a MD.
- **engines: imports pesados tardios**; `dry_run` planeja o cmdline; toolchain ausente → erro claro (nunca degradar em silêncio).
- **HADDOCK3 gera .cfg+ambig.tbl** de protocolo+sites (lógica pura, testável sem o binário).
- **Requisitos de UX registrados** (T-0013/T-0019), NÃO implementados nesta sessão: py3Dmol/3Dmol.js no dash; botão para disparar o runner e mostrar progresso/gates; revisar amigabilidade.

## 5. Pendências / próximos passos
- [ ] **Fase 5.2-5.4** (T-0013): `/api/config?campaign=` servindo o registry ao front; migrar o JS uma aba por commit; seletor de campanha. IC50/toxicidade do 2-SCN-NQ por último (paper submetido).
- [ ] **UX/3D** (T-0019): integrar 3Dmol.js/py3Dmol no dashboard (receptor+pose inline); botão "rodar pipeline" chamando `runner.execute` com progresso do `run.json`.
- [ ] **Fase 7** (SUB02 actina/Rho_OCRL): 7.1 ingestão read-only dos 13 JSONs; 7.2 execução nativa.
- [ ] **Fase 8-9**: `verify_paper_vs_pipeline.py`, `export_web.py` → `SUB01.paper.json` frozen; deleções.
- [ ] Decidir se **commita** o branch `feat/guarda-chuva-multialvo` (nada commitado ainda).

## 6. Como retomar
- Gates (do repo `C:\github\docking_suite_platform`, python do env `datascience`):
  `python tests/test_registry_parity.py`, `python tests/test_runner_provenance.py`, `python tests/test_md_target_safety.py`, `python tests/test_engines.py`, `python registry.py`, `python tools/audit_hardcode.py --check`.
- Dashboard: `python dashboard_v2.py` → http://localhost:5001. Runner: `python runner.py plan --campaign SUB01`.
- Ao reabrir (após /compact, /clear, /branch), **releia este arquivo** e a memória PMM (`.project-memory/`), continue pela seção 5.

## 7. Instruções de compactação (o que preservar)
- Preserve: objetivo (1), estado da arquitetura (2), decisões (4) e pendências (5).
- Pode descartar: saídas de ferramentas antigas e trechos longos já resumidos aqui.
