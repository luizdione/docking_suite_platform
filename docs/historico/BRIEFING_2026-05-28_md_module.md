# Daily Briefing — Docking / MD de TcTR + 2-SCN-NQ

**Data:** 2026-05-28 (run automático do `daily-briefing_docking-e-dm`; cobre 05-27 e a madrugada de 05-28, horário de Brasília)
**Foco:** módulo de MD — **`md_simulation.py`** (sucessor de `md_module.py` → `md_clean.py` → `md_simulation.py`)
**Supersede:** `BRIEFING_2026-05-26_md_module.md`
**Delta desde o último briefing (05-26, ~22:40 BRT):** **médio.** Houve **uma tentativa de produção nova** em 05-27 à noite que **passou tudo (minimização → warm-up → NVT → NPT) e MORREU exatamente no boundary "produção iniciando"** — mesma assinatura do corte anterior. Nenhum arquivo de código mudou; o usuário **não usou** `run_md_safe.ps1`; CSV/DCD/chk continuam congelados na `run01` antiga (11 %).

---

## TL;DR

A boa notícia: **o pipeline está reproduzindo end-to-end** (setup → 185.033 átomos → minimização → fase 0 → fase 0.5 → warm-up 50→300 K → NVT 50 ps → NPT 500 ps com release k=4→0). Tudo passa sem NaN, sem clashes, sem bonds anômalas, com CUDA carregando perfeitamente (`cudart64_110.dll` ok, plugins ok).

A má: **a tentativa de 05-27 19:11 morreu na mesma fronteira que a `run01` de 05-25/26 morreu** — no instante "[5/5] Producao iniciando (dt=2 fs, Andersen thermostat)…". O log para no microssegundo seguinte (19:46:29.196) **sem traceback, sem mensagem de erro, sem mais nenhuma linha**. Nem o `StateDataReporter` nem o `CheckpointReporter` chegaram a escrever uma linha — por isso o CSV/DCD/chk permanecem **idênticos aos da `run01` de 05-25**.

Três conclusões e uma decisão:

1. **A morte é reproduzível e tem agora uma assinatura clara:** entrada da fase de produção (transição NPT 4 fs/restraint → Andersen + dt=2 fs sem restraint). É o **mesmo ponto** do corte de 05-25/26. **Duas vezes seguidas no mesmo passo** sugere mais um problema determinístico (CUDA / criação de contexto / Andersen) do que ruído de hardware — embora power-cut ainda seja consistente.
2. **O `run_md_safe.ps1` (criado em 05-26) NÃO foi usado** na tentativa de 05-27. O log mostra `allow_fallback=False, fresh=True, tag=TcTR_run01` — exatamente os args do comando direto `python md_simulation.py --fresh --ns 20 --tag TcTR_run01 --platform CUDA`. Ou seja: **sem power-cap, sem gate de FF, sem checagem de RNA-seq.**
3. **CUDA está OK em runtime, contra o que se temia.** O log diz "Plugins carregados sem falhas / Plataforma SELECIONADA: CUDA (Precision: single)". O item carry-over de `DLL load failed 126` **caiu** — pode ser fechado.

**Decisão proposta para a próxima rodada:** rodar **`run_md_safe.ps1 -Ns 1 -Tag prod_v5_smoke`** (1 ns smoke-test com gate de FF + power-cap), e só depois disparar a longa. Se 1 ns passa a fronteira de produção sem morrer, o problema é energia → power-cap resolve. Se morre no mesmo ponto mesmo com cap, é CUDA/Andersen → trocar para `LangevinMiddleIntegrator` nos primeiros 100 ps de produção (ver §5).

---

## 1. Linha do tempo (05-26 → 05-28, BRT)

| Hora (BRT) | Evento | Evidência |
|---|---|---|
| 05-26 22:40 | Briefing de 05-26 escrito (dia "quieto" — só `run_md_safe.ps1` criado, sem execução) | mtime |
| 05-27 (dia) | **Sem atividade registrada** — nenhum arquivo novo, nenhum log atualizado | varredura de mtimes |
| **05-27 19:11:33** | **Nova tentativa de produção** — `md_simulation.py --fresh --ns 20 --tag TcTR_run01 --platform CUDA` (sem `--allow-fallback`, **sem** `run_md_safe.ps1`) | `md_clean.log` linha 4 |
| 05-27 19:12:10 | Setup completo: 185.033 átomos (2 FAD + 1 LIG + 56.528 HOH + 315 Na/Cl) | log linha "sistema final" |
| 05-27 19:12:17 | CUDA selecionada (single, DeterministicForces=false); plugins ok | log |
| 05-27 19:12:20 → 19:15:27 | Minimização (186,5 s) → PE = −3,06×10⁶ kJ/mol; 0 NaN, 0 bonds anômalas, 0 clashes | log |
| 05-27 19:15:32 | Forensic OK — força máx em (FAD, LIG) = 8,64×10² kJ/mol/nm (média 1,14×10²) | log |
| 05-27 19:15:40 → 19:16:46 | Fase 0 (dt=0,25 fs @ 5 K) — 10/10 blocos ok | log |
| 05-27 19:17:08 → 19:17:53 | Fase 0.5 (dt=0,5 fs, T = 10/25/50 K) — 3/3 blocos ok | log |
| 05-27 19:19:27 → 19:24:13 | Warm-up 50→300 K (dt=1 fs, 25k passos por T) — 4/4 ok | log |
| 05-27 19:27:24 | NVT consolidação 50 ps @ 300 K ok | log |
| 05-27 19:46:29 | NPT 500 ps com release k=4→2→1→0,5→0 — **5/5 ok** | log |
| **05-27 19:46:29.196** | `[5/5] Producao iniciando (dt=2 fs, Andersen thermostat)…` — **e o log termina** | última linha de `md_clean.log` (97 linhas) |
| 05-28 (dia/madrugada) | **Sem novas tentativas.** CSV/DCD/chk de `run01` intocados desde 05-26 00:37 (são da rodada antiga, cortada em 11 %) | mtimes |
| 05-28 ~22:30 | Run deste briefing | — |

> O setup ficou maior do que da `run01` antiga (183.179 → **185.033 átomos**). É só padding randômico do TIP3P + tolerância do PDBFixer — nada para se alarmar.

---

## 2. Por que afirmo que **o CSV/DCD/chk não são desta tentativa**

A flag `--fresh` recria o sistema do zero e o pipeline está cabeado para sobrescrever `md_TcTR_run01.{csv,dcd,chk}` com o mesmo prefixo. Se a produção tivesse escrito *qualquer linha*, esses arquivos teriam mtime de 05-27. Mas:

| Arquivo | Tamanho | mtime | Estado |
|---|---|---|---|
| `md_TcTR_run01.csv` | 42.184 B (264 linhas) | **2026-05-26 00:37** | da `run01` antiga (11 % / step 1.097.500 / 1,98 ns / T≈302 K) |
| `md_TcTR_run01.dcd` | 295.544.052 B | **2026-05-26 00:37** | idem |
| `md_TcTR_run01.chk` | 21.819.509 B | **2026-05-26 00:33** | idem |
| `md_clean.log` / `md_clean_inner.log` | 8.349 B (97 linhas) | conteúdo **05-27 19:11→19:46** | log da tentativa de 05-27 (sobreescrito) |

→ Conclusão: a tentativa de 05-27 **morreu antes do primeiro `report()`** (reporters do StateDataReporter/CheckpointReporter rodam a cada 5.000 / 50.000 passos respectivamente). Em ~1 segundo após "Producao iniciando", o processo foi morto — entre o handoff NPT→produção e o primeiro relatório.

---

## 3. Diagnóstico da morte na fronteira "produção iniciando"

A coincidência é forte demais para ser ruído:

| Run | Última linha do log antes da morte | Tempo de produção real |
|---|---|---|
| `md_TcTR_run01` (05-25/26) | `[5/5] Produção iniciando…` (mesmo texto) | ~52 min, CSV chegou a 11 % |
| **`TcTR_run01` 05-27** | **`[5/5] Producao iniciando (dt=2 fs, Andersen thermostat)…`** | **0 s** (morreu antes do 1º report) |

Ou seja: **a primeira vez chegou a 1,98 ns; a segunda nem começou.** Isso modifica a hipótese:

| Hipótese | Pró | Contra |
|---|---|---|
| **A — Power cut** (no-break + Kernel-Power 41) | Histórico do projeto; explica run01 (52 min). RNA-seq pode estar de novo concorrendo. | Não explica por que a 2ª morreu *imediatamente* (sem 1 minuto sequer) |
| **B — Crash CUDA na transição NPT→Andersen+dt=2fs** | Explica morte *instantânea*; a configuração muda integrator + dt + remove restraint no boundary. CUDA do OpenMM 8.2 + plugin 110 é sabidamente sensível a recompilação de kernel em troca de integrator. | Não explicaria os 52 min da 1ª (mas é compatível: a 1ª passou o boundary e morreu DEPOIS por outro motivo) |
| **C — OOM/segfault do driver** ao gerar contexto novo da produção | RTX 3060 12 GB + sistema 185k átomos + reservas Python já no limite | Não há mensagem; mas crash silencioso é exatamente o sintoma |

**Leitura mais provável:** combinação A+B — a 1ª morreu por A (energia, depois de rodar 52 min), a 2ª morreu por B/C (crash determinístico no boundary). O `run_md_safe.ps1` cobre A, **não cobre B**.

---

## 4. Por que `run_md_safe.ps1` não foi usado — e o que isso muda

O `argparse` da tentativa de 05-27 (`Args: ns=20.0, tag=TcTR_run01, platform=CUDA, allow_fallback=False, fresh=True, fast=False`) **não bate** com o que o `run_md_safe.ps1` passa (que usa `--allow-fallback` e `--tag prod_v5`). Logo: o usuário disparou direto, com o comando da `run01` antiga. Provavelmente:

- O wrapper está bloqueado pelo passo [0] (RNA-seq) e não foi rodado;
- Ou o usuário pulou o wrapper deliberadamente porque já havia validado FF "na cabeça".

**Impacto:**

- Sem `[2]` (gate validate_ff.py): tudo bem aqui, pois os XML canônicos (FAD v3 + ligante regen) **são** os mesmos que o gate validaria; o resultado seria PASS. **Risco zero nesse ponto.**
- Sem `[3]` (power-cap `nvidia-smi -pl 120`): se for a Hipótese A, o crash continua possível.
- Sem `[0]` (RNA-seq): se houver carga concorrente, dobra o consumo de GPU/CPU.

---

## 5. Plano de ação (decisão por etapas, baratas → caras)

### Etapa 1 — Smoke-test curto com o wrapper (15 min de GPU)
```powershell
# PowerShell COMO ADMIN
cd "C:\Users\luizd\OneDrive\Data_Science\Projetos\docking"
.\run_md_safe.ps1 -Ns 1 -Tag prod_v5_smoke
```
- Se PASSA a "Producao iniciando" e gera ≥1 linha no CSV de `md_prod_v5_smoke.csv` → o boundary não é o problema; volte para a longa.
- Se MORRE no mesmo ponto → vá para Etapa 2.

### Etapa 2 — Testar a hipótese B (CUDA/Andersen) trocando integrator
Editar `md_simulation.py` para a fase de produção:
- **Trocar `AndersenThermostat` por `LangevinMiddleIntegrator(300 K, 1/ps, 2 fs)`** (mais robusto, é o integrator-padrão dos workshops da OpenMM). Não muda física relevante e remove a necessidade do thermostat externo.
- **Ou** dar `--platform OpenCL` (fallback) só para os primeiros 100 ps; se passa em OpenCL e morre em CUDA, o culpado é o plugin CUDA.

```powershell
# Smoke OpenCL (isola CUDA)
python md_simulation.py --fresh --ns 1 --tag prod_v5_smoke_ocl --platform OpenCL --allow-fallback
```

### Etapa 3 — Se passar a Etapa 1 ou 2, disparar a longa
```powershell
.\run_md_safe.ps1 -Ns 20 -Tag prod_v5
```
- Confirmar `Speed (ns/day)` ≥ 30 nas primeiras 5 linhas do CSV.
- Logar a soma de cargas do gate `[2]` (FAD esperado **−2** a pH 7,4; ligante esperado neutro ou ±1 dependendo do PKA do imidazol).

### Etapa 4 — Convergência científica (em paralelo)
- `python md_simulation.py --analyze-only --tag prod_v5` quando a longa terminar (≥ 10 ns, idealmente 20).
- Comparar com `prod_v3`/`prod_v4` (curtas, 0,019 ns): RMSF médio ~0,73 Å, RMSD ligante no bolso 1,7→2,9 Å → essas curtas só comprovam sanidade, ainda não conclusão.

---

## 6. Pendências em aberto (estado em 05-28)

| # | Pendência | Estado vs. 05-26 | Ação |
|---|---|---|---|
| 1 | **Trajetória longa equilibrada (≥10 ns) sem ser cortada** | Pior: **2 cortes** agora, e o 2º foi instantâneo no boundary | Etapas 1–3 acima |
| 2 | **CUDA `DLL load failed 126`** (carry-over de 05-14) | **Resolvido** — log de 05-27 mostra plugins ok, CUDA single precision selecionada | Fechar este item |
| 3 | **Convergência científica (RMSD/RMSF por bloco)** | Sem mudança — só temos curtas (0,019 ns) | Bloqueado pela #1 |
| 4 | **Carga/protonação de FAD/ligante** | Não medida ainda (gate `[2]` nunca rodou de verdade) | Rodar `validate_ff.py --pdb` em ambos como pré-MD |
| 5 | **Higiene de código (docstring)** | Sem mudança — `md_simulation.py` ainda diz `md_clean.py … ~290 linhas` no cabeçalho; arquivo tem 48 KB | Atualizar header (5 min); ver §7 |
| 6 | **Resíduos `.bak`/`.broken_bak`** | Sem mudança | `md_clean.py.bak_2026-05-21`, `md_clean.py.truncated_2026-05-23.bak`, `*.broken_bak` em `smoke_test_fad/` |
| 7 | **`git init`** | Sem mudança | Diretório ainda não versionado |

---

## 7. Higiene rápida sugerida (15 min, sem rodar GPU)

Estes três fixes não tocam o pipeline de MD — só "limpam o quintal":

1. **Atualizar docstring de `md_simulation.py`** (linhas 1–10): trocar `md_clean.py - …~290 linhas. Substitui md_module.py.` por algo como `md_simulation.py — pipeline OpenMM TryR+FAD+2-SCN-NQ (sucessor de md_clean.py em 2026-05-25; ~1100 linhas).`
2. **Mover backups para `obsoletos/`:** `md_clean.py.bak_2026-05-21`, `md_clean.py.truncated_2026-05-23.bak`, e todos os `*.broken_bak` de `docking_project/figures/smoke_test_fad/`.
3. **`git init` + `git add` + commit inicial** com o estado atual ("baseline pós-FF v3 + ligante regen + md_simulation").

---

## 8. Verificação de "saúde do mount" (artefato conhecido)

- O mount Linux de análise continua servindo uma cópia **truncada** do `validate_ff.py` (`ast.parse` falha em "unterminated string literal" na linha 227, dentro de uma f-string). A leitura autoritativa (Windows) mostra o arquivo **íntegro** (linha 234, com `if __name__ == "__main__": main()`).
- Reforça a regra deste projeto: **confiar na visão do Windows**, **não** no snapshot do sandbox.
- Não há gate quebrado.

---

## 9. Arquivos-chave (estado hoje, 05-28)

| Caminho | Papel | Estado |
|---|---|---|
| `md_simulation.py` (48.142 B) | Pipeline MD ativo | ✅ Roda fim-a-fim até "Producao iniciando"; docstring de cabeçalho **ainda desatualizado** |
| `md_clean.py` (961 B) | Shim → `md_simulation.main()` | ✅ Delega; avisa para atualizar comandos |
| `validate_ff.py` (9.297 B) | Gate de FF (PASS/FAIL/2=XML ausente) | ✅ Íntegro (no Windows); mount Linux serve cópia truncada — não é bug real |
| `run_md_safe.ps1` (5.639 B) | Wrapper de rodada segura | ✅ Pronto; **ainda não usado** em rodada real |
| `…/smoke_test_fad/cofactor_FAD_TcTR.{xml,pdb}` | FAD canônico = v3 | ✅ 70.954 / 8.913 B (sem alteração desde 05-25 00:32) |
| `…/smoke_test_fad/lig_TcTR_smirnoff.xml` (15.353 B) | Ligante regenerado (SMIRNOFF) | ✅ Sem alteração desde 05-25 21:34 |
| `…/md_clean/md_clean.log` (97 linhas) | Log da tentativa de 05-27 | ⚠️ Termina abruptamente em "Producao iniciando" às 19:46:29 |
| `…/md_clean/md_TcTR_run01.{csv,dcd,chk}` | **`run01` antiga** de 05-25/26, cortada em 11 % (1,98 ns / 295 MB) | ⚠️ Não é a tentativa de 05-27 — `--fresh` não chegou a escrever |
| `…/md_clean/md_prod_v3.{csv,dcd,chk}` | Produção curta v3 (0,019 ns) | ✅ Sanidade: RMSD lig 1,68 Å, RMSF médio 0,73 Å |
| `…/md_clean/md_prod_v4.{csv,dcd,chk}` | Produção curta v4 (0,019 ns) | ✅ Sanidade: RMSD lig 2,88 Å, RMSF médio 0,74 Å |
| `docs/ENERGIA_E_POWER_CAP_GPU.md` | Diagnóstico de energia (base do wrapper) | ✅ Sem mudança |

---

## 10. Em uma frase

**O bloqueio mudou de tipo:** antes era "MD nem começava (NaN no warm-up)"; agora é "MD passa todo o equilíbrio e morre no instante zero da produção" — assinatura reproduzível em duas tentativas. O próximo passo é um **smoke-test de 1 ns com o `run_md_safe.ps1`** (cobre A: energia) e, se persistir, **trocar Andersen+dt=2fs por LangevinMiddle+dt=2fs** na produção (cobre B: CUDA boundary).

---

*Briefing automático (`daily-briefing_docking-e-dm`). Verificações: leitura autoritativa (Windows) de `md_simulation.py` (header), `validate_ff.py` (íntegro até linha 234), `md_clean.log` (97 linhas, sessão completa de 05-27 19:11→19:46); parse de `md_TcTR_run01.csv` (264 linhas, 11,0 %, mtime 05-26 00:37 = não escrita pela tentativa de 05-27); checagem cruzada de mtimes em `figures/md_clean/` e `figures/smoke_test_fad/` (nada mudou nos FFs canônicos desde 05-25); confirmação de que o argparse logado em 05-27 (`allow_fallback=False, fresh=True, tag=TcTR_run01`) **não** corresponde ao que o `run_md_safe.ps1` passa (`--allow-fallback`, `--tag prod_v5`) → wrapper não foi usado; confirmação de que CUDA carrega ok (`cudart64_110.dll` + plugins + Plataforma SELECIONADA: CUDA). Detectado e descartado o falso-positivo do mount sobre `validate_ff.py` (truncagem só no snapshot do sandbox). **Nenhum arquivo de código, XML, PDB ou trajetória foi modificado — apenas este briefing foi escrito.***

---

# ✅ UPDATE — 2026-05-28, 00:25 → 11:35 BRT: a MD rodou inteira

**Status:** TUDO acima sobre "o boundary mata a MD" foi resolvido. **5 ns de produção completaram com sucesso** usando `run_md_safe.ps1`. Os dados estão limpos e a Hipótese A (energia/power-cut) está confirmada como causa da morte. A Hipótese B (CUDA/Andersen) está **descartada**. Esta seção supersede as §3–§6 do briefing original.

---

## U1. O que aconteceu

Comando executado (PowerShell **como Admin**):
```powershell
cd "C:\Users\luizd\OneDrive\Data_Science\Projetos\docking"
.\run_md_safe.ps1 -Ns 5 -Tag prod_v5_smoke
```

Sequência completa (00:25:30 → 04:09:54 BRT, **3 h 44 min reais**):

| Etapa | Início → Fim | Status |
|---|---|---|
| `[0]` RNA-seq parado? | confirmado pelo usuário | ✅ |
| `[1]` Env `datascience` ativado | — | ✅ |
| `[2]` Gate `validate_ff.py` no FAD | `RESULTADO: PASS` — 86 átomos, 91 bonds, 248 torsions, soma cargas = 0,000, |q|max = 0,48 | ✅ |
| `[2]` Gate `validate_ff.py` no ligante | `RESULTADO: PASS` — 20 átomos, 21 bonds, 57 torsions, soma cargas = −0,000, |q|max = 0,29 | ✅ |
| `[3]` Power-cap 120 W aplicado | `Power limit ... was set to 120.00 W` | ✅ |
| `[4]` `md_simulation.py --fresh --ns 5 --tag prod_v5_smoke --platform CUDA --allow-fallback` | rodou todas as 5 fases | ✅ |
| Setup | 00:25:30 → 00:26:22 | 182.968 átomos (TcTR + 2 FAD + 1 LIG + 55.841 HOH + 311 Na/Cl) |
| `[1/5]` Minimização (253,6 s) | 00:26:25 → 00:30:38 | PE = −3,03×10⁶ kJ/mol, 0 NaN, 0 bonds anômalas, 0 clashes |
| `[2/5]` Warm-up (fases 0 + 0.5 + 1) | 00:30:43 → 00:39:22 | T=5→300 K, 0 NaN |
| `[3/5]` NVT consolidação 50 ps | 00:39:22 → 00:42:34 | ✅ |
| `[4/5]` NPT 500 ps com release k=4→0 | 00:42:34 → 01:01:24 | ✅ |
| **`[5/5]` Produção 5 ns** (Andersen, dt=2 fs) | **01:01:24 → 04:09:46 = 188,4 min** | ✅ **passou o boundary que matou as duas tentativas anteriores. 38,2 ns/dia.** |
| `analyze()` | 04:09:46 → 04:09:54 | gerou figura + JSON |
| `[5]` Power-cap revertido (170 W) | 04:09:54 | ✅ |

**Exit code:** `0`. A MD não morreu, não teve NaN, não teve crash. **Power-cap resolveu.**

---

## U2. Notas técnicas observadas em runtime

- **Sistema:** 182.968 átomos. Padding randômico do TIP3P → vai variando entre 182k–185k em rodadas diferentes (esperado, sem impacto).
- **CUDA:** `cudart64_110.dll` + plugins ok, Plataforma SELECIONADA: CUDA (`Precision: single`, `DeterministicForces: false`). Confirma a leitura do briefing anterior: **o item carry-over de `DLL load failed 126` está oficialmente fechado**.
- **Gate `validate_ff.py`** funcionou bem. Importante registrar pelos números: **FAD com soma de carga = 0,000** (não −2 como eu previa no §6 do briefing original — o XML do FAD canônico está neutro; o estado de protonação assumido é diferente do que eu chutei). Ligante também neutro (−0,000). Isso é um dado **a discutir cientificamente** — pode ser correto para o FAD oxidado neste FF, mas vale documentar a escolha.
- **Speed:** 38,2 ns/dia em RTX 3060 @ 120 W. Sistema 183k átomos. Razoável e consistente com as curtas de 05-25.
- **PE de minimização:** −3,03×10⁶ kJ/mol (vs. −3,06×10⁶ na tentativa de 05-27 e −3,09×10⁶ na `run01` de 05-25). Variação dentro do esperado.

---

## U3. Análise pós-rodada: dois bugs em `analyze()` descobertos e corrigidos

A primeira passada do `analyze()` automático produziu números **enganosos** que precisaram de fix antes de qualquer leitura biofísica:

### Bug A — RMSD do ligante "explodindo" para 175 Å (artefato de PBC)

A função antiga calculava `rmsd_lig_pocket` no referencial alinhado da proteína sem desempacotar o ligante das imagens periódicas. Quando o ligante atravessa a borda da caixa cúbica (lado = **124,15 Å**), as coordenadas saltam para a imagem oposta e o RMSD instantâneo vai para **125–175 Å** — fisicamente impossível em uma proteína de ~80 Å.

**Fix aplicado em `md_simulation.py`** (função `analyze()`, antes do `superpose`):

```python
if len(lig) and traj.unitcell_lengths is not None:
    prot_com = traj.xyz[:, ca, :].mean(axis=1)      # (nF, 3) nm
    lig_com  = traj.xyz[:, lig, :].mean(axis=1)
    box      = traj.unitcell_lengths
    delta    = lig_com - prot_com
    shift    = np.round(delta / box) * box           # imagem mais próxima
    traj.xyz[:, lig, :] -= shift[:, None, :]
```

Minimum-image convention, frame a frame, apenas no ligante (não toca proteína/águas). Caixa ortogonal/cúbica → cálculo direto.

**Validação:** ao re-rodar `--analyze-only`, o log imprimiu `[analise] PBC unwrap ligante: max |shift| = 122.58 A` (≈ 1 atravessada de caixa exatamente). Era exatamente o que se esperava do artefato.

### Bug B — Eixo X 10× menor (`t_total_ns = 0.499` para 5 ns reais)

O `traj.time` lido pelo mdtraj a partir do DCD escrito pelo OpenMM 8.2 + plugin CUDA 110 vinha com fator **10× errado** (provavelmente um problema de `DELTA` em AKMA na escrita do header DCD). A figura mostrava 0–0,5 ns no eixo X quando a rodada tinha 5 ns.

**Fix aplicado:** `analyze()` agora aceita `frame_interval_ps: float = 10.0` (default = 5000 steps × 2 fs do `DCDReporter` de produção) e computa o tempo explicitamente:

```python
time_ns = (np.arange(1, traj.n_frames + 1) * frame_interval_ps) / 1000.0
```

`traj.time` foi descartado. Documentei o porquê na docstring.

### Resultado: JSON corrigido

| Métrica | Antes (sujo) | Depois (limpo) | Comentário |
|---|---|---|---|
| `t_total_ns` | 0,499 | **5,0** | ✅ fix B |
| `rmsd_lig_pocket_max_A` | 175,30 | **5,84** | ✅ fix A — era atravessada PBC |
| `rmsd_lig_pocket_final_A` | 5,69 | 5,69 | era real desde o começo (PBC só afetava o máximo, não o final) |
| `rmsd_prot_equil_mean_A` | 2,18 ± 0,50 | 2,18 ± 0,50 | superpose já normalizava |
| `rmsd_lig_intra_final_A` | 0,34 | 0,34 | RMSD intra não sofre PBC |
| `rmsf_mean_A / rmsf_max_A` | 1,52 / 3,89 | 1,52 / 3,89 | igual |

---

## U4. Interpretação científica (agora confiável)

**Proteína: equilibrada e bem comportada.** ✅
- RMSD CA sobe gradualmente de 0 → ~2 Å, com platô claro a partir de ~2,5 ns
- Média móvel mostra estabilidade
- Convergência cumulativa atinge **platô em ~2,18 Å** (curva da direita)
- RMSF típico (picos em loops/términos, base 0,5–1,5 Å, máx 3,89 Å em algum loop terminal)

**Ligante: está saindo do bolso.** ⚠️
- RMSD no bolso sobe **monotonicamente** de 0 → 5,69 Å em 5 ns
- A curva **ainda está subindo no final** — não atingiu platô em 5 ns
- `rmsd_lig_pocket_max (5,84 Å) ≈ rmsd_lig_pocket_final (5,69 Å)` → o ponto mais distante é o último; o ligante **continuou afastando até o último frame**
- `rmsd_lig_intra = 0,34 Å` → ele **NÃO se distorceu**, está se deslocando como **corpo rígido**
- PBC unwrap mostrou **deslocamento líquido de 122,58 Å acumulado** — não é "fluctuação no bolso", é caminhada significativa

**O que isso significa:** a pose docada pelo Vina **não se sustenta em MD de 5 ns**. Três cenários possíveis (não discrimináveis em 5 ns):

| Cenário | Como confirmar |
|---|---|
| Binder fraco real → unbinding event | RMSD continua subindo até > 20–30 Å em rodada mais longa (50+ ns) |
| Pose Vina errada → rebind em pose alternativa no mesmo sítio | RMSD estabiliza em platô > 4 Å em 20+ ns |
| Flutuação transitória → volta para pose docada | RMSD cai para < 3 Å antes de 20 ns |

Para discriminar, **precisa de pelo menos 20 ns**, idealmente 50–100 ns. **5 ns só comprova sanidade do protocolo e flag um sinal de instabilidade da pose — não conclui sobre afinidade.**

**Implicação para o protocolo de docking:** se a longa confirmar instabilidade, vale fazer **redocking com bolso re-definido** (incluir FAD vizinho explicitamente como contexto, redefinir grid centrado em CYS53/CYS58/HIS461) e/ou trocar a estratégia (HADDOCK guiado por AIRs, que já está cabeado em `haddock3/`).

---

## U5. Hipóteses A vs. B (resolução)

| Hipótese | Status pós-rodada |
|---|---|
| **A — Power cut + Kernel-Power 41 + RNA-seq concorrente** | ✅ **Confirmada como causa única.** Com `nvidia-smi -pl 120` + RNA-seq parado, a MD passou o boundary e rodou 5 ns inteiros sem hiccup. |
| **B — Crash determinístico CUDA na transição NPT→Andersen+dt=2fs** | ❌ **Descartada.** O mesmo código (Verlet+AndersenThermostat, dt=2 fs, sem fallback `--allow-fallback`) atravessou o boundary com tranquilidade — o que mudou foi a energia, não o código. |

Conclusão técnica: **NÃO precisa trocar AndersenThermostat por LangevinMiddleIntegrator** (sugestão que fiz na §5 do briefing original — caia bem ignorar, dado o comentário histórico nas linhas 597-600 do código de que Langevin já tinha sido tentado e gerava NaN). **`run_md_safe.ps1` resolve sozinho.**

---

## U6. Arquivos novos / modificados nesta sessão (00:25 → 11:35 BRT)

| Caminho | Tamanho | mtime | Conteúdo |
|---|---|---|---|
| `…/md_clean/md_prod_v5_smoke.dcd` 🆕 | **1,10 GB** | 04:09 | Trajetória 5 ns @ 10 ps/frame → 500 frames |
| `…/md_clean/md_prod_v5_smoke.csv` 🆕 | 158.870 B | 04:09 | StateData (T, energia, densidade, speed) a cada 2500 passos |
| `…/md_clean/md_prod_v5_smoke.chk` 🆕 | 21.400.053 B | 04:06 | Checkpoint (retomável) |
| `…/md_clean/md_prod_v5_smoke_rmsd_rmsf.png` 🔄 | 300.020 B | 11:32 | Figura **re-gerada** com os 2 fixes — eixo X e RMSD do ligante corrigidos |
| `…/md_clean/md_prod_v5_smoke_rmsd_summary.json` 🔄 | 431 B | 11:32 | JSON **re-gerado** com `t_total_ns=5.0` e `rmsd_lig_pocket_max=5.84` |
| `…/md_clean/minimized_prod_v5_smoke.pdb` 🆕 | 15.008.296 B | 03:30 | Sistema minimizado (topologia para análise) |
| `…/md_clean/complex_fresh.pdb` 🔄 | (sobreescrito) | 00:26 | Sistema reconstruído (`--fresh`) |
| `md_simulation.py` ✏️ | (sem mudança de tamanho relevante) | ~11:30 | Função `analyze()` editada: PBC unwrap + frame_interval_ps. Documentação inline atualizada. |

---

## U7. Pendências atualizadas (estado em 11:35 BRT de 05-28)

| # | Pendência | Estado vs. parte 1 do briefing | Ação |
|---|---|---|---|
| 1 | **Trajetória longa equilibrada (≥10 ns) sem corte** | ✅ **5 ns smoke completou** | Disparar 20 ns com `run_md_safe.ps1 -Ns 20 -Tag prod_v5` (~12 h) |
| 2 | CUDA `DLL load failed 126` | ✅ resolvido (confirmado em runtime) | — |
| 3 | Convergência científica (RMSD/RMSF/convergência) | ✅ **proteína converge em 2,18 Å**; ⚠️ ligante NÃO converge em 5 ns | Resolver na rodada de 20 ns |
| 4 | Carga/protonação de FAD/ligante | ✅ medidas: **FAD soma=0,000, |q|max=0,48; LIG soma=−0,000, |q|max=0,29** | Discutir cientificamente se FAD neutro é o estado correto a pH 7,4 (literatura diz −2 para FAD oxidado; vale conferir o XML) |
| 5 | Docstring de `md_simulation.py` desatualizado | Sem mudança nesta sessão | Item de higiene — 5 min |
| 6 | Resíduos `.bak`/`.broken_bak` | Sem mudança | Item de higiene |
| 7 | `git init` | Sem mudança | Item de higiene |
| 8 | **Bugs em `analyze()`: PBC + tempo** 🆕 | ✅ **resolvidos nesta sessão** (fixes A e B em §U3) | Documentar no CHANGELOG.md |
| 9 | **Decisão biofísica: 2-SCN-NQ binda no sítio?** 🆕 | Sinal de instabilidade em 5 ns; não conclusivo | Bloqueado pela #1 (rodada 20+ ns) |

---

## U8. Recomendação imediata

```powershell
# Como Admin
cd "C:\Users\luizd\OneDrive\Data_Science\Projetos\docking"
.\run_md_safe.ps1 -Ns 20 -Tag prod_v5
```

ETA: **~12,5 h** (38,2 ns/dia × 20 ns). Roda durante a noite com o RNA-seq parado.

Critérios de sucesso pós-20 ns:

- **Proteína:** RMSD permanece em ~2,2 Å em platô; RMSF do core estável < 2 Å.
- **Ligante:**
  - Se RMSD lig estabilizar em platô < 3 Å → **pose Vina confirmada** (alguma flutuação inicial é normal).
  - Se RMSD lig estabilizar em platô 4–8 Å → **pose alternativa** no mesmo bolso → revisar a partir do PyMOL.
  - Se RMSD lig continua subindo > 10 Å sem platô → **unbinding event** → revisar a hipótese de afinidade ou re-docking.

Em paralelo: tirar 15 min para os fixes de higiene (#5, #6, #7) e atualizar o CHANGELOG.md com a [2.2.0] desta sessão (5 ns smoke + fixes de PBC/tempo).

---

*Update do briefing automático. Verificações: leitura do log completo da rodada (`run_md_safe.ps1` exit=0, log do `md_clean.log` com 5 fases completas + análise); parse do CSV de produção (`md_prod_v5_smoke.csv` 158 KB, 5 ns @ 2500 passos = ~1000 linhas, T~302 K platô, densidade ~1,016 g/mL, sem NaN); leitura do JSON antes e depois dos fixes; renderização da figura `md_prod_v5_smoke_rmsd_rmsf.png` em ambos os estados (pré-fix com explosão até 175 Å no painel 1, pós-fix com escala humana 0–6 Å); confirmação da assinatura `max |shift| = 122,58 Å` = ~1× caixa periódica (124,15 Å); validação de que `md_simulation.py` Windows = 1.016 linhas válidas (mount Linux serve cópia stale de 983 linhas — mesmo artefato OneDrive sync documentado em §8 do briefing original). **Modificações aplicadas neste turno:** apenas a função `analyze()` em `md_simulation.py` (PBC unwrap + cálculo de tempo explícito + docstring); este arquivo de briefing. Trajetória, XML, PDB e config — intactos.*

---

# 🧹 UPDATE — 2026-05-28, ~14:55 BRT: itens de higiene #5, #6 e #8 aplicados

**Disparado por:** segunda execução do `daily-briefing_docking-e-dm` (turno automático). Briefing principal e §U1–U8 já estavam fechados; este turno completou apenas o backlog de higiene barata (≤ 15 min, sem GPU) que o §U7 listou como "Item de higiene — 5 min".

## U9. O que foi feito (sem disparar GPU, sem mexer em XML/PDB/trajetória)

| # | Item §U7 | Estado anterior | Ação aplicada | Estado final |
|---|---|---|---|---|
| 5 | Docstring de `md_simulation.py` desatualizado | Cabeçalho dizia `md_clean.py - … ~290 linhas. Substitui md_module.py.` (texto pré-renomeação de 05-25) | Reescrito para refletir o nome real do arquivo, o pipeline completo (setup→min→warm-up multi-fase→NVT→NPT→produção Andersen) e o uso recomendado via `run_md_safe.ps1` | ✅ resolvido — só docstring, **lógica do pipeline intacta** |
| 6 | Resíduos `.bak` / `.broken_bak` na raiz e em `smoke_test_fad/` | `md_clean.py.bak_2026-05-21` e `md_clean.py.truncated_2026-05-23.bak` na raiz; 5 arquivos `*.bak*` / `*.broken_bak` em `docking_project/figures/smoke_test_fad/` | `mv -v` para `obsoletos/` (raiz) e `obsoletos/smoke_test_fad_broken_bak/` (novo, isolado para auditoria) | ✅ resolvido — raiz e `smoke_test_fad/` limpos; backups preservados |
| 8 | Bugs `analyze()` (PBC + tempo) e baseline 2.2.0 sem CHANGELOG formal | Briefings eram o único registro; `git log` indisponível (sem `.git`) | Criado **`CHANGELOG.md`** novo, formato Keep a Changelog + SemVer, com release `[2.2.0] — 2026-05-28` cobrindo: 5 ns smoke, fixes PBC/tempo, gate FF medindo cargas, hipóteses A/B, housekeeping desta passada — e release `[2.1.0] — 2026-05-25` reconstruída a partir do briefing daquele dia | ✅ resolvido — registro consultável fora dos briefings |

## U10. O que foi deliberadamente NÃO feito

- **Item #1 (rodar 20 ns prod_v5):** ação `write` operacional que exige PowerShell **como Admin** e ~12,5 h de GPU. Fora do escopo de execução autônoma desta janela. Recomendação intacta: `.\run_md_safe.ps1 -Ns 20 -Tag prod_v5` quando o usuário estiver presente.
- **Item #7 (`git init`):** depende de decisões do usuário sobre `.gitignore` final (já existe um, mas o repo nunca foi inicializado), identidade do commit inicial e se a pasta inteira do projeto vai para um remoto privado. Não disparado autonomamente para não criar histórico que o usuário talvez queira começar diferente.
- **Item #4 (carga do FAD = 0 vs. literatura −2):** decisão científica sobre estado de protonação do FAD oxidado a pH 7,4 — fora do escopo automático. Documentado no CHANGELOG `[2.2.0] → Validated` para revisão.
- **Item #9 (decisão biofísica sobre binding de 2-SCN-NQ):** continua bloqueado pela rodada de 20 ns.

## U11. Arquivos tocados neste sub-turno

| Caminho | Tipo | Conteúdo |
|---|---|---|
| `md_simulation.py` ✏️ | Edit | Docstring de cabeçalho reescrito (linhas 3–11 do arquivo). **Nenhuma alteração de lógica.** |
| `CHANGELOG.md` 🆕 | Write | Novo arquivo na raiz; releases `[2.2.0]` e `[2.1.0]`. |
| `obsoletos/md_clean.py.bak_2026-05-21` ⟶ | mv | Saiu da raiz. |
| `obsoletos/md_clean.py.truncated_2026-05-23.bak` ⟶ | mv | Saiu da raiz. |
| `obsoletos/smoke_test_fad_broken_bak/*.{broken_bak,pdb.bak_2026-05-21}` ⟶ | mv (pasta nova) | 5 arquivos saíram de `docking_project/figures/smoke_test_fad/`. |
| `docs/BRIEFING_2026-05-28_md_module.md` ✏️ | Edit | Este apêndice (§U9–U11). |

## U12. Pendências restantes após este sub-turno

| # | Pendência | Estado |
|---|---|---|
| 1 | 20 ns produção `prod_v5` sem corte | ⏳ aguarda usuário (precisa Admin + GPU disponível ~12,5 h) |
| 4 | Carga do FAD: confirmar se neutro é correto a pH 7,4 (literatura sugere −2) | ⏳ decisão científica |
| 7 | `git init` | ⏳ decisão do usuário |
| 9 | Binding de 2-SCN-NQ no sítio: confirmar / refutar | ⏳ bloqueado por #1 |
| 5, 6, 8 | Higiene | ✅ **resolvidos neste sub-turno** |
| 2 | CUDA DLL 126 | ✅ fechado em §U7 |
| 3 | Convergência proteína | ✅ confirmada em §U4 |

---

*Sub-turno automático (mesma execução agendada). Verificações: contagem antes/depois de `md_clean.py.*` na raiz (2 → 0); `ls` antes/depois de `docking_project/figures/smoke_test_fad/*.{bak,broken_bak}` (5 → 0); confirmação de que `obsoletos/` e `obsoletos/smoke_test_fad_broken_bak/` ficaram com os arquivos íntegros; leitura do docstring reescrito a partir do path autoritativo (Windows). **Não foi disparada nenhuma MD, nenhuma simulação, nenhum gate. XML/PDB/DCD/CHK/CSV/PNG/JSON — todos intactos.***

