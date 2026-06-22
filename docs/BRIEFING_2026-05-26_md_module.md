# Daily Briefing — Docking / MD de TcTR + 2-SCN-NQ

**Data:** 2026-05-26 (run automático do `daily-briefing_docking-e-dm`; cobre o dia 05-26, horário de Brasília)
**Foco:** módulo de MD — **`md_simulation.py`** (sucessor de `md_module.py` → `md_clean.py` → `md_simulation.py`)
**Supersede:** `BRIEFING_2026-05-25_md_module.md`
**Delta desde o último briefing (05-25, ~23:59 BRT):** **pequeno — dia "quieto".** Um único artefato novo de código (`run_md_safe.ps1`); **nenhuma MD nova rodou.**

---

## TL;DR

Depois do salto de ontem (a MD destravou: 0 NaN, FAD v3 e ligante regenerados), **hoje foi um dia de preparação, não de execução.** O único trabalho registrado é a criação do **`run_md_safe.ps1`** (01:03 BRT) — um wrapper que automatiza exatamente as recomendações de energia + portão de FF do briefing de ontem. **Nenhuma produção nova foi disparada** (não há saída `prod_v5`), e a `run01` (20 ns) continua parada em ~2 ns.

Três conclusões para hoje:

1. **A `run01` (20 ns) NÃO terminou — foi cortada.** Resolve a pendência nº 1 de ontem. O log para em "Produção iniciando" sem marcador de conclusão; o CSV chega a **11,0 % (≈1,98 ns de 20)** e estima **~11,5 h restantes**. A física estava sã o tempo todo (T≈302 K, energia plana, 0 NaN) — foi **interrupção de hardware** (consistente com o corte de energia / Kernel-Power 41), não falha do modelo. **`md_TcTR_run01.chk` existe → retomável.**
2. **O pipeline está "armado e pronto", mas não foi executado.** O `run_md_safe.ps1` está corretamente cabeado (verifiquei: todas as flags que ele passa existem no `md_simulation.py`; os 3 arquivos de FF que ele exige existem; o `validate_ff.py` compila). Só falta **rodá-lo** — provavelmente está esperando o RNA-seq terminar (o próprio passo [0] do script trava nisso).
3. **Decisão a tomar antes de rodar:** o script usa `--fresh` (cria uma `prod_v5` do zero) em vez de **retomar a `run01` do checkpoint**. As duas rotas são válidas — ver §4.

> **Nota de verificação importante:** o mount Linux de análise serviu uma cópia **truncada** do `validate_ff.py` (cortada em 9297 B); a leitura autoritativa (Windows) mostra o arquivo **completo e válido**. Reforça a regra deste projeto: confiar na visão do Windows, **não** no snapshot do sandbox. Não há gate quebrado.

---

## 1. Linha do tempo (05-26, BRT)

| Hora (BRT) | Evento | Evidência |
|---|---|---|
| 05-25 23:59 | Briefing de 05-25 escrito (fecha o "dia do salto") | mtime do `.md` |
| 05-26 01:03 | **`run_md_safe.ps1` criado** — wrapper de rodada segura (energia + gate de FF) | mtime; 121 linhas |
| 05-26 (dia) | **Nenhuma MD executada** — sem novos `prod_v5.*`; última saída segue sendo a `run01` (00:37 UTC = 21:37 BRT 05-25) | varredura de mtimes em `figures/md_clean/` |
| 05-26 ~22:40 | Run deste briefing | — |

Em outras palavras: o dia foi de **engenharia de processo** (deixar a próxima produção segura e reprodutível), não de **gerar trajetória**.

---

## 2. O artefato novo: `run_md_safe.ps1`

Wrapper PowerShell que materializa o §7 do briefing de ontem (serializar cargas + portão de FF + power-cap). Sequência:

| Passo | O que faz | Comportamento |
|---|---|---|
| **[0]** | Pergunta se o **RNA-seq já terminou** (serialização das cargas) | Aborta (exit 3) se não; `-SkipConfirm` pula |
| **[1]** | Prepara o env `datascience` (python + DLLs CUDA no PATH, igual ao dashboard) | Aborta (exit 2) se não achar `python.exe` |
| **[2]** | **PORTÃO de FF:** `validate_ff.py` no **FAD** (`--pdb`) e no **ligante** | **Aborta a MD (exit 1)** se qualquer um reprovar |
| **[3]** | **Power-cap** da GPU (`nvidia-smi -pl 120`, configurável) | Avisa se não for admin; pergunta se segue sem cap |
| **[4]** | `md_simulation.py --fresh --ns 20 --tag prod_v5 --platform CUDA --allow-fallback` | A produção em si |
| **[5]** | Reverte o cap para 170 W (padrão RTX 3060) | — |

**Parâmetros:** `-Ns` (default 20), `-Tag` (`prod_v5`), `-PowerLimit` (120), `-Platform` (CUDA), `-SkipConfirm`.
**Uso:** `.\run_md_safe.ps1` (20 ns); `.\run_md_safe.ps1 -Ns 5 -Tag teste5` (sanidade curta); `-PowerLimit 100` (cap mais agressivo).

### Verificações que fiz neste script (tudo confere)

- ✅ **Flags compatíveis:** `--fresh`, `--ns`, `--tag`, `--platform`, `--allow-fallback`, `--analyze-only` existem todas no `argparse` do `md_simulation.py`.
- ✅ **Arquivos de FF do gate existem:** `cofactor_FAD_TcTR.xml` (70.954 B = v3), `cofactor_FAD_TcTR.pdb` (8.913 B), `lig_TcTR_smirnoff.xml` (15.353 B).
- ✅ **`validate_ff.py` é válido** (compila; sai 0=PASS / 1=FAIL / 2=XML ausente — exatamente o contrato que o `[2]` consome).
- ⚠️ **Requer Administrador** para o `nvidia-smi -pl` funcionar; sem admin ele apenas avisa e segue (sem cap → maior risco de corte). Rodar o PowerShell **como admin**.

**Veredito:** o script está pronto e correto. O gargalo agora é **operacional** (rodá-lo com o RNA-seq parado), não de software.

---

## 3. Destino da `run01` (20 ns): **cortada, não concluída** (fecha a pendência nº 1)

Evidência cruzada (CSV + log, leitura autoritativa):

| Sinal | Valor | Leitura |
|---|---|---|
| Última linha do CSV | **11,0 %**, step **1.097.500 / 10.000.000**, t≈**1.983 ps (~1,98 ns)** | só ~10 % dos 20 ns-alvo |
| `Time Remaining` (última linha) | **~11 h 34 min** | faltava a maior parte |
| Estado físico até o corte | T≈302–303 K, energia ≈ −2,01×10⁶ kJ/mol plana, **0 NaN** | MD **saudável** até parar |
| Log (`md_clean.log`) | termina em `[5/5] Produção iniciando … 20.0 ns (10.000.000 steps)…` | **sem** marcador de conclusão |
| Janela de produção | ~20:45 → 21:37 BRT (≈52 min) | parou cedo |
| `md_TcTR_run01.chk` | existe (≈21,8 MB) | **retomável** |

→ Conclusão: **interrupção externa** (corte de energia é a hipótese coerente com o histórico de Kernel-Power 41 e com o no-break apitando). A MD em si não falhou. É exatamente o cenário que o `run_md_safe.ps1` foi feito para evitar (power-cap + RNA-seq serializado).

---

## 4. Decisão a tomar: `--fresh` (prod_v5) **vs** retomar a `run01`

O `run_md_safe.ps1` usa **`--fresh`** → começa uma **`prod_v5` do zero** e **descarta** o checkpoint da `run01` (equilibração + ~2 ns já calculados).

| Rota | Prós | Contras |
|---|---|---|
| **`--fresh` (prod_v5)** — como o script está | Rodada limpa **sobre os FFs já validados** (FAD v3 + ligante regen); mais fácil de defender cientificamente; sem dúvida sobre estado herdado | Joga fora a equilibração + ~2 ns da `run01` (~1 h de GPU) |
| **Retomar `run01.chk`** | Aproveita ~1 h já gasta; chega antes aos 20 ns | Exige flag de retomada (o `--fresh` é o contrário); herda o estado de uma rodada cortada |

**Recomendação:** se a equilibração da `run01` já usou os FFs canônicos atuais (FAD v3 + ligante regen — o que parece ser o caso), **`--fresh` prod_v5 é a escolha mais limpa** e foi a que o script já assume. Só vale retomar a `run01` se o tempo de GPU estiver muito caro. **Escolha consciente antes de disparar.**

---

## 5. Pendências em aberto (estado em 05-26)

1. **Rodar a produção longa.** Tudo está pronto: encerrar o RNA-seq → `.\run_md_safe.ps1` **como Administrador**. É o único passo que falta para voltar a gerar trajetória. *(Era a nº 1 de ontem; hoje a `run01` foi confirmada como cortada — a ação agora é simplesmente executar o script já pronto.)*
2. **Convergência científica.** As curtas (v3/v4, ~0,02 ns) só provam sanidade. Continua faltando a trajetória longa equilibrada + curva de convergência (RMSD/RMSF por bloco) para **qualquer** conclusão de binding.
3. **CUDA do OpenMM (carry-over).** Plugins falham com `DLL load failed 126` → cai para OpenCL. Independe de FAD/ligante. `conda install -c conda-forge openmm-cuda` e conferir `cudart64_*.dll`. (O `run_md_safe.ps1` já põe `Library\bin` no PATH, o que pode ajudar — confirmar.)
4. **Carga/protonação de FAD e ligante.** Confirmar no output do `validate_ff.py` a soma de cargas (FAD oxidado a pH 7,4 costuma ser **−2**). Agora **está embutido** no gate `[2]` do `run_md_safe.ps1` — basta ler o que ele imprime na primeira rodada.
5. **Higiene de código/repo (inalterado):**
   - `md_simulation.py` **ainda** tem o docstring antigo ("`md_clean.py` … ~290 linhas / Substitui md_module.py") — o arquivo tem ~48 KB. Atualizar o cabeçalho.
   - `git init` na pasta (segue **não** versionada).
   - Resíduos para limpar: `md_clean.py.truncated_2026-05-23.bak`, `md_clean.py.bak_2026-05-21`, `*.broken_bak`.

---

## 6. Próximos passos priorizados (env `datascience`, Windows, PowerShell **como Admin**)

1. **Confirmar RNA-seq parado** e disparar a produção segura:
   ```powershell
   cd "C:\Users\luizd\OneDrive\Data_Science\Projetos\docking"
   .\run_md_safe.ps1                      # 20 ns, tag prod_v5, cap 120 W
   # (ou para sanidade rápida primeiro:)
   .\run_md_safe.ps1 -Ns 5 -Tag teste5
   ```
   - Conferir que o gate `[2]` imprime **RESULTADO: PASS** nos dois XML (e anotar a soma de cargas).
   - Acompanhar `Speed (ns/day)` e T no CSV / dashboard. Checkpoint a cada 50k passos → retomável se cortar de novo.
2. **Decidir** explicitamente `--fresh prod_v5` (recomendado) **vs** retomar `run01.chk` (§4) **antes** de rodar.
3. Concluída a longa: `python md_simulation.py --analyze-only --tag prod_v5` (RMSD-no-bolso + RMSF + convergência por bloco).
4. **Em paralelo (não bloqueia a MD):** consertar CUDA (item 5.3) e a higiene de código/repo (item 5.5).

> ⚠️ Mantido: não reiniciar ciclos de MD sem o gate `validate_ff.py` = PASS e sem o RNA-seq serializado/power-cap — senão volta a cair (agora por hardware, não por NaN). O `run_md_safe.ps1` impõe isso automaticamente; **use-o** em vez de chamar `md_simulation.py` na mão.

---

## 7. Arquivos-chave (estado hoje)

| Caminho | Papel | Estado |
|---|---|---|
| `run_md_safe.ps1` 🆕 | Wrapper de rodada segura (RNA-seq off + gate FF + power-cap + MD) | ✅ pronto e cabeado; **ainda não executado**; rodar como Admin |
| `md_simulation.py` (~48 KB) | Pipeline MD ativo | ✅ roda fim-a-fim; docstring de cabeçalho **ainda desatualizado** |
| `md_clean.py` (961 B) | Shim → `md_simulation.main()` | ✅ delega; avisa para atualizar comandos |
| `validate_ff.py` | Portão de FF (segundos) | ✅ **completo e válido** (a "truncagem" vista no mount era artefato de sync) |
| `regen_fad_xml.py` / `regen_ligand_xml.py` | Regeneram FF de FAD / ligante | ✅ FAD v3 e ligante regen são os canônicos |
| `…/smoke_test_fad/cofactor_FAD_TcTR.{xml,pdb}` | Par FAD **canônico = v3** | ✅ 70.954 / 8.913 B |
| `…/smoke_test_fad/lig_TcTR_smirnoff.xml` | Ligante regenerado | ✅ 15.353 B |
| `…/md_clean/md_TcTR_run01.{dcd,csv,chk}` | Produção 20 ns **cortada em ~2 ns** | ⚠️ incompleta (11 %); **chk existe** → retomável |
| `…/md_clean/md_prod_v3 / v4.*` | Produções curtas + análise | ✅ 0 NaN; só sanidade (escala curta) |
| `…/md_clean/md_clean.log` | Log da `run01` | termina em "Produção iniciando" (sem conclusão) |
| `docs/ENERGIA_E_POWER_CAP_GPU.md` | Diagnóstico de energia + power-cap | ✅ base do `run_md_safe.ps1` |

---

*Briefing automático (`daily-briefing_docking-e-dm`). Verificações: leitura autoritativa (Windows) de `run_md_safe.ps1`, `md_simulation.py` (header), `validate_ff.py` (completo/válido); confirmação de que as flags do wrapper existem no `argparse` e de que os 3 XML/PDB do gate existem (tamanhos conferidos); parse do `md_TcTR_run01.csv` (11,0 %, step 1.097.500, 0 NaN) e do `md_clean.log` (sem marcador de conclusão); varredura de mtimes (nenhuma MD nova em 05-26). Detectado e descartado um falso-positivo: o mount de análise servia uma cópia truncada do `validate_ff.py` — o arquivo real está íntegro. **Nenhum arquivo de código, XML, PDB ou trajetória foi modificado — apenas este briefing foi escrito.***
