# Daily Briefing — Docking / MD de TcTR + 2-SCN-NQ

**Data:** 2026-05-25 (run automático do `daily-briefing_docking-e-dm`; cobre a sessão de 05-25 e a madrugada de 05-26)
**Foco:** módulo de MD — agora **`md_simulation.py`** (sucessor de `md_module.py` → `md_clean.py` → `md_simulation.py`)
**Supersede:** `BRIEFING_2026-05-23_md_module_UPDATE.md`
**Delta desde o último briefing (05-23/24):** **enorme — o maior salto do projeto**

---

## TL;DR

**A MD DESTRAVOU.** O bloqueio que durava desde 14/05 (`Particle coordinate is NaN` no warm-up) **acabou**. Hoje rodaram **três produções reais** que atravessaram minimização → warm-up (5→300 K) → NVT → NPT → produção **sem um único NaN**, com a proteína e o ligante estáveis. Não é mais "a um comando de distância": é MD rodando e gerando trajetória.

O que tornou isso possível foram **DUAS causas-raiz da mesma natureza**, descobertas e corrigidas:

1. **FAD (cofator)** — force-field parametrizado para a identidade atômica errada + torsions vazias (já diagnosticado em 05-21/23). O par **v3** foi promovido a canônico hoje (00:32) — os arquivos antigos viraram `.broken_bak`.
2. **Ligante 2-SCN-NQ (NOVO achado de hoje)** — o SDF/PDB usado na MD tinha **todas as 37 ligações como ordem 1**. O RDKit `AddHs` então **super-protonou** o esqueleto `C11NO2S` para 21 H (`C11H21NO2S`) — uma naftoquinona saturada quimicamente impossível. O SMIRNOFF gerava um FF **sem rigidez de anel** → a MD rodava sobre a **molécula errada**. Corrigido por `regen_ligand_xml.py` (separa POSE de QUÍMICA, reatribui ordens de ligação, reprotonação correta, SMIRNOFF).

**Conclusão conceitual:** os dois resíduos custom (cofator E ligante) sofriam do **mesmo defeito** — química perdida → FF sem torsions → física errada / NaN. Os dois foram regenerados pela **mesma rota robusta** (`regen_fad_xml.py` / `regen_ligand_xml.py`). Esse é o avanço da semana.

**Novo gargalo (mudou de software para hardware):** rodadas longas estão **desligando o PC** (Kernel-Power Event 41; no-break apitando). É **energia**, não bug — agravado por **RNA-seq rodando em paralelo** com a MD. Mitigação documentada: não sobrepor cargas + power-cap da GPU (`nvidia-smi -pl 120`).

> Também houve **refatoração**: `md_clean.py` agora é só um **shim** que delega para **`md_simulation.py`** (o módulo real). Comandos antigos seguem funcionando.

---

## 1. Linha do tempo (2026-05-25 → madrugada 05-26)

| Hora (aprox.) | Evento | Evidência |
|---|---|---|
| 05-25 00:32 | **Par FAD v3 promovido a canônico** (`cofactor_FAD_TcTR.xml/.pdb` = v3; antigos → `.broken_bak`) | mtime 00:32; tamanhos 70954/8913 == v3 |
| 05-25 01:08 | `diag_warmup_nan.py` criado — localiza passo-a-passo o **1º átomo que explode** no warm-up (nome/resíduo/cadeia) | novo arquivo |
| 05-25 ~04:17–04:56 | **1ª produção bem-sucedida** (`md_prod_v3`): minimizado + DCD (44 MB) + CSV | mtimes; 0 NaN no CSV |
| 05-25 04:30 | `check_ligand_binding.py` — pré-voo: confirma que o LIG está **ancorado no bolso** antes de gastar GPU | novo arquivo |
| 05-25 05:08 | `analyze_dcd_numpy.py` — análise de trajetória sem MDAnalysis (numpy puro) | novo arquivo |
| 05-25 ~16:05 | Diagrama do pipeline (`pipeline_md_clean.png/svg`) | mtime |
| 05-25 21:24 | `md_prod_v3` analisado (RMSD/RMSF `.png` + `_summary.json`) | mtime |
| 05-25 21:30 | **`md_clean.py` vira shim**; `regen_ligand_xml.py` criado (fix do ligante); `regen_fad_xml.py` tocado; `.bat` novos | mtimes |
| 05-25 21:34 | **Ligante regenerado** (`lig_TcTR_smirnoff.xml` + `lig_TcTR.pdb`, e cópias `_v2`) | mtime; 15353 bytes |
| 05-25 21:44 | `dashboard_v2.py` atualizado | mtime |
| 05-25 ~21:54–22:33 | **2ª produção** (`md_prod_v4`): minimizado + DCD + CSV + análise | mtimes; 0 NaN |
| 05-25 20:11 → 05-26 00:37 | **`md_TcTR_run01`** — produção de **20 ns** iniciada; log mostra pipeline completo OK; DCD chega a **295 MB** | `md_clean.log`; CSV 265 linhas |
| 05-26 02:08–02:10 | `validate_ff.py` criado (portão de FF); **`md_simulation.py`** finalizado (rename/refator) | mtimes |

---

## 2. Evidência de que a MD está saudável (não é mais NaN)

Do `md_clean.log` (run `md_TcTR_run01`, lido pela visão autoritativa do Windows):

```
[diag/post-min] atoms with NaN positions: 0   |   NaN forces: 0
[diag/post-min] PotentialEnergy = -3.09e+06 kJ/mol
[forensic] resumo: 0 bonds anomalas, 0 clashes <1.0 A     <-- o "BOND ANOMALA 2.06 A N3A-C4A" de 14/05 SUMIU
[fase0] bloco 1..10/10 (dt=0.25fs @ 5K) ok                <-- onde antes explodia, agora passa
[fase0.5] T=10/25/50 K ok
warm-up T=50/100/200/300 K ok
[3/5] NVT consolidacao 50 ps @ 300 K ok
[4/5] NPT 500 ps @ 1 bar (release gradual restraint k=4->0) ok
[5/5] Producao 20.0 ns (dt=2 fs, 10,000,000 steps)...
```

Os CSVs confirmam estabilidade em produção (sem NaN em nenhuma linha):

| Run | Linhas CSV | Tempo coberto | Temperatura | Energia total | Speed | NaN? |
|---|---|---|---|---|---|---|
| `md_prod_v3` | 41 | ~0,2 ns (curta) | 303–305 K | ~−1,99×10⁶ kJ/mol | ~38 ns/dia | **0** |
| `md_prod_v4` | 41 | ~0,2 ns (curta) | ~302–304 K | ~−1,99×10⁶ kJ/mol | ~38 ns/dia | **0** |
| `md_TcTR_run01` | 265 | ~2 ns de 20 ns (11%) | ~302–303 K | ~−2,01×10⁶ kJ/mol | ~37 ns/dia | **0** |

Temperatura travada em ~300 K, energia total plana, densidade ~1,01 g/mL. É exatamente o perfil de uma MD **bem comportada**.

---

## 3. Qualidade estrutural das rodadas (RMSD / RMSF)

Dos `_summary.json` (medidos com o ligante **no bolso** — superpõe a proteína e mede a deriva do LIG, capturando saída do sítio):

| Métrica (Å) | `prod_v3` | `prod_v4` | Leitura |
|---|---|---|---|
| RMSD proteína (final) | 1,66 | 1,78 | estável (típico < 2 Å) |
| RMSD proteína (equilíbrio, média±dp) | 1,54 ± 0,10 | 1,36 ± 0,22 | platô consistente |
| RMSD ligante no bolso (final / máx) | 1,68 / 3,51 | 2,88 / 2,88 | permanece no sítio, com deriva limitada |
| RMSD ligante intramolecular (final) | 1,02 | 0,38 | geometria interna do LIG **estável** |
| RMSF médio / máx | 0,73 / 1,76 | 0,74 / 2,53 | flexibilidade normal (loops nas pontas) |

> ⚠️ **Caveat de escala.** As análises de `prod_v3/v4` cobrem só ~0,02 ns (20 frames) — são **smoke/curtas**, ótimas para sanidade, **insuficientes para conclusão científica**. A produção de verdade é a `run01` (alvo 20 ns). Os números acima provam que a física está sã, **não** que há convergência. A convergência só será aferida na trajetória longa (curva de RMSD/RMSF + bloco).

---

## 4. Novos guard-rails criados hoje (para nunca mais perder semanas)

| Ferramenta | O que faz | Por que importa |
|---|---|---|
| `validate_ff.py` | **Portão de FF em segundos**: torsions não-vazias, bonds k>0, soma de cargas inteira, σ-LJ coerente com o elemento, elemento×nome, e **smoke-test de `createSystem`** | Pega *antes* da MD os exatos bugs que custaram semanas (NaN críptico só aparecia depois de minutos de GPU) |
| `diag_warmup_nan.py` | Roda o warm-up **passo a passo** e aponta o **1º átomo NaN/maior força** (nome/resíduo/cadeia) | Trocou "adivinhar a origem do NaN" por "medir" — foi a chave para fechar o caso |
| `check_ligand_binding.py` | Pré-voo: o LIG está **ancorado no bolso**? (dist mín, resíduos em contato, dist ao FAD; veredito) | Evita gastar horas de GPU numa pose dessorvida |
| `regen_ligand_xml.py` | Regenera FF+template de **qualquer ligante** separando POSE (docking) de QUÍMICA (SMILES/SDF/mol2/CCD) | Generaliza o conserto do FAD para ligantes — **reutilizável para compostos novos** |
| `md_simulation.py` | Renomeação/refator do pipeline (`md_clean.py` → shim) | Nome do módulo agora condiz com a função; comandos antigos preservados |

---

## 5. O novo gargalo: ENERGIA / hardware (não é software)

Documentado em `docs/ENERGIA_E_POWER_CAP_GPU.md` (05-25):

- **Sintoma:** PC **desliga abruptamente** durante a MD; no-break **apitando**. Visualizador de Eventos: `Kernel-Power ID 41` recorrente (06/05, 15–17/05, 24/05).
- **Causa:** pico GPU+CPU (~350–400 W) **excede o no-break** (ou bateria degradada) → proteção corta a energia. **Sem** WHEA/BugCheck → não é crash de CPU/RAM/driver. A MD em si está sã.
- **Gatilho provável:** **RNA-seq rodando em paralelo** (CPU 100%) **junto** com a MD (GPU 100%) soma os dois maiores consumos.
- **Mitigação nº 1 (de graça):** **não rodar RNA-seq e MD ao mesmo tempo** — serialize.
- **Paliativo:** `nvidia-smi -pl 120` (RTX 3060; cai de ~170 W). Custo em velocidade ~10–25%; **uma rodada capada que termina >> uma rodada a 170 W cortada**.
- **Durável (produção rotineira):** trocar bateria do no-break / PSU ≥650 W, **ou** rodar a produção (ns) em **GPU na nuvem** (sistema já montado e validado; leva-se só os arquivos).

> Provável reflexo disso: a `run01` (20 ns) tem CSV/DCD/checkpoint parando em ~2 ns (11%) às 00:37. Não dá para confirmar daqui se ela **terminou** ou foi **interrompida** (corte de energia é a hipótese mais coerente com o histórico). **Existe `md_TcTR_run01.chk`** → é **retomável**.

---

## 6. Pendências em aberto

1. **Confirmar o destino da `run01` (20 ns).** Ver se completou ou caiu por energia; se caiu, **retomar do checkpoint** (`md_TcTR_run01.chk`) — de preferência **sem RNA-seq concorrente** e com `nvidia-smi -pl 120`.
2. **Convergência científica.** As curtas (v3/v4) só provam sanidade. Falta a trajetória longa equilibrada + curva de convergência (RMSD/RMSF por bloco) para qualquer conclusão de binding.
3. **CUDA do OpenMM (carry-over).** Plugins falham com `DLL load failed 126`; cai para OpenCL. Independe do FAD/ligante; `conda install -c conda-forge openmm-cuda` e conferir `cudart64_*.dll`.
4. **Carga/protonação do FAD e do ligante.** Confirmar no log da `regen` a fórmula/carga total (FAD oxidado a pH 7,4 costuma ser **−2**) — rodar `validate_ff.py` nos dois XML (FAD_v3 + lig) como portão padrão antes de produção.
5. **Higiene de código/repo:**
   - `md_simulation.py` ainda tem o **docstring antigo** ("`md_clean.py` … ~290 linhas") — atualizar o cabeçalho (o arquivo tem ~48 KB).
   - `git init` na pasta (continua **não** sendo repositório) para versionar de verdade.
   - resíduos para limpar: `md_clean.py.truncated_2026-05-23.bak`, `cofactor_FAD_TcTR_v3_elemfix.xml(.INVALID.txt)`, `*.broken_bak`.

---

## 7. Próximos passos priorizados (env `datascience`, Windows)

1. **Antes de qualquer MD, rodar o portão de FF:**
   ```powershell
   conda activate datascience
   cd "C:\Users\luizd\OneDrive\Data_Science\Projetos\docking"
   python validate_ff.py docking_project\figures\smoke_test_fad\cofactor_FAD_TcTR_v3.xml --pdb docking_project\figures\smoke_test_fad\cofactor_FAD_TcTR_v3.pdb
   python validate_ff.py docking_project\figures\smoke_test_fad\lig_TcTR_smirnoff.xml
   # esperar "RESULTADO: PASS" nos dois
   ```
2. **Serializar cargas** (encerrar RNA-seq) e **capar a GPU**, depois retomar/rodar a produção:
   ```powershell
   nvidia-smi -pl 120
   python md_simulation.py --fresh --ns 20 --tag prod_v5 --platform CUDA --allow-fallback
   ```
   Monitorar a coluna `Speed (ns/day)` do CSV e temperatura (HWiNFO). Checkpoint a cada 50k passos → retomável.
3. **Pré-voo do bolso** em complexos novos: `python check_ligand_binding.py docking_project\figures\md_clean\complex_fresh.pdb`.
4. Concluída a longa: `python md_simulation.py --analyze-only --tag prod_v5` (RMSD-no-bolso + RMSF + convergência).
5. **Em paralelo:** consertar CUDA (item 6.3) e a higiene de código/repo (item 6.5).

> ⚠️ Mantido de `NEXT_SESSION.md`: não reiniciar ciclos automáticos de MD sem o portão `validate_ff.py` PASS e sem resolver a energia — senão volta a cair (agora por hardware, não por NaN).

---

## 8. Arquivos-chave (estado hoje)

| Caminho | Papel | Estado |
|---|---|---|
| `md_simulation.py` (~48 KB) | **Pipeline MD ativo** (renomeado de `md_clean.py`) | ✅ roda fim-a-fim; docstring de cabeçalho **desatualizado** |
| `md_clean.py` (961 B) | Shim de compatibilidade → `md_simulation.main()` | ✅ delega; avisa para atualizar comandos |
| `regen_fad_xml.py` | Regenera o FF do FAD (NAGL) | ✅ usado; par v3 promovido a canônico |
| `regen_ligand_xml.py` 🆕 | Regenera FF+template de ligante (POSE×QUÍMICA) | ✅ corrigiu o ligante super-protonado |
| `validate_ff.py` 🆕 | Portão de validação de FF (segundos) | ✅ usar antes de toda MD |
| `diag_warmup_nan.py` 🆕 | Localiza o 1º átomo NaN no warm-up | ✅ ferramenta de diagnóstico |
| `check_ligand_binding.py` 🆕 | Pré-voo: LIG ancorado no bolso? | ✅ pré-MD |
| `…/smoke_test_fad/cofactor_FAD_TcTR.xml/.pdb` | Par FAD **canônico = v3** | ✅ 70954/8913 B; antigos em `.broken_bak` |
| `…/smoke_test_fad/lig_TcTR_smirnoff.xml` + `lig_TcTR.pdb` | Ligante regenerado | ✅ 05-25 21:34 (com `_v2`) |
| `…/md_clean/md_TcTR_run01.dcd/.csv/.chk` | Produção 20 ns (em curso/retomável) | ✅ DCD 295 MB; CSV ~2 ns; **chk existe** |
| `…/md_clean/md_prod_v3/v4.*` | Produções curtas + análise RMSD/RMSF | ✅ 0 NaN; só sanidade (escala curta) |
| `…/md_clean/md_clean.log` | Log da última rodada | ✅ pipeline completo OK |
| `docs/ENERGIA_E_POWER_CAP_GPU.md` 🆕 | Diagnóstico de energia + power-cap | ✅ ler antes de produção longa |
| `md_run_v2_prod.log` | Último run com **NaN** (14/05) | histórico — **superado** |

---

*Briefing automático (`daily-briefing_docking-e-dm`). Verificações: leitura autoritativa via Windows do `md_clean.log`, `md_simulation.py`, `validate_ff.py`, `regen_ligand_xml.py`, `diag_warmup_nan.py`, `check_ligand_binding.py`; parse numérico dos CSVs de produção (v3/v4/run01 → 0 NaN, T≈300 K, energia plana); leitura dos `_summary.json` (RMSD/RMSF); varredura de mtimes para a linha do tempo. **Nenhum arquivo de código, XML, PDB ou trajetória foi modificado — apenas este briefing foi escrito.** O destino exato da `run01` (terminou × cortada por energia) não é determinável daqui; rodar/retomar conforme §7 resolve.*
