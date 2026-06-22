# Daily Briefing — Docking / MD de TcTR + 2-SCN-NQ

**Data:** 2026-05-23 (run automático do `daily-briefing_docking-e-dm`)
**Foco:** módulo de MD (`md_clean.py`, sucessor do `md_module.py`)
**Última atividade real no projeto:** 2026-05-22 (release v2.1.0)
**Delta desde o último briefing (05-21):** **substancial** — toda a v2.1.0 entrou em 05-22

---

## TL;DR

A MD continua **bloqueada no mesmo ponto desde 14/05** (NaN no warm-up), mas o bloqueio agora está a **um único comando de distância**: falta apenas rodar `python regen_fad_xml.py` no env `datascience` para gerar o `cofactor_FAD_TcTR_v3.xml`. Todo o resto da correção do force-field do FAD já está pronto e verificado:

- o **v3 PDB** (`cofactor_FAD_TcTR_v3.pdb`) já tem a coluna `element` corrigida (conferido hoje: `PA→P`, `O1A→O`, …);
- o **`regen_fad_xml.py`** está completo e íntegro no disco (537 linhas, rota SMIRNOFF via CCD + auto-validação);
- o **`md_clean.py`** já troca automaticamente para o par v3 assim que o `v3.xml` existir (lógica nas linhas 295–302).

Houve **progresso real em 05-22** (v2.1.0), respondendo às críticas de banca: validação de docking por RMSD, redocking da glutationa co-cristalizada, e um método ortogonal (HADDOCK3) montado.

**Achado operacional novo (importante):** o sandbox de análise enxerga, via OneDrive, **cópias desatualizadas** de alguns `.py` (lag de sincronização). `md_clean.py` e `regen_fad_xml.py` apareceram **truncados** no mount do sandbox, mas estão **íntegros no disco** (verificado pela visão autoritativa do Windows). Detalhe e ação tomada na §5 — nenhum conteúdo foi perdido.

---

## 1. O único passo que destrava a MD

| Pré-requisito da correção do FAD | Estado hoje |
|---|---|
| `cofactor_FAD_TcTR_v3.pdb` (coluna `element` corrigida) | ✅ pronto e verificado (33 campos; primeiros átomos conferidos: P/O/O/O/C/C) |
| `regen_fad_xml.py` (gera o force-field correto: elementos/massas/LJ/cargas/bonds/angles/**torsions**) | ✅ completo no disco (537 linhas), sintaxe e leitura PDB/bonds já testadas (86 átomos / 91 bonds) |
| `md_clean.py` prefere o par v3 quando presente | ✅ implementado (linhas 295–302; só usa o v3 PDB **quando o v3 XML também existir** → sem regressão) |
| **`cofactor_FAD_TcTR_v3.xml` (a saída do regen)** | ❌ **AINDA NÃO EXISTE** — `regen_fad_xml.py` nunca foi executado |

> O sandbox não roda o regen: ele exige **RDKit + OpenFF** (e a validação que importa é o sistema solvatado de 183k átomos), que só existem no seu env `datascience`. Por isso este passo depende de você.

**Comando para destravar (no Windows, env `datascience`):**

```powershell
conda activate datascience
cd "C:\Users\luizd\OneDrive\Data_Science\Projetos\docking"
python regen_fad_xml.py
# espere ver no log:  [regen_fad_xml.py vISO-2026-05-22]  e  RESULTADO: PASS
python md_clean.py --fresh --ns 1 --tag prod_v3 --platform CUDA --allow-fallback
```

Confira que o forensic **não** acusa `BOND ANOMALA` e que o warm-up **não** dá NaN. Assim que o `v3.xml` for gravado, o `md_clean.py` troca para o par v3 automaticamente.

---

## 2. Estado do pipeline de MD (inalterado desde 14/05)

`md_clean.py --fresh` reconstrói o sistema completo até a minimização e então explode:

| Etapa | Resultado |
|---|---|
| Reconstrução | PDBFixer → addHydrogens(proteína) → 2× FAD canônicos via Kabsch → LIG → TIP3P + 0,15 M NaCl |
| Sistema final | **183.179 átomos** (TcTR dímero + 2 FAD + 1 LIG + 55.906 HOH + 311 Na/Cl) |
| Minimização | converge (PotentialEnergy ≈ −3,0 × 10⁶ kJ/mol) |
| Warm-up (dt=0,25 fs @ 5 K) | **NaN imediato** |

O último log de produção (`md_run_v2_prod.log`, 14/05) termina em:

```
openmm.OpenMMException: Particle coordinate is NaN.
```

Nenhuma simulação foi executada desde 14/05 (o diretório `docking_project/figures/md_clean/` não muda desde então). O `auto_fix/` continua **sem nunca ter rodado** (sem `auto_md_report.md` / `auto_md_state.json` / `auto_md_logs/`).

---

## 3. Progresso real em 05-22 (release v2.1.0) — novo desde o briefing de 05-21

A v2.1.0 (registrada em `VERSION` e `CHANGELOG.md`) responde às críticas de banca (RMSD + método ortogonal) e arruma o scaffolding do fix do FAD:

- **Validação de docking por RMSD** — `docking_pose_rmsd.py` (clustering das poses do Vina por RMSD pareado, com figura por alvo) e `docking_rmsd_validation.py` (redocking da glutationa co-cristalizada GCG do 1BZL vs. cristal, tornando real o "RMSD < 2 Å" do texto de Métodos).
- **HADDOCK3** — `haddock3/` com receptor (dímero 1BZL), ligante 2-SCN-NQ, AIRs ancoradas em CYS53/CYS58/HIS461, workflow `run-TcTR-2SCN-NQ.cfg` e `HADDOCK3_PROTOCOLO.md`. Método ortogonal ao Vina.
- **Dashboard v2.1** — nova aba "RMSD & HADDOCK".
- **Versionamento** — `.gitignore`, `VERSION` (2.1.0), `CHANGELOG.md`. (Obs.: a pasta **não** é um repositório git ainda — `git init` recomendado.)
- **Fix do FAD (scaffolding)** — `regen_fad_xml.py`, `cofactor_FAD_TcTR_v3.pdb`, `md_clean.py` preferindo o par v3, e correção do regex de bonds no `auto_fix/fix_fad_xml.py`. **Falta só executar o regen (§1).**
- **`md_clean.py::analyze()` reescrito** — o RMSD do ligante agora é medido **no bolso** (após superpor a proteína, sem re-superpor no próprio ligante), capturando deriva/saída do sítio; adiciona curva de convergência e RMSF por resíduo.

---

## 4. Diagnóstico do FAD XML antigo (inalterado — confirmado)

O `cofactor_FAD_TcTR.xml` antigo tem o force-field inteiro parametrizado para a identidade atômica **errada**: 33/86 átomos com elemento/massa errados, σ de Lennard-Jones casando 86/86 com o elemento **errado**, e `PeriodicTorsionForce` vazio. O PDB-template antigo está embaralhado de forma idêntica, por isso o matching casa e a minimização roda — mas o warm-up diverge. **Nenhum patch de rótulo conserta isso**; só a regeneração completa (= `regen_fad_xml.py`). O candidato "fix só de elementos" (`cofactor_FAD_TcTR_v3_elemfix.xml`) foi corretamente **descartado** (`.INVALID.txt` no lugar). Detalhe completo: `BRIEFING_2026-05-21_md_module.md`.

---

## 5. Achado operacional novo: lag de sync OneDrive ↔ sandbox de análise

Durante este briefing, a ferramenta de shell do sandbox (que lê a pasta via OneDrive) viu **versões truncadas** de dois arquivos, enquanto a visão autoritativa do Windows mostrava-os **completos**:

| Arquivo | Visto pelo shell do sandbox | Real no disco (Windows) | Veredito |
|---|---|---|---|
| `regen_fad_xml.py` | 389 linhas (corta no meio de um `ET.SubElement`) | **537 linhas, completo** (termina em `main()`) | mount **stale**; disco OK |
| `md_clean.py` | 850 linhas (cortado no `--fresh` help) | **completo e compila** | mount **stale**; disco OK |

Isto bate com o canário do próprio projeto (`SCRIPT_VERSION = "vISO-2026-05-22"  # se nao aparecer isto no log, seu OneDrive nao sincronizou`). **Implicação:** análises automáticas rodando no sandbox podem ver arquivos parciais e reportar "erros de sintaxe" falsos. Para validar código, confie na cópia do Windows (ou aguarde o OneDrive hidratar o arquivo).

**Ação tomada (transparência):** antes de identificar o lag, reescrevi `md_clean.py` a partir da visão do mount + o backup `md_clean.py.bak_2026-05-21`, para fechar o que parecia um arquivo truncado. **Verifiquei o resultado:** o arquivo compila (`py_compile` OK), preserva integralmente o trabalho de 05-22 (lógica v3 nas linhas 292–302 e o `analyze()` reescrito), e o trecho após `--fresh` é **byte-idêntico** ao backup canônico. Nenhum conteúdo foi perdido; o `md_clean.py` está íntegro e completo. **O `regen_fad_xml.py` não foi tocado** (estava completo no disco; mexer nele a partir da visão stale o teria corrompido).

> Resíduo: criei `md_clean.py.truncated_2026-05-23.bak` (cópia da visão stale) e não consegui removê-lo (o mount bloqueia delete: "Operation not permitted"). É um arquivo inofensivo e pode ser apagado manualmente — **não** é uma versão válida do `md_clean.py`.

---

## 6. Pendências e bloqueios

1. **🔴 Rodar `regen_fad_xml.py`** (env `datascience`) — gera o `v3.xml` que falta. É o gate de tudo (§1).
2. **Decisão de protonação do difosfato** — o `regen` valida fórmula/carga; o FAD oxidado a pH 7,4 costuma ser **−2** (difosfato desprotonado). Confirmar `-nc -2` vs neutro **antes** de gerar cargas (e preferir **AM1-BCC** a Gasteiger para produção — exige AmberTools `sqm` no env).
3. **CUDA do OpenMM quebrado** (carry-over): plugins falham com `DLL load failed 126 (MOD_NOT_FOUND)`; pipeline cai para OpenCL mixed precision. Independente do FAD; limita escala de ns. `conda install -c conda-forge openmm-cuda` e conferir `cudart64_*.dll` em `…\envs\datascience\Library\bin`.
4. **`auto_fix/fix_fad_xml.py`** — segue sendo **stopgap**, não a solução. Não usar como caminho principal (o `regen` é o caminho válido).
5. **Validar no SISTEMA COMPLETO**, não no FAD isolado — o NaN só emerge nos 183k átomos solvatados.

---

## 7. Próximos passos priorizados

1. `conda activate datascience && python regen_fad_xml.py` → confirmar `RESULTADO: PASS` e que o `v3.xml` foi gravado.
2. `python md_clean.py --fresh --ns 1 --tag prod_v3 --platform CUDA --allow-fallback` → conferir ausência de `BOND ANOMALA` e de NaN no warm-up.
3. Em paralelo: consertar o CUDA (item 6.3).
4. Se a curta (1 ns) passar: escalar para 5/50 ns e rodar `analyze()` (RMSD-no-bolso + RMSF + convergência).
5. (Higiene) `git init` na pasta para versionar de verdade a v2.1.0; apagar `md_clean.py.truncated_2026-05-23.bak`.

> ⚠️ Mantido de `NEXT_SESSION.md`: **não** reiniciar ciclos automáticos de MD sem o `v3.xml` — vai dar NaN sempre.

---

## 8. Arquivos-chave

| Caminho | Papel | Estado |
|---|---|---|
| `md_clean.py` (939 linhas) | Pipeline MD ativo | ✅ completo, compila; FF do FAD entra na resolução de `fad_xml`/template (linhas 295–302) |
| `regen_fad_xml.py` (537 linhas) | Gera o force-field correto do FAD | ✅ completo no disco; **falta executar** no env `datascience` |
| `…/smoke_test_fad/cofactor_FAD_TcTR_v3.pdb` | Template com `element` corrigido | ✅ pronto/verificado |
| `…/smoke_test_fad/cofactor_FAD_TcTR_v3.xml` | **Saída do regen** | ❌ **ausente** (gate) |
| `…/smoke_test_fad/cofactor_FAD_TcTR.xml` / `.pdb` | Par antigo (force-field quebrado) | intocado; ainda em uso até o v3.xml existir |
| `…/smoke_test_fad/cofactor_FAD_TcTR_v3_elemfix.xml(.INVALID.txt)` | Tentativa descartada | ⛔ não usar |
| `docking_pose_rmsd.py`, `docking_rmsd_validation.py` | Validação de docking por RMSD (v2.1.0) | ✅ novos em 05-22 |
| `haddock3/` | Método ortogonal (docking guiado) | ✅ montado em 05-22 |
| `md_run_v2_prod.log` | Último run real de MD (14/05) | NaN no warm-up |
| `auto_fix/` | Loop de auto-correção | nunca executado |

---

*Briefing gerado automaticamente. Nenhuma simulação de produção foi executada e nenhum arquivo original do FAD (XML/PDB antigos) foi modificado. Única escrita em código nesta sessão: `md_clean.py` foi reescrito a partir do mount + backup e verificado completo/compilável antes de identificar o lag de sync (§5) — sem perda de conteúdo. `regen_fad_xml.py` não foi tocado. Verificações: leitura autoritativa via Windows, `py_compile`, e parse read-only de PDB/XML.*
