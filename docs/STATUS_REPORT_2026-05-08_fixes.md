# Status Report — 2026-05-08 (segunda execução)
**Foco:** análise do output da MD `2026-05-08_11-25-29` + fixes aplicados

---

## 1. Análise do output que você forneceu

### Pontos POSITIVOS (do briefing matinal)

| Item | Estado |
|---|---|
| Pose docked TcTR | ✓ Encontrada (`result_TcTR.pdbqt` em `docking_project/docking/`) |
| Gate Docking→MD | ✓ Funcionou (`[OK] Todas as poses docked encontradas`) |
| OpenMM 8.2 + CUDA RTX 3060 | ✓ Detectado |
| Detecção do FAD | ✓ `[Cofator] Detectados: FAD` |
| Parametrização SMIRNOFF do FAD | ✓ `OK: 86 átomos, 91 bonds, 157 angles` |
| XML do FAD carregado no ForceField | ✓ `[OK] Cofator FAD parametrizado e registrado (86 átomos)` |
| Conversão do ligante PDBQT→SDF | ✓ `36 atoms` |
| Whitelist com FAD mantido | ✓ `Cofatores parametrizados mantidos: FAD` |

### Bugs/avisos CONHECIDOS (não bloqueadores)

1. `AmberTools/OpenEye não encontrados` — usa Gasteiger (RDKit) em vez de AM1-BCC. Aceitável para preview; para publicação, instalar `ambertools`.
2. `openff.interchange ausente` → pulou Rota 0.
3. `antechamber ausente no PATH` → pulou Rotas A/B/B2/C.
4. `Found multiple definitions for atom type: LGT0` na Rota D — caiu na Rota E (RDKit/MMFF94 PREVIEW). **Para publicação**, instalar `openff-toolkit`:
   ```powershell
   conda install -c conda-forge openff-toolkit
   ```

### BUG CRÍTICO identificado (era o bloqueador)

```
[ERRO] OpenMM simulation failed: addHydrogens falhou:
No template found for residue 970 (FAD).
The set of atoms is similar to FAD, but it is missing 33 hydrogen atoms.
```

**Causa-raiz diagnosticada:**

- O FAD no PDB cristalográfico (`1BZL.pdb`) tem só **53 heavy atoms** (sem hidrogênios — padrão de cristalografia de proteínas).
- O template SMIRNOFF gerado pelo `_extract_and_parametrize_cofactors` cobre **86 átomos** (53 heavy + 33 H, vindos do `RDKit.AddHs + EmbedMolecule`).
- Quando `Modeller.addHydrogens()` busca o template do FAD para inserir os Hs faltantes no PDB, o template registrado **não tem `<Hydrogens>` definitions** (seção do XML que diz "para cada heavy atom X, adicionar Hs ligados Y1, Y2"). Sem isso, OpenMM detecta que faltam Hs mas não sabe como adicioná-los.
- Resultado: `addHydrogens` aborta, MD cai em `simulate_md_data` (sintético), o gráfico final mostra `[sintetico]`.

---

## 2. Fixes aplicados nesta sessão

### Fix #1 — `md_module.py`: estratégia "remove → addHydrogens → re-add com Hs"

**Localização:** `run_openmm_real`, em torno do passo 6 (linhas ~1482-1583).

**Lógica nova:**

1. **6a (NOVO)** — Antes de `addHydrogens`, salvar topologia+posições dos cofatores cristalográficos do modeller, depois removê-los temporariamente via `modeller.delete()`.
2. **6 (já existia)** — `modeller.addHydrogens(ff, pH=7.4)` agora roda só na proteína padrão AMBER14 (que sempre funciona).
3. **6a-bis (NOVO)** — Para cada cofator salvo, ler o PDB protonado já gerado pelo `_extract_and_parametrize_cofactors` (`cofactor_<NAME>_<TID>.pdb`, com 86 átomos), substituir as coords dos heavy atoms pelas cristalográficas originais do receptor (preservando geometria do sítio ativo), e re-adicionar ao modeller via `modeller.add(topology, positions)`.

**Por que funciona:** o `addHydrogens` nunca vê o cofator, então o erro "No template found" não acontece. Os Hs do cofator vêm do RDKit (`AddHs + EmbedMolecule + MMFFOptimizeMolecule`), com geometria 3D razoável; a parte heavy mantém as coords cristalográficas, preservando interações com sítio ativo.

**Output esperado nos logs após o fix:**
```
[Cofator->H] Salvando FAD (53 heavy atoms)
[OpenMM] Adding hydrogens (pH 7.4) ...
[OpenMM] Hydrogens added: ... (total H: ...)
[Cofator->H] FAD: 53 heavy alinhados ao receptor, 33 Hs do RDKit
[OK] Cofatores re-adicionados ao modeller (com Hs preservando geometria cristalografica)
```

### Fix #2 — `docking_v2.py`: hardening preventivo (rota RDKit fallback)

**Localização:** `_convert_to_pdbqt`, rota B.

**Mudanças:**

1. Trocado `pdb_lines[i].rstrip("\n")` por `rstrip("\r\n")` na L582-584.
   - Razão: arquivos PDB com `\r\n` (Windows) deixavam `\r` residual que invadia o slice `[:66]`, deslocando o campo charge.
2. Adicionado **self-test** ao final da função: lê a primeira linha ATOM/HETATM gerada e valida `float(line[70:76])`. Se falhar → retorna falha imediata em vez de deixar Vina rejeitar 9× depois.

**Status:** ambos no arquivo `docking_v2.py` atual. Bug do briefing matinal não deve mais voltar.

### Fix #3 — Smoke test criado

**Arquivo:** `smoke_test_fad_fix.py`

Roda MD curta (5 ns) para TcTR e verifica se:
- AST de `md_module.py` e `docking_v2.py` estão OK
- `result_TcTR.pdbqt` existe
- MD chega em `real_sim=True` (não cai em fallback)
- FAD aparece em `cofactors_registered`

Uso:
```powershell
conda activate datascience
python smoke_test_fad_fix.py
```

---

## 3. Plano de validação

### Passo 1 — Limpar __pycache__ (forçar reload do módulo)

```powershell
cd C:\Users\luizd\OneDrive\Data_Science\Projetos\docking
Remove-Item __pycache__\md_module.cpython-*.pyc -ErrorAction SilentlyContinue
Remove-Item __pycache__\docking_v2.cpython-*.pyc -ErrorAction SilentlyContinue
```

### Passo 2 — Validar AST dos arquivos modificados

```powershell
python -c "import ast; ast.parse(open('md_module.py', encoding='utf-8').read()); print('md_module OK')"
python -c "import ast; ast.parse(open('docking_v2.py', encoding='utf-8').read()); print('docking_v2 OK')"
```

### Passo 3 — Smoke test (5 ns TcTR)

```powershell
python smoke_test_fad_fix.py
```

Esperado: `[PASSOU] FAD foi parametrizado e a MD real completou.`

### Passo 4 — Se o smoke test passar, rodar o restante via dashboard

```powershell
python dashboard_v2.py
# Abrir http://localhost:5001 → aba MD → modo real → 50 ns multi-alvo
```

---

## 4. Saídas/arquivos do projeto

### Arquivos modificados nesta sessão

- [md_module.py](computer://C:\Users\luizd\OneDrive\Data_Science\Projetos\docking\md_module.py) — fix FAD: passos 6a (salvar+remover) e 6a-bis (re-add com Hs)
- [docking_v2.py](computer://C:\Users\luizd\OneDrive\Data_Science\Projetos\docking\docking_v2.py) — rstrip `\r\n` + self-test PDBQT

### Arquivos novos

- [smoke_test_fad_fix.py](computer://C:\Users\luizd\OneDrive\Data_Science\Projetos\docking\smoke_test_fad_fix.py) — validação rápida do fix
- [STATUS_REPORT_2026-05-08_briefing.md](computer://C:\Users\luizd\OneDrive\Data_Science\Projetos\docking\STATUS_REPORT_2026-05-08_briefing.md) — briefing matinal
- [STATUS_REPORT_2026-05-08_fixes.md](computer://C:\Users\luizd\OneDrive\Data_Science\Projetos\docking\STATUS_REPORT_2026-05-08_fixes.md) — este relatório

---

## 5. Riscos remanescentes

1. **Rota E é PREVIEW** — para publicação, instalar `openff-toolkit` (Rota D ou 0). Linha de comando do output sugere: `conda install -c conda-forge openff-toolkit`.
2. **Hs do FAD pelo RDKit** — geometria razoável mas não cristalográfica. Para precisão máxima de FAD em sítio ativo, usar `pdb4amber` ou pré-protonar com `propka` antes de carregar.
3. **Outros cofatores não testados** — fix funciona genericamente para `_cofactors_registered` (FAD, FMN, NADP+, NADPH, NAD+, HEM). Para LmDHODH (FMN+FAD) e hDHODH (FMN), testar separadamente.
4. **Cache `.pyc`** — sempre limpar antes de re-rodar para garantir que a versão nova do `md_module.py` é carregada.

---

## 6.1 Iteração v2 — após primeiro smoke test

### O que aconteceu

O smoke test rodou no env correto (`datascience` com OpenMM 8.2 + CUDA RTX 3060 + RDKit 2025.09.5 + openff-toolkit 0.18.0). **Avancou muito mais longe que antes:**

```
[Cofator->H] Salvando FAD (53 heavy atoms)        ← FAD chain A
[Cofator->H] Salvando FAD (53 heavy atoms)        ← FAD chain B (1BZL é dimero)
[OpenMM] Adding hydrogens (pH 7.4) ...
[OpenMM] Hydrogens added: 7,483                   ← addHydrogens NÃO falhou mais!
[Cofator->H] FAD: 53 heavy alinhados ao receptor, 33 Hs do RDKit  ← 2x
[OK] Cofatores re-adicionados ao modeller
[OK] Ligand added: residuo={'LIG'} | 36 atomos
[OpenMM] System: 166,268 atoms (protein+ligand+solvent)
[OpenMM] Energia pos-minimizacao: -3139431 kJ/mol  ← OK
[OpenMM] NVT equilibration (fase 1: 50 ps, dt=1 fs) ...  ← passou
[OpenMM] NVT equilibration (fase 2: 150 ps, dt=4 fs) ...
[ERRO] OpenMM simulation failed: Particle coordinate is NaN.
```

**Bug novo identificado:** explosão de NaN na transição 1 fs → 4 fs da NVT.

### Causa-raiz suspeita

O fix v1 substituiu **as coords heavy** do PDB protonado pelas cristalográficas, mas **deixou os Hs nas coords originais do RDKit Embed3D**. Resultado: Hs flutuando dezenas de Å fora dos heavy atoms aos quais estão covalentemente ligados → forças astronômicas → integrator HMR a 4 fs não suporta → NaN.

### Fixes v2 aplicados nesta sessão

#### Fix v2.1 — Kabsch alignment (geometria H-heavy preservada)

Em vez de substituir só os heavy, agora calculamos a **rotação+translação ótima** (algoritmo Kabsch via SVD) que mapeia os heavy atoms do PDB protonado para as coords cristalográficas, e aplicamos essa transformação em **todos os átomos** (heavy + Hs). Resultado: Hs ficam a ~1 Å dos heavy aos quais estão ligados (geometria física), preservando a topologia molecular intacta.

**Output esperado nos logs (com o fix):**
```
[Cofator->H] FAD #1: 53 heavy alinhados (RMSD=0.012 A), 33 Hs em geometria local
[Cofator->H] FAD #2: 53 heavy alinhados (RMSD=0.018 A), 33 Hs em geometria local
[OK] 2 cofator(es) re-adicionado(s) ao modeller (Kabsch align — geometria H-heavy preservada)
```

RMSD < 0.02 Å significa o alinhamento heavy-heavy é praticamente perfeito (limite de precisão das coords PDB), e os Hs vêm "carregados" nessa transformação rígida.

#### Fix v2.2 — Rampa de timestep 1 → 2 → 4 fs

Inserida fase intermediária na NVT:
- **Fase 1:** 50 ps @ 1 fs (resolve instabilidades iniciais)
- **Fase 1.5 (NOVA):** 25 ps @ 2 fs (transição suave)
- **Fase 2:** 125 ps @ 4 fs (HMR full speed)

Total ainda = 200 ps de equilibração.

#### Fix v2.3 — NaN check sub-step

Cada fase NVT é dividida em sub-steps com checagem de energia (`isnan`/`isinf`) entre eles. Se explodir, o erro aponta exatamente qual sub-fase falhou (10 ps de fase 1, 25 ps de fase 1.5, ou 25 ps de fase 2). Isso permite diagnóstico cirúrgico em iterações futuras.

**Output esperado (caso de sucesso):**
```
[OpenMM] NVT equilibration (fase 1: 50 ps, dt=1 fs) ...
    [OK] dt=1fs estavel (-3.0e6 kJ/mol apos 10 ps)
[OpenMM] NVT equilibration (fase 1.5: 25 ps, dt=2 fs) ...
[OpenMM] NVT equilibration (fase 2: 125 ps, dt=4 fs) ...
```

**Output esperado (caso de falha cirúrgica):**
```
RuntimeError: Energia NaN/inf apos NVT fase 2 sub-step 3/5 (25ps cada @ 4fs)
```

### Como testar v2

```powershell
conda activate datascience
cd C:\Users\luizd\OneDrive\Data_Science\Projetos\docking
python smoke_test_fad_fix.py
```

Resultado esperado: `[PASSOU] FAD foi parametrizado e a MD real completou.` em ~5 min.

Se ainda explodir, o log agora dirá exatamente qual sub-fase, o que permite:
- Se for em fase 1 sub-step 1 (10 ps @ 1 fs) → problema de Hs ainda mal posicionados
- Se for em fase 1.5 (2 fs) → energia ainda alta, precisa minimização extra
- Se for em fase 2 (4 fs) → HMR muito agressivo, considerar dt=3 fs final

---

## 6. Resumo executivo

- **Bug bloqueador identificado:** `Modeller.addHydrogens()` falhando com FAD porque o template SMIRNOFF não tem `<Hydrogens>` definitions.
- **Fix aplicado:** estratégia "remove cofatores → addHydrogens só na proteína → re-adiciona cofatores com Hs do RDKit + coords cristalográficas".
- **Hardening preventivo:** `docking_v2.py` agora resiste a CRLF residual e trava cedo se o PDBQT sair mal-formado.
- **Smoke test pronto:** `python smoke_test_fad_fix.py` valida tudo em ~5 min.
- **Próximo passo:** rodar smoke test → se passar, escalar para 50-100 ns multi-alvo.
