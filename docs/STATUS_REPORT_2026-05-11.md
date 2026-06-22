# Status Report — 2026-05-11 (sessão automatizada)

**Tarefa:** retomar do limite de créditos da sessão anterior — foco em
`md_module.py` e na cadeia de protonação de cofatores (`cofactor_prep.py`).

---

## 1. Estado encontrado no disco

| Arquivo | Tamanho | Última edição | Estado |
|---|---|---|---|
| `md_module.py` | 170 806 B | 2026-05-11 22:45 | OK (AST válido), fingerprint `2026-05-11.v4.9.pdbfixer-cofactor` |
| `cofactor_prep.py` | 26 112 B | 2026-05-11 22:45 | **TRUNCADO** — terminava em `de` no meio de uma definição |
| `smoke_test_fad_fix.py` | 10 958 B | 2026-05-11 22:45 | **TRUNCADO** — última linha `sys.exit` (sem parênteses) |
| `test_cofactor_fad.py` | 4 218 B | 2026-05-11 15:38 | OK |
| `md_module.py.broken_2026-05-07_backup` | 141 412 B | — | backup mantido |

A inspeção confirmou que a sessão anterior **parou no meio da escrita** de
`cofactor_prep.py` (parou ao escrever o wrapper retro-compat
`prepare_cofactor_openff_mol → _LegacyMolWrapper._Conf`) e provavelmente
também no fechamento do `smoke_test_fad_fix.py`.

O `cofactor_prep.py` por acaso continuava sendo importável (o `de` solto
era interpretado como uma referência a identificador num corpo de método),
mas qualquer chamada a `prepare_cofactor_openff_mol(...)` resultaria em
`NameError: name 'de' is not defined`. O `md_module.py` v4.9 não chama
o wrapper legado — só a API nova `prepare_cofactor_topology_positions`.
Portanto **a MD funcionava se nenhum código antigo tocasse o wrapper**,
mas o módulo estava semanticamente quebrado.

---

## 2. Correções aplicadas nesta sessão

### Fix #1 — fechar truncamento do `cofactor_prep.py`

Substituí o `_LegacyMolWrapper` parcial por um **wrapper minimalista
explicitamente DEPRECATED**: ele expõe `topology`, `positions`, `n_atoms`,
`n_bonds`, `n_conformers` e um método `to_topology_positions()` que retorna
a tupla nova. Emite `DeprecationWarning` no uso. Removida toda a tentativa
incompleta de mimar a API do `openff.Molecule` (`.atoms[i].element`,
`.conformers[0].m_as("nanometer")`), que estava cortada no meio.

Adicionei também ao final do módulo uma função `_module_smoke_check()`
chamada automaticamente no `import` — verifica que os símbolos críticos
(`prepare_cofactor_topology_positions`, `_openff_mol_to_openmm`,
`__version__`) existem. Isso protege contra cache `.pyc` stale de uma
versão truncada anterior.

### Fix #2 — fechar truncamento do `smoke_test_fad_fix.py`

Última linha era `sys.exit` (referência à função, não chamada). Trocado
por `sys.exit(9)` na branch `[FALHOU]` + impressão do `data["error"]` se
presente, para diagnóstico imediato em caso de falha.

### Validação AST

```
md_module.py            → AST OK (3 402 linhas)
cofactor_prep.py        → AST OK (581 linhas, depois do fix)
smoke_test_fad_fix.py   → AST OK (261 linhas, depois do fix)
test_cofactor_fad.py    → AST OK (136 linhas)
```

---

## 3. Arquitetura atual (v4.9 — `pdbfixer-cofactor`)

A pipeline de protonação de cofator é hoje um **cascade em camadas** dentro
de `cofactor_prep.prepare_cofactor_topology_positions()`:

1. **PDBFixer** (`addMissingHydrogens(pH=7.0)`) — primeira tentativa. Para o
   FAD produz 86 átomos (53 heavy + 33 H), exato match com o template
   SMIRNOFF do `md_module`. É a rota canônica para Windows sem AmberTools.
2. **SMILES + Hungarian + Kabsch** — para cofatores em que PDBFixer não tem
   template interno. Faz match heavy↔heavy por elemento via
   `scipy.optimize.linear_sum_assignment`, alinha rigidamente, e propaga a
   transformação para todos os átomos (heavy + H) do template embed-3D do
   RDKit. RMSD esperado < 2.0 Å heavy-heavy.
3. **OpenBabel pybel → PDB → openmm.app.PDBFile** — terceira tentativa, com
   `correctForPH=False` (evita os 2 H extras nos fosfatos do FAD que o
   OpenBabel adicionaria com `correctForPH=True`).
4. **openff.toolkit.Topology.from_pdb** — só se o PDB de entrada já tiver Hs.
5. **RDKit AssignBondOrdersFromTemplate + AddHs** — fallback final em 4
   variantes (CONECT vs proximity, sanitize on/off).

Se todas falharem, é gravado `cofactor_<NAME>_prep_errors.log` no mesmo
diretório do PDB de entrada com o motivo de cada estratégia.

O `md_module.py` v4.9 consome o resultado em `run_openmm_real`
(linhas ≈ 1730–1806), executando o ciclo:

```text
1) Salvar cofatores cristalográficos + remover do modeller
2) modeller.addHydrogens(ff, pH=7.4)              ← agora só na proteína
3) Para cada cofator salvo:
     cofactor_prep.prepare_cofactor_topology_positions(...)
     modeller.add(top, pos)                       ← re-injeção com Hs
4) Adicionar ligante
5) Solvatar (TIP3P, 8 Å) e simular
```

---

## 4. Riscos/regressões a observar na próxima execução

- **Equilibração simplificada.** O v4.9 voltou a ter só **2 fases**
  (50 ps @ 1 fs + 150 ps @ 4 fs). O `warm-up dt=0.5 fs @ 100 K + T-ramp`
  da v3.2 e a `fase 1.5 dt=2 fs` da v2.2 foram **removidos**. Isso provavelmente
  é intencional (PDBFixer dá Hs em posição canônica, com energia muito mais
  baixa que o ICP/Kabsch da v3.x). Se reaparecer `Particle coordinate is NaN`
  na fase 1, reintroduzir a fase 0 conforme `BRIEFING_2026-05-08_md_module_v3.2.md`.
- **Rota E ainda PREVIEW.** Para publicação, instalar `ambertools` para
  AM1-BCC + GAFF2 (Rotas A/B/C ativas).
- **OneDrive sync lag.** Os arquivos editados nesta sessão vão para
  `C:\Users\luizd\OneDrive\Data_Science\Projetos\docking\`. O cliente OneDrive
  do Windows reescreve com latência de segundos. Antes de rodar o smoke test,
  conferir com `Get-Item md_module.py | Select Length, LastWriteTime`.

---

## 5. Plano de validação — comandos prontos

Execute em um PowerShell limpo, com **todos os processos Python fechados**
(dashboard, jupyter):

```powershell
Get-Process python -ErrorAction SilentlyContinue | Stop-Process -Force
Get-ChildItem -Path C:\Users\luizd -Recurse -Include 'md_module*.pyc','cofactor_prep*.pyc' -ErrorAction SilentlyContinue | Remove-Item -Force

conda deactivate
conda activate datascience

cd C:\Users\luizd\OneDrive\Data_Science\Projetos\docking

# Sanity: arquivos parseáveis
python -c "import ast; ast.parse(open('md_module.py',encoding='utf-8').read()); print('md_module OK')"
python -c "import ast; ast.parse(open('cofactor_prep.py',encoding='utf-8').read()); print('cofactor_prep OK')"
python -c "import ast; ast.parse(open('smoke_test_fad_fix.py',encoding='utf-8').read()); print('smoke OK')"

# Import isolado do cofactor_prep (dispara _module_smoke_check)
python -c "import cofactor_prep; print('cofactor_prep v', cofactor_prep.__version__, '- OK')"

# Smoke test completo (5 ns TcTR, ~3-5 min na RTX 3060)
python smoke_test_fad_fix.py
```

**Critérios de sucesso (saída esperada do smoke test):**

```
[OK] OpenMM 8.2 ...
[OK] CUDA operacional (NVIDIA GeForce RTX 3060)
[AST OK] md_module.py
[AST OK] docking_v2.py
[OK] result_TcTR.pdbqt
md_module versao   = 2026-05-11.v4.9.pdbfixer-cofactor
[OK] Versao correta importada
...
[cofactor_prep] FAD: PDBFixer (pH 7) OK (86 atoms = 53 heavy + 33 H, 91 bonds)
[Cofator->H] FAD #1: 53 heavy + 33 H via cofactor_prep v2.0
[Cofator->H] FAD #2: 53 heavy + 33 H via cofactor_prep v2.0
[OK] 2 cofator(es) re-adicionado(s) ao modeller via cofactor_prep
[OK] Ligand added: residuo={'LIG'} | 36 atomos
[OpenMM] Energia pos-minimizacao: -3.0e6 kJ/mol
[OpenMM] NVT equilibration (fase 1: 50 ps, dt=1 fs) ...
[OpenMM] NVT equilibration (fase 2: 150 ps, dt=4 fs) ...
... (production run)
[PASSOU] FAD foi parametrizado e a MD real completou.
```

**Cenários de falha residual e mitigação:**

| Sintoma | Causa provável | Mitigação |
|---|---|---|
| `cofactor_prep incompleto` no import | Cache `.pyc` da versão truncada | `Remove-Item __pycache__\cofactor_prep*.pyc` |
| `[cofactor_prep] FAD: ATOM COUNT MISMATCH` (Estratégia 0) | OpenBabel super-protonou | Estratégias 1/2 assumem; sem ação |
| `Particle coordinate is NaN` na fase 1 NVT | Cofator com Hs ainda fora do bolso | Reintroduzir warm-up dt=0.5fs (v3.2) — código já versionado |
| `No template found for residue FAD` voltou | XML do FAD não foi registrado no ForceField | Verificar `_cofactors_registered`; whitelist em `run_openmm_real` |
| `Versao do md_module em memoria != EXPECTED` | Outro processo Python já tem v1 carregado | Matar python.exe, limpar `__pycache__\` |

---

## 6. Arquivos modificados nesta sessão

- [cofactor_prep.py](computer://C:\Users\luizd\OneDrive\Data_Science\Projetos\docking\cofactor_prep.py) — fechado truncamento + smoke-check no import
- [smoke_test_fad_fix.py](computer://C:\Users\luizd\OneDrive\Data_Science\Projetos\docking\smoke_test_fad_fix.py) — `sys.exit(9)` na branch [FALHOU]
- [STATUS_REPORT_2026-05-11.md](computer://C:\Users\luizd\OneDrive\Data_Science\Projetos\docking\STATUS_REPORT_2026-05-11.md) — este relatório

---

## 7. Resumo executivo

A sessão anterior terminou no meio da escrita de `cofactor_prep.py` (último
arquivo a ser modificado pelo limite de créditos). O `md_module.py` v4.9
estava completo e parseável, mas a chamada de `prepare_cofactor_openff_mol`
(API legada) explodiria com `NameError`. Os dois truncamentos foram fechados,
o smoke test agora retorna código de saída corretamente, e o módulo
auto-verifica integridade no import. **Próximo passo do usuário:** rodar
`python smoke_test_fad_fix.py` no env `datascience` para validar a MD real
do TcTR de ponta a ponta.
