# Status Report — 2026-05-13 (sessão automatizada / hotfix FAD template)

**Tarefa:** retomar do limite de créditos da sessão anterior. **Resultado:**
identificado e corrigido o bug que travava `md_clean.py --fresh` durante
`addHydrogens` (crash `No template found for residue 970 (FAD)`).

---

## 1. Mudança arquitetural detectada (contexto)

Em 2026-05-12, o usuário aposentou toda a stack `md_module + cofactor_prep
+ smoke_test_fad_fix` para `obsoletos/` e migrou para `md_clean.py`
(antes 570 linhas) acionado pela rota `/run_md_clean` do `dashboard_v2.py`.

---

## 2. Crash observado às 19:09 (run do usuário, 2026-05-13)

```
2026-05-13 19:09:08  INFO  FAD #1 adicionado (protonado): 88 atomos
                            (cofactor_FAD_0_xtal_protonated_fixed.pdb)
2026-05-13 19:09:08  INFO  FAD #2 adicionado (protonado): 88 atomos
                            (cofactor_FAD_1_xtal_protonated_fixed.pdb)
2026-05-13 19:09:08  INFO  ligante adicionado: 36 atomos
2026-05-13 19:09:08  INFO  addHydrogens(pH=7.4)...
Traceback ...
  File "...openmm/app/modeller.py", line 1035, in addHydrogens
    system = forcefield.createSystem(newTopology, rigidWater=False,
                                      nonbondedMethod=CutoffNonPeriodic, ...)
ValueError: No template found for residue 970 (FAD). ...
```

`Modeller.addHydrogens` chama `ff.createSystem` internamente para validar
templates. Com o FAD já na topologia ANTES do `addHydrogens`, o matcher
do ForceField precisa casar o resíduo FAD contra o template do XML — e
falha porque os átomos não correspondem.

---

## 3. Causa raiz (diagnóstico programático)

Comparação dos átomos do `cofactor_FAD_TcTR.xml` (template SMIRNOFF) com
os PDBs disponíveis em `docking_project/figures/smoke_test_fad/`:

| Arquivo | Átomos | Nomes batem com XML? | Compatível? |
|---|---|---|---|
| `cofactor_FAD_TcTR.xml` (template) | 86 (53 heavy + 33 H) | — | — |
| `cofactor_FAD_TcTR.pdb` | 86 | **86/86** | **sim** (template ref) |
| `cofactor_FAD_0_xtal.pdb` | 53 (só heavy) | **53/53** | parcial (sem H) |
| `cofactor_FAD_1_xtal.pdb` | 53 (só heavy) | **53/53** | parcial (sem H) |
| `cofactor_FAD_0_xtal_protonated.pdb` | 88 | 53/54 (todos 35 H chamados `H`) | **não** (Hs duplicados) |
| `cofactor_FAD_0_xtal_protonated_fixed.pdb` | 88 | **12/88** (renomeados `P1/O1/C1`) | **não** (nomes errados) |
| `cofactor_FAD_1_xtal_protonated_fixed.pdb` | 88 | **12/88** | **não** |

O `md_clean.py --fresh` estava chamando `app.PDBFile` nos arquivos
`*_xtal_protonated_fixed.pdb` (88 átomos, nomes não-canônicos). Quando o
addHydrogens internamente roda createSystem, o matcher do FF não consegue
casar 88 átomos contra os 86 do template, e mesmo se contasse igual,
**74 dos 86 nomes do template não existem no PDB**.

---

## 4. Fix aplicado a `md_clean.py`

### 4.1 Helper `build_canonical_fad_topology(template, xtal, log_tag)`

Constrói uma topologia FAD com nomes canônicos (que casam com o XML) e
posições cristalográficas:

1. Lê `cofactor_FAD_TcTR.pdb` (86 átomos canônicos, IDÊNTICO ao XML).
2. Lê `cofactor_FAD_{i}_xtal.pdb` (53 átomos heavy, nomes canônicos).
3. Substitui as 53 posições heavy do template pelas do xtal (mesmo nome
   = mesmo átomo).
4. Aplica transformação rígida Kabsch (template_heavy → xtal_heavy) nos
   33 H do template, propagando-os para a localização correta do
   sítio ativo.
5. Devolve `(topology, positions)` prontos para `modeller.add(...)`.

Sanidade matemática validada em pure-numpy: para transformação identidade,
RMSD = 0 e erro de H = 3.55e-15 nm. Para transformação random,
RMSD heavy < 0.001 Å (Kabsch é exato para rígido).

### 4.2 Reorder em `build_system_fresh`

Ordem nova (fix 2026-05-13):

```
1) PDBFixer no receptor (remove HOH/GCG/FAD)
2) Modeller (proteína apenas)
3) modeller.addHydrogens(ff, pH=7.4)        ← ANTES de FAD/LIG
4) Para cada cofactor_FAD_{i}_xtal.pdb (53 atomos):
       build_canonical_fad_topology(template, xtal)
       modeller.add(top, pos)                ← 86 atomos canonicos
5) Adicionar ligante (lig_TcTR.pdb)
6) addSolvent (TIP3P, 1 nm, 0.15 M NaCl)
7) createSystem (PME)
```

A ordem é crítica: addHydrogens chama internamente `createSystem`, que
validaria TODOS os resíduos contra os templates do FF; FAD/LIG têm
templates rígidos com átomos canônicos que precisam casar exatamente.
Mantendo FAD/LIG fora da topologia até depois do addHydrogens, o
matcher só precisa lidar com aminoácidos padrão (que o
`amber14-all.xml` resolve naturalmente).

### 4.3 Wiring no `main()` — defaults atualizados

`main()` agora aponta os defaults para os arquivos canônicos:

```python
fad1 = args.fad1 or (sub / "cofactor_FAD_0_xtal.pdb")   # 53 heavy canonicos
fad2 = args.fad2 or (sub / "cofactor_FAD_1_xtal.pdb")   # 53 heavy canonicos
```

(antes apontavam para `cofactor_FAD_*_xtal_protonated_fixed.pdb`).

Auto-detecção dentro do loop: se o PDB tem 53 átomos, chama o helper;
se já tem 86 canônicos, usa direto. Os arquivos `_protonated_fixed.pdb`
(88 átomos) ainda podem ser apontados explicitamente via `--fad1`/`--fad2`
mas vão cair no fallback `app.PDBFile(...)` direto — que continuará
falhando se eles forem usados.

### 4.4 Import robusto de `Vec3`

```python
try:
    from openmm import Vec3        # 8.x
except ImportError:
    from openmm.vec3 import Vec3   # 7.x fallback
```

---

## 5. Estado do arquivo após o fix

| Métrica | Antes (03:25) | Após fix |
|---|---|---|
| Linhas | 570 (truncado em 559 no disco) | 709 |
| Funções top-level | 7 (sem helper) | 8 (com `build_canonical_fad_topology`) |
| AST | OK | **OK** |
| `__main__` guard | OK | OK |
| Helpers Kabsch | — | Validados em pure-numpy |

Verificação programática:
```
Functions: ['build_canonical_fad_topology', 'build_system_fast',
            'build_system_fresh', 'pick_platform', 'add_position_restraints',
            'diagnose_system', 'run_md', 'analyze', 'main']
Has __main__ guard:        True
Vec3 try/except fallback:  True
build_canonical_fad call:  True
addHydrogens([protein]):   True
main() defined:            True
Ends with main():          True
```

---

## 6. Plano de validação — comandos para o usuário

PowerShell limpo, sem Python rodando do run de 19:09:

```powershell
Get-Process python -ErrorAction SilentlyContinue | Stop-Process -Force
Get-ChildItem -Path 'C:\Users\luizd\OneDrive\Data_Science\Projetos\docking' `
  -Recurse -Include 'md_clean*.pyc','__pycache__' -ErrorAction SilentlyContinue |
  Remove-Item -Force -Recurse

conda activate datascience
cd C:\Users\luizd\OneDrive\Data_Science\Projetos\docking

# 1) Sanidade AST
python -c "import ast; ast.parse(open('md_clean.py',encoding='utf-8').read()); print('md_clean OK')"

# 2) Smoke test rapido (1 ns para validar pipeline ate a producao)
python md_clean.py --fresh --ns 1 --tag dryrun_fix --platform CUDA
```

**Logs esperados se o fix funcionar:**

```
[FRESH] receptor=1BZL_clean.pdb fad1=cofactor_FAD_0_xtal.pdb
        fad2=cofactor_FAD_1_xtal.pdb lig=lig_TcTR.pdb
...
addHydrogens(pH=7.4) [protein-only]...
proteina protonada: ~7600-7800 atomos (was 7459 heavy)
FAD #1: xtal heavy-only (53 atomos); canonicalizando via Kabsch...
[canonical] FAD#1: Kabsch RMSD heavy = 0.000 A (n_heavy=53, n_H_propagated=33)
FAD #1 adicionado (canonico): 86 atomos via template cofactor_FAD_TcTR.pdb
FAD #2: xtal heavy-only (53 atomos); canonicalizando via Kabsch...
[canonical] FAD#2: Kabsch RMSD heavy = 0.000 A (n_heavy=53, n_H_propagated=33)
FAD #2 adicionado (canonico): 86 atomos via template cofactor_FAD_TcTR.pdb
ligante adicionado: 36 atomos
addSolvent (TIP3P, padding 1.0 nm, 0.15 M NaCl)...
sistema final: ~177 000 atomos
cofator FAD: 2  |  ligante LIG: 1  |  HOH: ~161 000  |  Na/Cl: ~300
Plataforma SELECIONADA: CUDA (props={'Precision': 'single', ...})
...
[OpenMM] Producao OK em X min
RMSD/RMSF salvo em md_dryrun_fix_rmsd_rmsf.png
DONE.
```

---

## 7. Cenários de falha residual a observar

| Sintoma | Causa possível | Mitigação |
|---|---|---|
| `No template found for residue (FAD)` ainda aparece | `cofactor_FAD_TcTR.pdb` ausente ou tem != 86 atomos | Conferir arquivo; se necessário re-gerar com PDBFixer + addH no FAD canonico do RCSB |
| `Particle coordinate is NaN` no warm-up | Hs propagados via Kabsch caem em colisão com proteína | Aumentar minimizacao (default já é `LocalEnergyMinimizer` em `run_md`); se persistir, reintroduzir warm-up dt=0.5fs |
| FAD na topologia mas createSystem falha | XML usa types que o FF não conhece | Validar XML: `python -c "from openmm.app import ForceField; ForceField('cofactor_FAD_TcTR.xml')"` |
| `--fast` ainda usa complex sem FAD | `complex_TcTR_initial.pdb` é o arquivo bugado de 12/05 | Já mitigado: dashboard default é `--fresh`. Para `--fast`, gerar novo complex_initial pelo próprio --fresh primeiro e copiar |

---

## 8. Arquivos modificados nesta sessão

- [md_clean.py](computer://C:\Users\luizd\OneDrive\Data_Science\Projetos\docking\md_clean.py) — adicionado `build_canonical_fad_topology` (Kabsch align template→xtal); reordenado `build_system_fresh` para addHydrogens(protein) ANTES de adicionar FAD/LIG; main() agora usa `cofactor_FAD_*_xtal.pdb` por default. Vec3 import com fallback 7.x/8.x.
- [STATUS_REPORT_2026-05-13.md](computer://C:\Users\luizd\OneDrive\Data_Science\Projetos\docking\STATUS_REPORT_2026-05-13.md) — este relatório (atualizado).

---

## 9. Resumo executivo

O crash de 19:09 (`No template found for residue 970 (FAD)`) tinha duas
causas combinadas: (a) `Modeller.addHydrogens` chamava internamente
`createSystem` na topologia inteira, exigindo template match para o FAD
**antes** dos Hs estarem montados; (b) os arquivos
`cofactor_FAD_*_xtal_protonated_fixed.pdb` haviam sido renomeados
(P1/O1/C1...) de forma incompatível com o template
`cofactor_FAD_TcTR.xml` (PA/O1A/PA...).

**Fix:** adicionado helper `build_canonical_fad_topology(template, xtal)`
que faz Kabsch align entre o template canônico (86 átomos) e o xtal
heavy (53 átomos), produzindo um FAD com **nomes canônicos + posições
cristalográficas**. Reordenado `build_system_fresh` para rodar
`addHydrogens` somente na proteína, ANTES de adicionar FAD/LIG via
o novo helper.

**AST OK, 709 linhas, 9 funções top-level, lógica Kabsch validada em
pure-numpy.** Próximo passo do usuário: rodar
`python md_clean.py --fresh --ns 1 --tag dryrun_fix --platform CUDA`
(~3 min de produção em RTX 3060) para confirmar a MD completa
end-to-end com os dois FADs corretamente canonicalizados.
