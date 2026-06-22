# NEXT_SESSION — Docking (MD de TcTR + 2-SCN-NQ)

**Atualizado:** 2026-05-14
**Última sessão:** Tentativa de produção MD; root cause do bug isolada

---

## Status atual

### ✅ Pipeline FRESH end-to-end funcionando (até minimização)

- Sistema reconstrói: PDBFixer → addHydrogens(protein) → 2× FAD canônicos via Kabsch → LIG → addSolvent(TIP3P + 0.15 M NaCl) = **183.179 átomos**
- Inventário: TcTR dímero + 2 FAD + 1 LIG (2-SCN-NQ) + 55.906 HOH + 311 Na/Cl
- Minimização converge: PotentialEnergy = -3.00 × 10⁶ kJ/mol
- PDB minimizado salvo: `docking_project/figures/md_clean/minimized_dryrun_v2_prod.pdb` (15 MB)

### ❌ Produção MD não inicia — NaN no warm-up

Forensic check pós-minimização aponta:
```
[forensic] BOND ANOMALA 2.06 A: N3A/FAD1 - C4A/FAD1
[forensic] forca max em ('FAD', 'LIG'): 7.83e+02 kJ/mol/nm
[forensic] resumo: 2 bonds anomalas, 0 clashes <1.0 A
[forensic] >>> geometria ruim em FAD/LIG detectada (provavel causa do NaN no warm-up)
```

Após dt=0.25 fs @ 5 K, NaN imediato.

### 🎯 Root cause isolada — bug no `cofactor_FAD_TcTR.xml`

Não é bug do `md_clean.py` v2. O XML do FAD em `docking_project/figures/smoke_test_fad/cofactor_FAD_TcTR.xml` tem **mapeamento atom-element-tipo scrambled** para 9 dos 10 átomos do anel adenina:

| Atom name | Tipo atribuído | Element no Type | Correto |
|---|---|---|---|
| N9A | FAD12 | C | N |
| C8A | FAD13 | C | C ✓ |
| N7A | FAD14 | O | N |
| C5A | FAD15 | C | C ✓ |
| C6A | FAD16 | O | C |
| N6A | FAD17 | C | N |
| N1A | FAD18 | O | N |
| C2A | FAD19 | C | C ✓ |
| N3A | FAD20 | O | N |
| C4A | FAD21 | P | C |

E pior: `HarmonicBondForce` **não define bond** entre FAD20 (N3A) e FAD21 (C4A). Por isso o anel adenina não fecha — a "bond" N3A-C4A só é mantida pela 1-3 14 do dihedral e por forças repulsivas, daí os 2.06 Å.

### ⚠️ Bonus: CUDA plugin OpenMM falhou

Plugin CUDA do OpenMM falhou ao carregar (`DLL load failed 126 = MOD_NOT_FOUND`). Pipeline caiu para OpenCL com mixed precision. Independente do bug do FAD, vale conferir `cudart64*.dll` em `C:\Users\luizd\miniforge3\envs\datascience\Library\bin`.

---

## Plano de retomada — duas opções

### Opção A (preferida): regerar o XML do FAD via OpenFF Toolkit

```powershell
conda activate datascience
pip install openff-toolkit openff-interchange ambertools
```

```python
# regen_fad_xml.py
from openff.toolkit import Molecule, ForceField
from openff.units import unit
from openff.interchange import Interchange

# SMILES canônico do FAD oxidado
SMILES_FAD = ("c1nc2c(n1[C@@H]3O[C@H](COP(=O)(O)OP(=O)(O)OC[C@@H](O)"
              "[C@@H](O)[C@@H](O)CN4c5cc(C)c(C)cc5N=C6C4=NC(=O)NC6=O)"
              "[C@@H](O)[C@H]3O)c(N)nc2")

mol = Molecule.from_smiles(SMILES_FAD)
mol.generate_conformers(n_conformers=1)
ff = ForceField("openff-2.1.0.offxml")
interchange = Interchange.from_smirnoff(force_field=ff, topology=[mol])
interchange.to_openmm_xml("cofactor_FAD_TcTR_v3.xml")
```

Substituir `cofactor_FAD_TcTR.xml` pelo `_v3` e re-rodar:
```powershell
python md_clean.py --fresh --ns 1 --tag prod_v3 --platform CUDA --allow-fallback
```

### Opção B: regerar via AmberTools antechamber

```bash
# Dentro de container ou WSL com ambertools instalado
antechamber -i fad.mol2 -fi mol2 -o fad.prep -fo prepi -c bcc -nc 0
parmchk2 -i fad.prep -f prepi -o fad.frcmod
# Converter para OpenMM XML via parmed
python -c "import parmed; ..."
```

### Fix do CUDA (paralelo)

```powershell
conda list openmm-cuda 2>$null    # confirmar versão instalada
where.exe cudart64_12.dll          # ver se está no PATH
# Se faltar:
conda install -c conda-forge openmm-cuda
```

---

## Arquivos importantes

| Caminho | Conteúdo |
|---|---|
| `md_clean.py` | v2 749 linhas — pipeline OK até minimização |
| `STATUS_REPORT_2026-05-13_v2.md` | Diagnóstico anterior (parent-anchored H placement) |
| `docking_project/figures/smoke_test_fad/cofactor_FAD_TcTR.xml` | **XML quebrado** — atom typing scrambled |
| `docking_project/figures/smoke_test_fad/cofactor_FAD_TcTR.pdb` | Template (heavy + Hs) — coordenadas OK, conectividade OK |
| `docking_project/figures/md_clean/minimized_dryrun_v2_prod.pdb` | PDB minimizado (último estado válido pré-NaN) |
| `md_run_v2_prod.log` | Log completo da última execução (1.8 KB no início, ~4 KB total com forensic + traceback) |

---

## ⏸ NÃO retomar a iteração MD em ciclos automáticos sem o fix do XML

Próxima execução só faz sentido depois que o XML do FAD for regenerado. Tentar `md_clean.py --fresh` repetidamente sem corrigir vai dar NaN sempre.
