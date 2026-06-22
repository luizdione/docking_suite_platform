# FAD XML scramble diagnosis (cofactor_FAD_TcTR.xml)

**Diagnosticado em:** 2026-05-14
**Status:** confirmado — XML é fundamentalmente inválido como force-field para FAD

## Localização do arquivo

`docking_project/figures/smoke_test_fad/cofactor_FAD_TcTR.xml`

## Sumário dos defeitos

| Defeito | Magnitude | Impacto |
|---|---|---|
| Atom-type element mismatch | 33 dos 86 átomos pesados (38 %) | Massa atômica errada, parâmetros LJ/NB errados |
| HarmonicBondForce vazio | **0 bonds** definidos por classe | Sem termo harmônico para nenhuma ligação |
| Bonds da Residue (topologia) | 91 bonds presentes | Conectividade definida apenas em nível topológico |

## Detalhes — elementos errados (33 atom types)

Padrão: o element column do PDB original e os atributos `<Type element="...">` no XML estão dessincronizados em relação aos nomes de átomos. 9 dos 10 átomos do anel adenina têm element errado.

Saída de `_fad_diagnose.py`:

```
PDB wrong-element atoms: 33/86
  PA    declared=C, should be=P
  O1A   declared=C, should be=O
  O2A   declared=C, should be=O
  O5B   declared=C, should be=O
  C5B   declared=N, should be=C
  ...

XML wrong-element types: 33/86
```

Consequência imediata para a dinâmica: `<Type name="FAD20" element="O" mass="15.9994">` é atribuído a N3A (que deveria ter element=N, mass=14.0067). OpenMM cria a partícula como oxigênio com massa errada e parâmetros LJ de oxigênio (σ ≈ 0.30 nm) em vez de nitrogênio.

## Detalhes — bonds (91 missing, todos)

```
Residue bonds defined: 91
HBF bond classes defined: 0
```

A seção `<Residue name="FAD">` lista 91 bonds via `<Bond atomName1="..." atomName2="...">`. A seção `<HarmonicBondForce>` está **vazia** — não tem nenhum `<Bond class1="FADn" class2="FADm" length="..." k="...">`.

Como OpenMM monta o sistema:
1. A topologia da Residue define as conexões físicas entre átomos
2. O ForceField busca parâmetros harmônicos por par de classes (FADn-FADm)
3. Como nenhum par está em HBF, **nenhum termo harmônico é adicionado**
4. As ligações ficam mantidas apenas pela 1-3 e 1-4 das angle/torsion (também provavelmente vazias) e pelas forças não-bonded
5. Em minimização (otimização sem inércia), as não-bonded podem ainda convergir a uma geometria razoável
6. Em dinâmica a 5 K com `dt=0.25 fs`, qualquer pequeno gradiente sem mola harmônica para conter leva a deslocamento sem retorno → NaN

A combinação **"residue declarada com 91 bonds, força sem termo harmônico"** explica:
- Por que `forensic` aponta `BOND ANOMALA 2.06 A: N3A-C4A`: a "bond" topológica existe mas tem k=0
- Por que a minimização converge (-3 MJ/mol é uma energia razoável quando dominada por LJ/Coulomb)
- Por que warm-up explode em NaN imediatamente: sem termo harmônico, qualquer dt > 0 leva drift

## Origem do bug

O XML provavelmente foi gerado por uma ferramenta automática (ex.: PDB → XML ad-hoc, ou OpenMM Modeller + custom script) que:
1. Leu o element column do PDB (já com dessincronização)
2. Criou `<Type>` enumerados sequencialmente sem corrigir o element
3. Definiu `<Residue>` com a topologia certa, mas
4. **Não foi rodada a etapa que popula `<HarmonicBondForce>` com parâmetros** (provavelmente faltou rodar antechamber/parmchk2 ou OpenFF Toolkit + Interchange para extrair os parâmetros do force-field selecionado)

## Como corrigir

Existem três caminhos viáveis. Em ordem de robustez:

### Caminho A — Regenerar via OpenFF Toolkit (recomendado)

OpenFF tem o force-field `openff-2.1.0.offxml` (SMIRNOFF) e roda em pure Python.

```python
# scripts/regen_fad_xml_openff.py
from openff.toolkit import Molecule, ForceField
from openff.interchange import Interchange

# SMILES do FAD oxidado (sem o radical da forma reduzida)
SMILES = ("c1nc2c(N)nc3n(cnc3c2NC(=O)NC1=O)[C@@H]4O[C@H]"
          "(COP(=O)([O-])OP(=O)([O-])OC[C@@H](O)[C@@H](O)"
          "[C@@H](O)CN5c6cc(C)c(C)cc6N=C7C5=NC(=O)NC7=O)"
          "[C@@H](O)[C@H]4O")

mol = Molecule.from_smiles(SMILES, allow_undefined_stereo=True)
mol.generate_conformers(n_conformers=1)

ff = ForceField("openff-2.1.0.offxml")
inter = Interchange.from_smirnoff(force_field=ff, topology=mol.to_topology())

# Salvar como OpenMM System (não como ForceField XML)
inter.to_openmm_system_xml("cofactor_FAD_TcTR_v3.system.xml")
```

**Limitação:** OpenFF gera um *OpenMM System* (já com parâmetros aplicados a uma topologia específica), não um *ForceField* (template reusável). Para usar em `md_clean.py`, precisaria mudar o pipeline para carregar o sistema serializado em vez de aplicar ForceField em runtime.

### Caminho B — Antechamber + parmchk2 (AmberTools)

```bash
conda install -c conda-forge ambertools
# Converte FAD.mol2/sdf → AMBER prmtop com GAFF2
antechamber -i fad.mol2 -fi mol2 -o fad.prep -fo prepi -c bcc -nc -2 -at gaff2
parmchk2 -i fad.prep -f prepi -o fad.frcmod -s gaff2
# Converter para OpenMM XML via parmed
python -c "
import parmed as pmd
amber = pmd.amber.AmberParm('fad.prmtop')
amber.write_xml('cofactor_FAD_TcTR_v3.xml')
"
```

### Caminho C — Patch manual (rápido mas frágil)

Apenas para validação rápida:
1. Corrigir `<Type element>` baseado no primeiro caractere do atom name
2. Adicionar bonds manuais em HBF com valores típicos (aromatic C-N: 0.134 nm, k=290000 kJ/mol/nm²)

Não é fix definitivo — falta também angle/torsion/charges corretos.

## Próximos passos para o usuário

1. Decidir entre A/B/C
2. Após escolher, gerar `cofactor_FAD_TcTR_v3.xml` (ou system.xml para Caminho A)
3. Substituir o XML usado no `md_clean.py` (linha que carrega `cofactor_FAD_TcTR.xml`)
4. Re-rodar: `python md_clean.py --fresh --ns 1 --tag prod_v3 --platform CUDA --allow-fallback`
5. Verificar que minimização não produz `BOND ANOMALA` na seção forensic
6. Se OK, escalar para 5 ns / 50 ns

## Arquivos relacionados

- `_fad_diagnose.py` — script de diagnóstico que gerou esta análise
- `docking_project/figures/smoke_test_fad/cofactor_FAD_TcTR.pdb` — PDB template (coordenadas OK, element column scrambled)
- `docking_project/figures/smoke_test_fad/cofactor_FAD_TcTR.xml` — XML quebrado (esta análise)
- `md_clean.py` — pipeline OpenMM (linha 73-83 carrega o XML)
- `STATUS_REPORT_2026-05-13_v2.md` — análise anterior (parent-anchored H placement, antes do diagnóstico do XML)
