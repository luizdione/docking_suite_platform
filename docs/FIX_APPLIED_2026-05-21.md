# Correções aplicadas — FAD force-field / MD (2026-05-21)

Contexto e diagnóstico completo: `BRIEFING_2026-05-21_md_module.md`.

Os dois defeitos confirmados do `cofactor_FAD_TcTR.xml`:
1. 33/86 atom types com **elemento/massa/LJ do átomo errado** (scramble) — e a coluna `element` do PDB-template estava embaralhada do mesmo jeito;
2. `PeriodicTorsionForce` **vazio** (0 torsions).

Causa-raiz da geração: a molécula foi parametrizada numa **ordem de átomos diferente** da ordem dos nomes canônicos, e o FAD passou por um caminho de fallback que largava torsions (a rota SMIRNOFF que já funciona para o ligante, `_build_lig_xml_smirnoff`, emite torsions e elementos corretos).

---

## O que foi feito agora (verificável, não-destrutivo)

| # | Mudança | Arquivo | Status |
|---|---|---|---|
| 1 | Coluna `element` corrigida (determinístico, 33 campos) | **novo** `docking_project/figures/smoke_test_fad/cofactor_FAD_TcTR_v3.pdb` | ✅ feito e verificado (0 mismatches) |
| 2 | Script de regeneração do force-field (rota SMIRNOFF, ordem canônica garantida, auto-validação) | **novo** `regen_fad_xml.py` (raiz) | ✅ escrito; sintaxe OK; leitura de PDB/bonds testada (86 átomos / 91 bonds) |
| 3 | `fix_fad_xml.py` deixou de ser perigoso: regex agora lê `type1/type2` (lia 0 das 91 bonds e injetava 91 **duplicatas**); injeção de bonds desativada quando o HBF já tem bonds; banner de stopgap | `auto_fix/fix_fad_xml.py` | ✅ feito; dry-run agora reporta `HBF bonds (existing): 91` e 0 injetadas |
| 4 | Pipeline passa a **preferir** o par v3 (XML+PDB) quando existir, com fallback consistente para o antigo | `md_clean.py` (resolução de `fad_xml` e `canonical_template`) | ✅ feito; `py_compile` OK |

Backups criados: `cofactor_FAD_TcTR.pdb.bak_2026-05-21`, `md_clean.py.bak_2026-05-21`.
Originais `cofactor_FAD_TcTR.xml` e `cofactor_FAD_TcTR.pdb` **não foram modificados**.

> Segurança: o `md_clean.py` só usa o `cofactor_FAD_TcTR_v3.pdb` **quando o `cofactor_FAD_TcTR_v3.xml` também existir**. Como o v3 XML ainda não foi gerado, o pipeline continua usando o par antigo (estado atual preservado, sem regressão). Assim que você rodar o passo abaixo, ele troca automaticamente para o par v3.

---

## Passo que falta (precisa do seu ambiente `datascience`)

A regeneração do force-field exige RDKit + OpenFF, que **não couberam no sandbox** (sem espaço em disco) — e, além disso, a única validação que importa (o sistema solvatado de 183k átomos) só roda no seu ambiente. Então:

```powershell
conda activate datascience
cd "C:\Users\luizd\OneDrive\Data_Science\Projetos\docking"
python regen_fad_xml.py
```

Isso gera `docking_project/figures/smoke_test_fad/cofactor_FAD_TcTR_v3.xml` com elementos/massas/LJ/cargas/bonds/angles/**torsions** corretos e nomes de átomo canônicos, e faz uma auto-validação (createSystem + warm-up do FAD isolado). Em seguida:

```powershell
python md_clean.py --fresh --ns 1 --tag prod_v3 --platform CUDA --allow-fallback
```

e confira que o forensic **não** acusa `BOND ANOMALA` e que o warm-up **não** dá NaN.

### Observações importantes

- **Cargas:** o `regen_fad_xml.py` usa **AM1-BCC** se houver backend (AmberTools `sqm` / OpenEye); senão cai para **Gasteiger** (estabiliza, mas para produção/publicação prefira AM1-BCC — instale AmberTools no env). O log diz qual foi usado.
- **Protonação:** o template tem 86 átomos = FAD neutro `C27H33N9O15P2` (difosfato protonado). O script valida a fórmula/carga e avisa se divergir.
- **SMILES de fallback:** o script tenta perceber ordens de ligação pela geometria 3D (`rdDetermineBonds`); se falhar, usa um SMILES embutido do FAD. Se a auto-validação reclamar da fórmula, ajuste `FAD_SMILES` no script.
- **CUDA:** continua quebrado (plugins OpenMM `DLL load failed 126`); independente do FAD. `conda install -c conda-forge openmm-cuda` e conferir `cudart64_*.dll` no `Library\bin` do env.

---

## Arquivo descartado

`cofactor_FAD_TcTR_v3_elemfix.xml` (e `.INVALID.txt`): tentativa de "fix só de elementos", **reprovada** (quebra template matching e não corrige LJ). Não usar — substituída pelo `regen_fad_xml.py`.
