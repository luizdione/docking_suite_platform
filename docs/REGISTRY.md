# Registry declarativo

Fonte única de verdade sobre **quem** é docado, **contra o quê**, e **sob que
protocolo**. Código lê o registry; código não contém dicionário de domínio.

Antes disso a plataforma tinha três dicionários de alvos que divergiam entre si
— um deles atribuía a *TcCruzain* o PDB `3GRS`, que é glutationa redutase
humana. O registry existe para que isso não possa voltar a acontecer.

## Hierarquia

```
registry/
  program.json                              guarda-chuva: identidade, defaults químicos
  campaigns/<ID>.json                       campanha (= subprojeto, = um paper)
  targets/<id>.json                         ALVO      — global, reutilizável
  ligands/<id>.json                         LIGANTE   — global, reutilizável
  pairs/<CAMPANHA>__<alvo>__<ligante>.json  PAR       — escopado por campanha
  protocols/<id>.json                       PROTOCOLO — docking, estrutura, MD
```

**ALVO e LIGANTE são globais; o PAR é da campanha.** Um mesmo hGR é contra-alvo
de seletividade num estudo e poderia ser alvo primário em outro: o alvo não
muda, o par muda.

Foi a ausência do conceito de par que permitiu `md_simulation.py` reusar
`cofactor_FAD_TcTR.xml` para alvos que não são flavoenzimas.

## O par é onde mora a peculiaridade

Tudo que só existe *na combinação* alvo × ligante, e que antes ficava implícito
no script, é campo obrigatório aqui:

| Campo | Por que é do par, não do alvo nem do ligante |
|---|---|
| `assembly_state.required` | O sítio da jasplaquinolida só existe na junção de 3 protômeros do filamento; o da latrunculina existe no monômero. Mesmo alvo, exigências opostas. |
| `experimental_kd[].state` | **Obrigatório.** Citocalasina D é 4,6 µM em G-actina e ~2 nM na extremidade barbada. Kd sem estado é erro de categoria. |
| `engines.affinity_estimator_reliable` | Um estimador pode ser confiável num sistema e inútil noutro. É configuração, não suposição. |
| `validation.training_leakage` | Se o PDB de referência é anterior ao cutoff de treino do preditor, o acerto não prova nada. |
| `validation.verdict_rule` | Ler só o `model_0` de um co-dobramento é estatística de ordem disfarçada de resposta. |
| `cofactor_requirement` | FAD na tripanotiona redutase e ATP/Mg na actina não são acessórios: sem eles o sítio deixa de existir. |

## Campos que mudam o comportamento do pipeline

**Alvo** — `structure.source` (`pdb` \| `alphafold` \| `local`). Modelo predito
exige `confidence.global` (pLDDT), e o registry recusa o arquivo sem isso.
`chemistry.cofactor_source` (`crystal` \| `transferred` \| `none`); `transferred`
exige `cofactor_transfer.{method,donor_pdb,align_scope}` — a superposição para
transferir cofator pode precisar ser **local ao bolso**, não global, quando o
alvo tem movimento interdomínio.

**Sítio** — `sites.<nome>.definition` aceita cinco modos. Três já eram
implementados por `_calculate_grid_center()` (`ligand_center`, `metal_center`,
`residue_center`); dois foram acrescentados para cobrir alvo sem co-cristal
próprio (`derived_from_reference_complex`, `explicit_center`). `box_policy`
diz se a caixa é fixa ou derivada do bounding box do ligante cristalográfico.

**Ligante** — `topology.is_macrocycle` obriga
`prep.pdbqt.rigid_macrocycles: false`: um anel de 19 membros tratado como rígido
produz pose errada por construção, sem erro. `prep.protonation.required` obriga
`post_check`, porque a protonação falha em silêncio.

## Uso

```python
import registry

registry.get_target("TcTR")                  # aceita id ou alias ("1BZL", "3GRS")
registry.resolve_site("TcTR")                # sítio default, com id resolvido
registry.get_pair("SUB01", "TcTR", "2-scn-nq")
registry.resolve_protocol("SUB01", "md", "TcTR", "2-scn-nq")   # par > campanha
registry.list_targets("SUB01")
registry.validate_all()                      # [] = tudo certo
```

`DOCKING_REGISTRY` sobrepõe o diretório; `DOCKING_CAMPAIGN` sobrepõe a campanha
corrente do dashboard.

## Regras de escrita

- **UTF-8 explícito, sempre.** Os campos `notes`/`notes_pt` carregam prosa
  científica acentuada e o default do Windows (cp1252) corrompe em silêncio.
- **Gerar, não transcrever.** Os alvos e ligantes vieram de
  `tools/export_target_db_to_registry.py` e
  `tools/export_compound_db_to_registry.py`. Transcrever 10 alvos × 12 campos à
  mão é como nasceram os PDBs errados.
- **Campo ausente fica ausente.** Nada de preencher com valor plausível.

## Validação

```bash
python registry.py                      # resumo + validação
python registry.py --campaign SUB01     # alvos da campanha
python tests/test_registry_parity.py    # gate: SUB01 inalterado
```

`validate_all()` é um validador próprio, sem dependência dura de `jsonschema`
(que não está em `requirements.txt`). Se `jsonschema` estiver instalado e houver
arquivos em `registry/_schema/`, eles são usados como checagem adicional.

## Ponte com o código legado

`registry.as_legacy_target_db(campaign)` projeta o registry na forma do antigo
`docking_v2.TARGET_DB`, com os mesmos tipos — inclusive `grid_size` como
**tupla**, porque uma lista desempacota igual e a regressão passaria
despercebida em runtime. `docking_v2.load_target_db()` a consome e cai num
fallback embutido, ruidosamente, se o registry sumir.
