# Daily Briefing (UPDATE da noite) — Docking / MD de TcTR + 2-SCN-NQ

**Data:** 2026-05-23 (~22:30, run automático do `daily-briefing_docking-e-dm`)
**Foco:** módulo de MD (`md_clean.py`, sucessor do `md_module.py`) e force-field do FAD
**Supersede:** `BRIEFING_2026-05-23_md_module.md` (gerado de manhã, 11:57) — **houve atividade real depois dele**
**Delta desde o briefing da manhã:** **grande** — o gate que travava a MD foi removido às 15:05 e o debug avançou até ~18:00

---

## TL;DR

O briefing da manhã dizia que a MD estava "a um único comando de distância": faltava rodar `python regen_fad_xml.py`. **Esse comando foi rodado hoje (15:05)** e o gate caiu: `cofactor_FAD_TcTR_v3.xml` **agora existe** e — verificado por mim na cópia autoritativa do Windows — o par v3 (PDB+XML) está **internamente consistente**:

- **91 bonds no PDB CONECT == 91 bonds no template do XML** (0 divergências nos dois sentidos);
- **248 torsions** no v3 (229 proper + 19 improper) **vs. 0** no XML antigo — o `PeriodicTorsionForce` vazio foi corrigido;
- 86 atom types, 4 forças presentes (Bond/Angle/Torsion/Nonbonded);
- `regen_fad_xml.py` foi **atualizado hoje** para a rota de cargas **NAGL** (`SCRIPT_VERSION = "vNAGL-2026-05-23"`, era `vISO-2026-05-22`).

Ou seja: o defeito original do FAD (scramble de elementos + torsions vazias) está, pela análise estática, **resolvido**. A hipótese de "bonds do PDB não batem com o XML" — exatamente o que o seu `diag_fad_bonds.py` foi escrito para caçar — está **descartada no nível estático PDB↔XML**.

**O que ainda NÃO sei (e não dá para saber daqui):** se o `createSystem` no sistema completo (183k átomos) passou e se o warm-up parou de dar NaN. **Nenhuma saída de produção da MD foi gravada hoje** (o diretório `docking_project/figures/md_clean/` não muda desde 14/05; nenhum log novo do run das ~18:00). O `diag_fad_bonds.py` imprime no stdout, não em arquivo, então o resultado dele não ficou registrado.

> ⚠️ **Armadilha de sync (de novo, e eu caí nela antes de corrigir).** O mount do sandbox via OneDrive mostrou **cópias truncadas** de `v3.xml`, `v3.pdb` e `md_clean.py` (justamente os arquivos escritos/tocados mais recentemente). Uma análise crua no sandbox "concluiu" um bug de bonds inexistente (65 em vez de 91). A verificação na **visão autoritativa do Windows** desmente isso. Detalhe na §4.

---

## 1. Linha do tempo de hoje (2026-05-23)

| Hora (local aprox.) | Evento | Evidência |
|---|---|---|
| 11:51 | `md_clean.py` reescrito (sessão da manhã, ver §5 do briefing das 11:57) | mtime |
| 11:57 | Briefing da manhã gerado | `BRIEFING_2026-05-23_md_module.md` |
| **15:05** | **`regen_fad_xml.py` rodado → `cofactor_FAD_TcTR_v3.xml` + `_v3.pdb` gravados** | mtime 15:05 dos dois arquivos; gate da manhã removido |
| 17:56 | `diag_fad_bonds.py` criado (localiza onde os bonds do FAD divergem entre PDBFile / build_canonical / modeller.add / XML) | novo arquivo |
| 17:57–18:00 | `md_clean` importado/compilado em **dois** interpretadores (`.cpython-311` e `.cpython-310`) → diag e/ou tentativa de run | `__pycache__` |
| (depois) | **Nenhuma saída de MD** (sem `prod.log` novo, sem `minimized_*` novo, sem DCD) | `find` pós-18:00 vazio |

A `regen` agora usa **NAGL** para cargas — isso responde à pendência nº 2 do briefing da manhã (AM1-BCC vs Gasteiger): NAGL (GNN da OpenFF) é uma escolha moderna e determinística, adequada para produção. Confirme no log da `regen` qual backend de carga foi de fato usado (e a carga formal total — ver §5).

---

## 2. O que EU verifiquei hoje (cópia autoritativa do Windows)

Tudo abaixo foi conferido lendo os arquivos pela visão do Windows (não pelo mount do sandbox, que estava defasado):

| Verificação | Resultado | Significado |
|---|---|---|
| `cofactor_FAD_TcTR_v3.pdb` íntegro | ✅ 86 átomos, CONECT para todos os 86 (inclusive cauda fosfato/ribitol), termina em `END` | PDB completo no disco |
| Bonds PDB CONECT × bonds do XML (`<Residue>`) | ✅ **91 == 91, 0 faltando, 0 sobrando** | topologia do FAD consistente; matching de template não deve falhar por bonds |
| Torsions no v3.xml | ✅ **248** (229 proper + 19 improper) | corrige o `PeriodicTorsionForce` vazio do XML antigo (0) |
| Atom types / forças no v3.xml | ✅ 86 types; Bond+Angle+Torsion+Nonbonded presentes | force-field completo |
| `regen_fad_xml.py` e `diag_fad_bonds.py` compilam | ✅ `py_compile` OK | sintaxe íntegra |

**Conclusão da §2:** a regeneração do FAD (v3) produziu um par PDB+XML quimicamente coerente e auto-consistente. O bloqueio diagnosticado em 05-21 (force-field parametrizado para a identidade atômica errada + torsions vazias) está, estaticamente, **corrigido**.

---

## 3. O que continua em aberto (não verificável daqui)

1. **`createSystem` no sistema completo (183k átomos):** o `md_clean.py` tem, desde hoje, um bloco `DIAG 2026-05-23` (linhas ~340–364) que compara bonds FAD↔XML, checa nomes duplicados e tenta `createSystem` PRE-solvente. **Não há log desse run** — só dá para confirmar rodando (ver §6).
2. **Warm-up NaN:** o último run real de MD (`md_run_v2_prod.log`, 14/05) parou em `Particle coordinate is NaN`. Como nenhuma MD foi rodada hoje até gravar saída, **não sei se o NaN sumiu** com o par v3. O forensic de 14/05 acusava `BOND ANOMALA 2.06 Å: N3A–C4A` (geometria) — isso é um eixo **diferente** do bond-matching (que agora está OK) e pode ou não persistir; só o run dirá.
3. **CUDA do OpenMM** (carry-over): plugins falham com `DLL load failed 126`; cai para OpenCL. Independe do FAD.

> **Importante:** o "stable 5 ns" que aparece em `md_progress_TcTR.log` é de **11/05** (caminho FAST, 177k átomos) — **não** é um run de hoje. Não confunda com sucesso da pipeline `--fresh` v3.

---

## 4. A armadilha de sync OneDrive↔sandbox (transparência)

Antes de ler pela visão do Windows, o shell do sandbox mostrou os três arquivos recém-escritos **truncados**, e uma análise estática crua tirou conclusão errada:

| Arquivo | Visto pelo mount (defasado) | Real no Windows | Veredito |
|---|---|---|---|
| `cofactor_FAD_TcTR_v3.pdb` | cortado em `CONECT 38 31`, sem `END`; só **65 bonds** | 86 átomos, **91 bonds**, `END` presente | mount **stale** |
| `cofactor_FAD_TcTR_v3.xml` | `unclosed token` no fim (`</NonbondedForce></F`) | seções Residue/Torsions íntegras (91 bonds, 248 torsions) | mount **stale** (só a cauda) |
| `md_clean.py` | "SyntaxError linha 936: `(` never closed" | compila no disco | mount **stale** |

A prova de que `md_clean.py` está íntegro no disco: os `.pyc` `md_clean.cpython-310/311` foram gravados às 17:57–18:00 de hoje — **impossível** se a fonte tivesse erro de sintaxe. O canário do próprio projeto (`SCRIPT_VERSION` "se nao aparecer no log, OneDrive nao sincronizou") aponta o mesmo. **Lição operacional:** análises automáticas no sandbox devem validar código/arquivos pela cópia do Windows (ou esperar o OneDrive hidratar), nunca pelo mount cru logo após uma escrita.

**Nenhum arquivo foi escrito ou modificado por mim nesta sessão, exceto este briefing.** Não toquei em código, XML ou PDB.

---

## 5. Pontos de atenção para o run de validação

- **Carga formal do FAD:** NAGL/AM1-BCC dão cargas parciais que **somam à carga formal**. O FAD oxidado a pH 7,4 costuma ser **−2** (difosfato desprotonado); o template tem 86 átomos `C27H33N9O15P2`. Confira no log da `regen` a fórmula/carga total validada — se divergir do esperado, ajustar antes de produção.
- **Nomes de átomo duplicados:** o bloco DIAG do `md_clean` também checa duplicados dentro do resíduo FAD (quebram o matching). No PDB v3 os 86 nomes são únicos — ok no isolado; o risco real é só **se** `modeller.add` colidir nomes entre os dois FADs/LIG (improvável, mas o DIAG vai dizer).
- **Se o `createSystem` AINDA falhar** mesmo com bonds consistentes (91=91): suspeite, nesta ordem, de (a) `modeller.add` perdendo/criando bond na junção (o `diag_fad_bonds.py` mede isso em `[B vs C]`), (b) o `lig_TcTR_smirnoff.xml`, ou (c) leitura de um arquivo ainda **não sincronizado** pelo OneDrive (rode só depois de confirmar o `SCRIPT_VERSION` no log).

---

## 6. Próximos passos priorizados (no env `datascience`, no Windows)

1. **Rodar o diagnóstico e LER a saída** (responde a §3.1):
   ```powershell
   conda activate datascience
   cd "C:\Users\luizd\OneDrive\Data_Science\Projetos\docking"
   python diag_fad_bonds.py
   ```
   Olhe especialmente `[C vs X]` (deve vir **vazio**, como minha análise estática prevê) e a última linha `createSystem (topologia C, FF completa): OK/FALHOU`.
2. **Se o createSystem passar, rodar a MD curta** e conferir warm-up:
   ```powershell
   python md_clean.py --fresh --ns 1 --tag prod_v3 --platform CUDA --allow-fallback
   ```
   Confirme no log: `SCRIPT_VERSION` correto, **sem** `BOND ANOMALA`, **sem** NaN no warm-up, e o backend de carga usado.
3. **Em paralelo:** consertar o CUDA (`conda install -c conda-forge openmm-cuda`; conferir `cudart64_*.dll` em `…\envs\datascience\Library\bin`).
4. Se a curta (1 ns) passar: escalar 5/50 ns e rodar `analyze()` (RMSD-no-bolso + RMSF + convergência).
5. **(Higiene)** `git init` na pasta para versionar a v2.1.0; apagar `md_clean.py.truncated_2026-05-23.bak` (resíduo da sessão da manhã, **não** é versão válida).

> ⚠️ Mantido de `NEXT_SESSION.md`: não reiniciar ciclos automáticos de MD sem confirmar o run de validação acima.

---

## 7. Arquivos-chave (estado hoje)

| Caminho | Papel | Estado |
|---|---|---|
| `…/smoke_test_fad/cofactor_FAD_TcTR_v3.xml` | Force-field regenerado do FAD | ✅ **gerado hoje 15:05**; 91 bonds, 248 torsions, 86 types |
| `…/smoke_test_fad/cofactor_FAD_TcTR_v3.pdb` | Template canônico (86 átomos) | ✅ íntegro; CONECT == XML (91=91) |
| `regen_fad_xml.py` | Gera o FF do FAD | ✅ atualizado p/ NAGL (`vNAGL-2026-05-23`); executado 15:05 |
| `diag_fad_bonds.py` | Localiza divergência de bonds FAD (PDBFile/canonical/modeller/XML) | 🆕 criado hoje 17:56; **rodar e ler saída** |
| `md_clean.py` | Pipeline MD ativo (com bloco DIAG 2026-05-23, ~L340–364) | ✅ íntegro no disco (mount mostra stale) |
| `…/smoke_test_fad/cofactor_FAD_TcTR.xml/.pdb` | Par antigo (FF quebrado) | intocado; substituído pelo v3 |
| `…/cofactor_FAD_TcTR_v3_elemfix.xml(.INVALID.txt)` | Tentativa descartada | ⛔ não usar |
| `docking_project/figures/md_clean/` | Saídas de MD | sem mudança desde 14/05 (nenhum run novo) |
| `md_run_v2_prod.log` / `md_progress_TcTR.log` | Últimos runs reais | 14/05 (NaN) / 11/05 (5 ns FAST, **antigo**) |

---

*Briefing automático. Verificações: leitura autoritativa via Windows (Read), comparação estática PDB-CONECT × XML-Residue (91=91, 0 divergências), contagem de torsions/types, `py_compile`, e varredura de mtimes. **Nenhuma simulação foi executada; nenhum arquivo de código/FAD foi modificado — apenas este briefing foi escrito.** O outcome do `createSystem`/warm-up no sistema completo permanece desconhecido porque o run das ~18:00 não gravou log; rodar `diag_fad_bonds.py` (§6.1) resolve isso.*
