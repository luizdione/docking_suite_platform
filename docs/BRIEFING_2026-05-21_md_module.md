# Daily Briefing — Docking / MD de TcTR + 2-SCN-NQ

**Data:** 2026-05-21 (run automático do `daily-briefing_docking-e-dm`)
**Foco:** módulo de MD (`md_clean.py`, sucessor do `md_module.py`)
**Última atividade real no projeto:** 2026-05-14 — **7 dias parado**
**Mudança desde o último briefing (05-20):** nenhuma no projeto (só este briefing foi adicionado)

---

## TL;DR

O pipeline de MD continua **parado no mesmo ponto desde 14/05**: minimiza, mas explode em NaN no warm-up. Nada foi executado desde então — o `auto_fix/` segue sem nunca ter rodado e o XML do FAD nunca foi regenerado.

A verificação independente de hoje **confirma o diagnóstico de 05-20 e o aprofunda em três pontos novos e decisivos**:

1. **Não é só "elemento/massa trocados". O force-field inteiro foi parametrizado para os átomos errados.** Os raios de Lennard-Jones (σ em `NonbondedForce`) batem com o elemento **errado** em **86/86** átomos. Ou seja: massa **e** LJ **e** (por consequência) cargas foram gerados para a identidade atômica embaralhada. **Nenhum patch de rótulo conserta isso** — só a regeneração completa do force-field resolve.
2. **O PDB-template está embaralhado de forma idêntica ao XML** (33/86 na coluna `element`). Eles são *mutuamente consistentes no erro* — é por isso que o pipeline casa templates e minimiza apesar de tudo. Corrigir só o XML **quebra** o template matching (testado: `No template found for residue FAD`). Qualquer correção tem de tocar **os dois** (PDB + XML) ou regenerar ambos.
3. **O `auto_fix/fix_fad_xml.py` é ativamente perigoso como está.** O parser dele procura bonds no formato `class1/class2`, mas o XML usa `type1/type2`. Resultado: ele "vê" 0 das 91 bonds existentes e **injeta 91 bonds duplicadas** → toda ligação com mola dobrada. Rodar `fix_fad_xml.py --inplace` deixaria o XML **pior**, não melhor. **Não rodar.**

**Caminho único viável:** regenerar o force-field do FAD a partir dos *nomes* de átomo (que estão corretos), via OpenFF Toolkit ou antechamber/GAFF2, gerando elementos, massas, LJ, cargas, bonds, angles **e** torsions de forma consistente — e remapeando os nomes de átomo para casar com `cofactor_FAD_TcTR.pdb`.

---

## 1. Estado do pipeline (inalterado desde 14/05)

`md_clean.py --fresh` (858 linhas, módulo ativo; `md_module.py` de 3.402 linhas aposentado em `obsoletos/`) reconstrói o sistema completo com sucesso até a minimização (log `md_run_v2_prod.log`, 14/05):

| Etapa | Resultado |
|---|---|
| Reconstrução | PDBFixer → addHydrogens(proteína) → 2× FAD canônicos via Kabsch → LIG → TIP3P + 0,15 M NaCl |
| Sistema final | **183.179 átomos** (TcTR dímero + 2 FAD + 1 LIG + 55.906 HOH + 311 Na/Cl) |
| Minimização | converge, PotentialEnergy ≈ **−3,0 × 10⁶ kJ/mol**, 0 NaN |
| Último estado válido | `docking_project/figures/md_clean/minimized_dryrun_v2_prod.pdb` |
| Warm-up (dt=0,25 fs @ 5 K) | **NaN imediato** |

Forensic pré-warm-up: `BOND ANOMALA 2.06 Å N3A–C4A` e força máx ~7,8×10² kJ/mol/nm no anel adenina do FAD.

---

## 2. Diagnóstico do FAD XML — verificado hoje (parse read-only independente)

Arquivo: `docking_project/figures/smoke_test_fad/cofactor_FAD_TcTR.xml` (mtime 12/05, intocado).

| Item | 14/05 (`FAD_XML_SCRAMBLE`) | 20/05 (briefing) | **21/05 (hoje)** | Veredito |
|---|---|---|---|---|
| Atom types totais | 86 | 86 | **86** | — |
| Element/massa errados | 33/86 | 33/86 | **33/86** (PA→C, N3A→O, C4A→P…) | ✅ bug real |
| `HarmonicBondForce` | "vazio (0 bonds)" | 91 bonds, k>0 | **91 bonds, k=1,6×10⁵–6,8×10⁵, nenhuma k=0** | ❌ 14/05 estava errado |
| Bond N3A–C4A | "não existe" | existe | **existe** | ❌ 14/05 estava errado |
| `HarmonicAngleForce` | n/a | 157 angles | **157 angles, k válido** | ✅ OK |
| `PeriodicTorsionForce` | n/a | 0 torsions | **0 torsions (vazio)** | ⚠️ defeito real |
| `NonbondedForce` cargas | n/a | soma = 0 | **86 átomos, soma de cargas = 0,000** | neutro (ver §3) |
| **LJ σ vs elemento (NOVO)** | — | — | **σ casa com elemento ERRADO em 86/86** | 🔴 force-field todo errado |

### 2.1 Por que isto fecha o caso (descoberta de hoje)

Os σ de Lennard-Jones vivem no `<NonbondedForce>`, **independentes do atributo `element`**. Testei a qual elemento cada σ corresponde:

- Concordância com o elemento **correto** (implícito pelo nome): **53/86**
- Concordância com o elemento **errado** (declarado no XML): **86/86**

Exemplos: `N3A` (deveria ser N, σ≈0,325) tem σ=0,3033 nm = oxigênio. `C4A` (C, σ≈0,339) tem σ=0,3742 = fósforo. `PA` (P, σ≈0,374) tem σ=0,3388 = carbono.

**Conclusão:** o gerador do XML aplicou massa, LJ e cargas usando a identidade atômica **errada** de ponta a ponta. Não é um rótulo solto — é o force-field inteiro parametrizado para os átomos trocados. Corrigir só `<Type element=…>` e a massa **deixa o LJ (e as cargas) ainda errados** → o sistema continuaria fisicamente inválido.

### 2.2 PDB e XML embaralhados juntos (descoberta de hoje)

A coluna `element` (cols 77-78) de `cofactor_FAD_TcTR.pdb` tem **o mesmo embaralhamento 33/86**, átomo por átomo, que o XML. Por isso o OpenMM consegue casar os templates (os dois estão errados igual) e a minimização roda.

Teste empírico (OpenMM 8.5.1, sandbox): gerei um candidato corrigindo só os elementos do XML e tentei `createSystem` no template `cofactor_FAD_TcTR.pdb` →

```
ValueError: No template found for residue 0 (FAD). The set of atoms matches FAD,
but the residue is missing 6 H-N bonds, 2 H-O bonds, 3 H-P bonds … too many.
```

Ou seja: mexer só no XML **dessincroniza** a dupla e quebra o matching. Qualquer correção precisa ser **consistente nos dois arquivos** (ou regenerar os dois).

### 2.3 O NaN é fenômeno do sistema solvatado, não do FAD isolado (descoberta de hoje)

Rodei o warm-up do pipeline (5 K, dt=0,25 fs, Langevin) no **FAD sozinho em vácuo**, com o XML original quebrado: **400 passos sem NaN**. Isso indica que o NaN de produção **não é reproduzível no FAD isolado** — ele emerge no sistema completo solvatado (183k átomos), onde massas erradas geram velocidades erradas no warm-up (v ∝ 1/√m) e os σ errados criam clashes com água/proteína. Implicação prática: **não dá para depurar isso no FAD isolado**; o teste de validação tem de ser o sistema inteiro, após a regeneração.

---

## 3. Pendências e bloqueios

- **`auto_fix/` nunca rodou** (criado 14/05 21:4x; sem `auto_md_report.md`/`auto_md_state.json`/`auto_md_logs/`).
- **🔴 `fix_fad_xml.py` não deve ser rodado.** Bug confirmado: o regex de leitura do HBF espera `class1/class2`, mas o XML grava `type1/type2`. Ele lê 0 bonds existentes e injetaria **91 bonds duplicadas** (toda mola dobrada). Além disso, mesmo a parte de elementos que ele acerta é **insuficiente** (não corrige LJ/cargas, e desincroniza do PDB — ver §2.1/2.2). Se quiser reaproveitá-lo, é preciso: (a) corrigir o regex para `type1/type2`, (b) **desligar** a injeção de bonds, e (c) aplicar a mesma correção de elemento ao PDB. Ainda assim, **não resolve LJ/cargas/torsions**.
- **Cargas neutras (soma = 0):** o FAD oxidado com dois fosfatos costuma ter carga líquida **−2** (difosfato desprotonado a pH 7,4). O XML traz soma 0 (fosfatos protonados). Não é necessariamente um bug, mas é uma **decisão de protonação a confirmar** na regeneração (`-nc 0` vs `-nc -2`).
- **CUDA do OpenMM quebrado** (carry-over de 05-14/05-20): 4 plugins CUDA falham (`DLL load failed 126 = MOD_NOT_FOUND`); pipeline caiu para OpenCL mixed precision. Independente do FAD, limita escala 0,1 ns → 5/50 ns.

---

## 4. Próximos passos recomendados (priorizados)

1. **Regenerar o force-field do FAD a partir dos nomes de átomo (único caminho válido).** Os *nomes* (PA, O1A, N3A, C4A…) estão corretos — eles definem a química certa. Pipeline sugerido:
   1. Corrigir a coluna `element` de `cofactor_FAD_TcTR.pdb` deterministicamente (primeira letra alfabética do nome). Isso dá uma estrutura com elementos corretos **e** nomes preservados.
   2. Gerar parâmetros com **OpenFF Toolkit** (`openff-2.1.0.offxml`) ou **antechamber + parmchk2** (GAFF2, AM1-BCC) a partir dessa estrutura corrigida.
   3. **Remapear os nomes de átomo** do force-field gerado de volta para os 86 nomes canônicos do `cofactor_FAD_TcTR.pdb` (por isomorfismo de grafo / coordenadas), para o template continuar casando no `md_clean.py`.
   4. Alternativa que evita o remap manual: usar o **`SMIRNOFFTemplateGenerator`** do `openmmforcefields` (que o `md_clean.py` já cita nos créditos, linha 9) alimentado com a `Molecule` OpenFF — ele registra o template de resíduo automaticamente.
   - Ponto de integração no `md_clean.py`: a linha **258-259** (`ff = app.ForceField("amber14-all.xml", "amber14/tip3p.xml", str(fad_xml), str(lig_xml))`) é onde o XML do FAD entra; o `createSystem` final é na **linha 347**. O FAD é construído via Kabsch a partir de `cofactor_FAD_TcTR.pdb` (linhas 292/303).
2. **Decidir a protonação do difosfato** (`-nc -2` recomendado para pH 7,4) **antes** de gerar cargas.
3. **Validar no sistema completo, não no FAD isolado** (ver §2.3): `python md_clean.py --fresh --ns 1 --tag prod_v3 --platform CUDA --allow-fallback` e conferir que o forensic **não** acusa `BOND ANOMALA` e que o warm-up não dá NaN.
4. **Consertar o CUDA** em paralelo (`conda install -c conda-forge openmm-cuda`; conferir `cudart64_*.dll` no `Library\bin` do env `datascience`).

> ⚠️ Mantido de `NEXT_SESSION.md`: **não** reiniciar ciclos automáticos de MD sem o force-field do FAD regenerado — vai dar NaN sempre. E **não** rodar `auto_fix/fix_fad_xml.py` como está (§3).

---

## 5. Arquivos-chave

| Caminho | Papel | Estado |
|---|---|---|
| `md_clean.py` (858 linhas) | Pipeline MD ativo | OK até minimização; FF do FAD entra na linha 258-259 |
| `obsoletos/md_module.py` (3.402 linhas) | Módulo MD antigo, monolítico | Aposentado |
| `docking_project/figures/smoke_test_fad/cofactor_FAD_TcTR.xml` | Force-field do FAD | **Quebrado** (elementos + LJ + torsions) — **intocado** |
| `docking_project/figures/smoke_test_fad/cofactor_FAD_TcTR.pdb` | Template canônico (86 átomos) | coordenadas/nomes OK, **coluna element embaralhada** (33/86) |
| `docking_project/figures/md_clean/minimized_dryrun_v2_prod.pdb` | Último estado válido pré-NaN | OK |
| `md_run_v2_prod.log` | Log da última execução (14/05) | NaN no warm-up |
| `auto_fix/fix_fad_xml.py` | Correção rule-based | **Não rodar** (bug de bonds duplicadas, §3) |
| `auto_fix/` (report/state/logs) | Saídas do loop | Nunca gerados |

---

*Briefing gerado automaticamente. Nenhum arquivo original do projeto foi modificado e nenhuma simulação de produção foi executada. As verificações foram parse read-only do XML/PDB e um smoke-test do FAD isolado em sandbox (OpenMM 8.5.1). Um candidato de "fix só de elementos" foi gerado, testado, reprovado e neutralizado — ver `cofactor_FAD_TcTR_v3_elemfix.xml.INVALID.txt`.*
