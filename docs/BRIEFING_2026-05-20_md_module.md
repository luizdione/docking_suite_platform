# Daily Briefing — Docking / MD de TcTR + 2-SCN-NQ

**Data:** 2026-05-20 (run automático do `daily-briefing_docking-e-dm`)
**Foco:** módulo de MD (`md_clean.py`, sucessor do `md_module.py`)
**Última atividade no projeto:** 2026-05-14 — **6 dias parado**

---

## TL;DR

O pipeline de MD não avançou desde 14/05. O `md_module.py` original (3.402 linhas, monolítico) já foi aposentado para `obsoletos/` e substituído pelo `md_clean.py` (858 linhas), que roda **end-to-end até a minimização** sem problema. O bloqueio **não está no código Python** — está no force-field do FAD (`cofactor_FAD_TcTR.xml`), que explode em NaN no warm-up.

Rodei hoje uma verificação independente do XML e ela **corrige um ponto central do diagnóstico anterior**: o `HarmonicBondForce` **NÃO está vazio** (tem 91 bonds válidas, incluindo a N3A–C4A que o forensic acusava de faltar). O problema real e confirmado é (1) **33/86 atom types com elemento/massa/LJ trocados** e (2) **`PeriodicTorsionForce` totalmente vazio (0 torsions)**.

Consequência prática: o `auto_fix/fix_fad_xml.py` — que foi preparado em 14/05 mas **nunca chegou a rodar** — corrige os elementos mas **não adiciona torsions**, então provavelmente ainda daria NaN. O caminho mais seguro continua sendo **regenerar o XML via OpenFF Toolkit** (gera elementos, bonds, angles e torsions consistentes de uma vez).

---

## 1. Estado do pipeline

O `md_clean.py --fresh` reconstrói o sistema completo com sucesso (log `md_run_v2_prod.log`, 14/05):

- PDBFixer → addHydrogens(proteína) → 2× FAD canônicos via Kabsch → LIG → solvatação TIP3P + 0,15 M NaCl
- Sistema final: **183.179 átomos** (TcTR dímero + 2 FAD + 1 LIG + 55.906 HOH + 311 Na/Cl)
- Minimização converge: PotentialEnergy = **−3,00 × 10⁶ kJ/mol**, 0 átomos com NaN pós-minimização
- PDB minimizado salvo (último estado válido): `docking_project/figures/md_clean/minimized_dryrun_v2_prod.pdb`

O **warm-up (dt=0,25 fs @ 5 K) explode em NaN imediatamente**. O forensic pré-warm-up aponta `BOND ANOMALA 2.06 Å: N3A–C4A` no FAD e força máxima de ~7,8×10² kJ/mol/nm concentrada em átomos do anel adenina do FAD.

---

## 2. Diagnóstico do FAD XML — verificado e corrigido hoje

Rodei um parse independente de `cofactor_FAD_TcTR.xml` (mtime 12/05, intocado). Resultado comparado ao que dizem `FAD_XML_SCRAMBLE_DIAGNOSIS.md` e `NEXT_SESSION.md` (ambos de 14/05):

| Item verificado | Diagnóstico anterior (14/05) | Verificação de hoje (20/05) | Veredito |
|---|---|---|---|
| Element/massa/LJ por atom type | 33/86 errados | **33/86 errados** (PA→C, N3A→O, C4A→P, N9A→C…) | ✅ confirmado — bug real |
| `HarmonicBondForce` | "vazio, 0 bonds" | **91 bonds, k de 1,6×10⁵–6,8×10⁵, nenhuma com k=0** | ❌ diagnóstico anterior incorreto |
| Bond N3A–C4A | "não existe na HBF, anel não fecha" | **existe** (length 0,163 nm, k=2,16×10⁵) | ❌ diagnóstico anterior incorreto |
| `HarmonicAngleForce` | não avaliado | 157 angles com k válido | ✅ OK |
| `PeriodicTorsionForce` | não enfatizado | **0 torsions (vazio)** | ⚠️ defeito real, não corrigido por nada ainda |
| `NonbondedForce` (cargas) | não avaliado | 86 átomos, soma de cargas = 0,000 (neutro) | ✅ OK |

### Causa-raiz corrigida do NaN

Não é "bond sem mola". As bonds e angles **têm** parâmetros. O que está scrambled é o mapeamento **atom type → (elemento, massa, sigma/epsilon)**:

- `C4A` recebe o type FAD21 com **element=P, massa=30,97** (deveria ser C, 12,01)
- `N3A` recebe FAD20 com **element=O, massa=15,99** (deveria ser N, 14,01)
- `PA` (o fósforo real) recebe FAD0 com **element=C, massa=12,01** (deveria ser P, 30,97)

Ou seja: massas erradas + raios de Lennard-Jones errados em 38% dos átomos pesados. Na minimização (dominada por LJ/Coulomb, sem inércia) o sistema acomoda numa geometria torta — daí a "bond" N3A–C4A a 2,06 Å lutando contra o comprimento harmônico de 1,63 Å. No warm-up, massas erradas geram velocidades erradas (v ∝ 1/√m) e a geometria tensionada libera energia → NaN. A **ausência total de torsions** agrava: nada mantém a planaridade dos anéis adenina/isoaloxazina.

---

## 3. Pendências e bloqueios

- **`auto_fix/` nunca rodou.** O pacote rule-based offline (`auto_md_loop.py`, `fix_fad_xml.py`, `validate_system.py`) foi criado em 14/05 21:44, mas não há `auto_md_report.md`, `auto_md_state.json` nem `auto_md_logs/`. Não foi executado.
- **`fix_fad_xml.py` provavelmente é insuficiente sozinho.** Ele corrige elementos e "adiciona bonds faltando" — mas bonds **não** estão faltando, e ele **não** adiciona torsions (que estão faltando de verdade). O próprio README do `auto_fix` antecipa isso: se após o fix de elementos ainda explodir, o gargalo é `PeriodicTorsionForce` vazio → regenerar via OpenFF.
- **CUDA do OpenMM quebrado.** Os 4 plugins CUDA falham ao carregar (`DLL load failed 126 = MOD_NOT_FOUND`); o pipeline caiu para OpenCL com mixed precision. Independente do bug do FAD, performance está comprometida.

---

## 4. Próximos passos recomendados (prioritizados)

1. **Regenerar `cofactor_FAD_TcTR.xml` via OpenFF Toolkit (Caminho A).** É o único caminho que entrega elementos/massas/LJ corretos **e** torsions **e** bonds/angles consistentes de uma vez. Atenção à limitação documentada: OpenFF emite um *OpenMM System* serializado, não um *ForceField* template — então `md_clean.py` (linhas 258-259, que aplicam `app.ForceField(...)` em runtime) precisaria de um pequeno ajuste para carregar o System serializado, ou usar `openmmforcefields` `SMIRNOFFTemplateGenerator` (que o `md_clean.py` já cita nos créditos) para gerar um template reusável. Recomendo o `SMIRNOFFTemplateGenerator`, que evita mexer no fluxo do pipeline.
2. **Em paralelo, rodar `auto_fix/auto_md_loop.py` uma vez** só para registrar empiricamente o comportamento pós-fix-de-elementos e confirmar a previsão de que torsions vazias ainda causam NaN. É barato e gera log estruturado — mas **não** trate como solução definitiva.
3. **Consertar o CUDA** (`conda install -c conda-forge openmm-cuda`; conferir `cudart64_*.dll` em `...\envs\datascience\Library\bin`). Necessário para escalar de 0,1 ns para 5/50 ns em tempo viável.
4. Após XML corrigido: re-rodar `python md_clean.py --fresh --ns 1 --tag prod_v3 --platform CUDA --allow-fallback` e confirmar que o forensic **não** acusa `BOND ANOMALA`. Se limpo, escalar.

> ⚠️ Como já anotado em `NEXT_SESSION.md`: **não** reiniciar ciclos automáticos de MD sem o fix do XML — vai dar NaN sempre.

---

## 5. Arquivos-chave

| Caminho | Papel | Estado |
|---|---|---|
| `md_clean.py` (858 linhas) | Pipeline MD ativo (substitui `md_module.py`) | OK até minimização |
| `obsoletos/md_module.py` (3.402 linhas) | Módulo MD antigo, monolítico | Aposentado |
| `docking_project/figures/smoke_test_fad/cofactor_FAD_TcTR.xml` | Force-field do FAD | **Quebrado** (elementos + torsions) |
| `docking_project/figures/md_clean/minimized_dryrun_v2_prod.pdb` | Último estado válido pré-NaN | OK |
| `md_run_v2_prod.log` | Log da última execução (14/05) | NaN no warm-up |
| `auto_fix/` | Pacote de correção offline | Criado, **nunca executado** |
| `_fad_diagnose.py` | Script de diagnóstico do XML | — |

---

*Briefing gerado automaticamente. Nenhum arquivo do projeto foi modificado e nenhuma simulação foi executada — esta é uma tarefa de leitura/relatório. As verificações do XML foram feitas por parse read-only em sandbox isolado.*
