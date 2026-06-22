# Relatório de Status — NQ-Docking Suite v2
**Data:** 2026-03-25 (atualização #2) | **Execução automática (scheduled task)**

---

## 1. Correções Realizadas Nesta Sessão

### Bug #1: SVG com nome hardcoded (admet_module.py)
**Problema:** Linha 836 gerava SVG sempre como `2SCNNQ_structure.svg`, independentemente do composto ativo.
**Correção:** Substituído por `{safe_abbr}_structure.svg` usando `COMPOUND_ABBR` sanitizado.
**Arquivo:** `admet_module.py` (linha 836)

### Bug #2: SMILES incorreto do Nifurtimox
**Problema:** SMILES `O=S1(=O)CCN(N=Cc2ccc([N+](=O)[O-])o2)CC1` gera C9H11N3O5S (MW=273.27), mas Nifurtimox é C10H13N3O5S (MW=287.29). Faltava o grupo metil no anel tiomorfolina.
**Correção:** SMILES atualizado para `CC1CS(=O)(=O)CCN1/N=C/c1ccc([N+](=O)[O-])o1` (MW=287.30 ✓, fórmula C10H13N3O5S ✓).
**Arquivos:** `populate_ligand_db.py` (linha 61), `docking_project/ligand/Nifurtimox/metadata.json`
**Nota:** Os arquivos .mol2 e .pdbqt do Nifurtimox precisam ser regenerados na máquina local com o SMILES corrigido.

### Bug #3: Globais não restauradas após run_admet()
**Problema:** `run_admet()` modificava `COMPOUND_NAME`/`COMPOUND_ABBR` globais mas nunca as restaurava, causando contaminação em chamadas subsequentes.
**Correção:** Adicionado bloco `try/finally` que restaura os valores originais ao fim da execução.
**Arquivo:** `admet_module.py`

### Melhoria #4: Parâmetros compound_name/compound_abbr em toxicity_module e cyclodex_module
**Problema:** `run_toxicity()` e `run_cyclodextrin()` não aceitavam `compound_name`/`compound_abbr` como parâmetros. O dashboard manipulava globais diretamente sem restaurá-las.
**Correção:**
- `toxicity_module.py`: Adicionados parâmetros `compound_name`, `compound_abbr`, `exp_data` a `run_toxicity()` com padrão try/finally.
- `cyclodex_module.py`: Adicionados parâmetros `compound_name`, `compound_abbr`, `ligand_smiles` a `run_cyclodextrin()` com padrão try/finally.
- `dashboard_v2.py`: Rotas `/run_toxicity` e `/run_cyclodextrin` agora passam `compound_name` via parâmetros em vez de mutação direta de globais.

---

## 2. Integridade dos Arquivos

Todos os 12 scripts Python passaram na verificação de sintaxe (AST parse) sem erros após as correções:
`admet_module.py`, `analysis.py`, `compound_db.py`, `cyclodex_module.py`, `dashboard_v2.py`, `diagnostic.py`, `drug_design_module.py`, `md_module.py`, `populate_ligand_db.py`, `setup_openmm_venv.py`, `target_selector.py`, `toxicity_module.py`

---

## 3. Status por Composto (sem alteração)

### 2-SCN-NQ (composto principal) — COMPLETO
- Docking: 10 alvos — resultados e logs presentes
- ADMET: perfil completo com figuras
- MD: 3 sessões sintéticas + 1 sessão OpenMM real
- Figuras de docking: PNG + SVG + TIFF
- Textos para manuscrito: docking, MD, ADMET (EN)

### Lapachol — PARCIAL (ADMET apenas)
### Nifurtimox — PARCIAL (ADMET apenas, SMILES corrigido nesta sessão)
### alpha-Lapachone — PARCIAL (ADMET apenas)
### NQ1-naphthoimidazole — INTERROMPIDO (nenhum resultado)
### 14 compostos restantes — APENAS LIGANTES PREPARADOS

---

## 4. Tarefas Pendentes (prioridade atualizada)

1. ~~Corrigir bug SVG~~ ✓ FEITO
2. ~~Corrigir SMILES do Nifurtimox~~ ✓ FEITO
3. **Regenerar ligantes .mol2/.pdbqt do Nifurtimox** com SMILES corrigido (máquina local)
4. **Re-executar ADMET do Nifurtimox** para atualizar figuras e perfil com MW correto
5. **Executar docking** para Lapachol, Nifurtimox, alpha-Lapachone (sessões criadas mas vazias)
6. **Completar pipeline NQ1-naphthoimidazole** — totalmente interrompido
7. **Executar pipeline completo** para os 14 compostos restantes
8. **MD, toxicidade, ciclodextrina** para todos os compostos além do 2-SCN-NQ

---

## 5. Nota sobre Execução

Correções de código e verificações foram realizadas com sucesso neste ambiente sandbox. Tarefas de docking molecular (AutoDock Vina) e MD (OpenMM) requerem execução na máquina local Windows com GPU.
