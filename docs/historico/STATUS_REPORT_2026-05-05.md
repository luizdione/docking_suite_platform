# Relatório de Status — md_module.py
**Data:** 2026-05-05 | **Execução automática (scheduled task)**

---

## 1. Sintaxe do Arquivo

O `md_module.py` (2143 linhas, ~103 KB, última modificação 2026-05-04) passa no parse AST sem erros de sintaxe. Contém 20 funções. Nenhum `bare except` encontrado.

---

## 2. ERROS CRÍTICOS IDENTIFICADOS

### ERRO #1: Arquivo truncado — funções essenciais ausentes

O arquivo `md_module.py` está **truncado na linha 2144**. A última linha é `            wi` (texto cortado no meio da palavra "window"). Isso significa que as seguintes funções estão **completamente ausentes**:

- **`fig_md_comparison()`** — incompleta (cortada no meio do loop)
- **`fig_md_summary_bar()`** — ausente
- **`generate_md_text()`** — ausente
- **`run_md_analysis()`** — ausente (entry point principal)
- **`if __name__ == "__main__":`** — ausente (CLI)

**Impacto:** O dashboard_v2.py chama `md_module.run_md_analysis()` e `md_module.generate_md_text()` em pelo menos 5 pontos (linhas 3869, 3878, 4144, 4368, 4422). Sem estas funções, toda a aba de MD do dashboard está **quebrada** — lançará `AttributeError` ao ser acessada.

**Causa provável:** Edição anterior do arquivo foi interrompida (provavelmente por limite de créditos/contexto), resultando em arquivo incompleto.

**Correção necessária:** Restaurar as funções ausentes a partir do backup `md_module(V1).py` (que contém versões funcionais dessas funções nas linhas 1130–1491), adaptando para as novas funções adicionadas na V2 (como `reanalyze_existing_trajectory`, `find_existing_trajectory`, `_build_lig_xml_smirnoff`, etc.)

### ERRO #2: Análise de trajetória falhou — RMSD/RMSF todos zeros

Os resultados da última simulação OpenMM real (sessão 2026-04-09_19-54-17) mostram:

| Alvo | Status | RMSD mean | RMSF mean | Ligante | Erro |
|------|--------|-----------|-----------|---------|------|
| TcTR | complete | 0.0 | 0.0 | Não incluído | Could not open file (path relativo) |
| hGR  | complete | 0.0 | 0.0 | Não incluído | Could not open file (path relativo) |

**Causa:** O MDTraj não conseguiu abrir os arquivos `.dcd` porque o path armazenado usa caminhos Windows (`docking_project\figures\2-SCN-NQ\2026-04-09_19-54-17\traj_TcTR.dcd`) que dependem do CWD estar correto no momento da análise. O código na linha 1635–1637 já tenta usar `.resolve()`, mas o path relativo original falha quando o script roda de outro diretório.

**Nota:** As trajetórias DCD existem e são válidas:
- `traj_TcTR.dcd` — 4.67 GB (na pasta figures/)
- `traj_hGR.dcd` — 2.12 GB (na pasta figures_backup/)

Os logs CSV confirmam simulações completas (5.1M steps, 10 ns, ~24.9 ns/day para TcTR e ~55.3 ns/day para hGR, plataforma OpenCL).

**Correção necessária:** Re-executar a análise MDTraj usando `reanalyze_existing_trajectory()` com paths absolutos.

### ERRO #3: Ligante não incluído na simulação

Ambas as simulações rodaram como **proteína-only** (`ligand_included: False`, `gaff2_used: False`). O ligante 2-SCN-NQ não foi parametrizado com GAFF2, provavelmente por falta de AmberTools/openmmforcefields configurados corretamente.

---

## 3. Funcionalidades Operacionais (OK)

- `get_md_status()` — diagnóstico do ambiente OpenMM ✓
- `_probe_openmm()` — detecção de GPU com fallback OpenCL ✓
- `_pdbqt_to_clean_pdb_text()` — conversão PDBQT→PDB ✓
- `_prepare_ligand_for_openmm()` — pipeline PDBQT→SDF→GAFF2 com 3 fallbacks ✓
- `_build_lig_xml_smirnoff()` — parametrização SMIRNOFF ✓
- `_build_lig_xml_rdkit()` — parametrização RDKit fallback ✓
- `run_openmm_real()` — simulação OpenMM completa (764 linhas) ✓
- `GPUTemperatureGuard` — proteção térmica da GPU ✓
- `find_existing_trajectory()` — busca de trajetórias em backup ✓
- `reanalyze_existing_trajectory()` — re-análise MDTraj ✓
- `calculate_rmsd_rmsf()` — análise MDAnalysis ✓
- `simulate_md_data()` — dados sintéticos baseados em literatura ✓
- `fig_md_rmsd()` — figura RMSD/RMSF por alvo ✓

---

## 4. Instruções para Correção (Próxima Sessão)

### Prioridade ALTA — Restaurar funções truncadas

1. Abrir `md_module.py` e localizar a última linha (2144: `            wi`)
2. Substituir a partir da linha truncada com o código completo de `fig_md_comparison()` (da V1, linha 1162 em diante)
3. Adicionar `fig_md_summary_bar()` (V1, linhas 1183–1227)
4. Adicionar `generate_md_text()` e helpers `_md_text_en()`, `_md_text_pt()` (V1, linhas 1230–1328)
5. Adicionar `run_md_analysis()` (V1, linhas 1330–1445) — adaptar para usar `find_existing_trajectory()` e `reanalyze_existing_trajectory()` que são novas na V2
6. Adicionar bloco `if __name__ == "__main__":` (V1, linhas 1447–1491)
7. Verificar sintaxe com `python -c "import ast; ast.parse(open('md_module.py').read())"`

### Prioridade MÉDIA — Re-analisar trajetórias existentes

Após restaurar as funções, executar na máquina local:
```python
from md_module import reanalyze_existing_trajectory
# TcTR
r1 = reanalyze_existing_trajectory(
    "TcTR",
    r"C:\Users\luizd\OneDrive\Data_Science\Projetos\docking\docking_project\figures\2-SCN-NQ\traj_TcTR.dcd",
    r"C:\Users\luizd\OneDrive\Data_Science\Projetos\docking\docking_project\figures\2-SCN-NQ\complex_TcTR_initial.pdb",
    ns=10.0
)
# hGR
r2 = reanalyze_existing_trajectory(
    "hGR",
    r"C:\Users\luizd\OneDrive\Data_Science\Projetos\docking\docking_project\figures_backup\2-SCN-NQ\2026-04-09_19-54-17\traj_hGR.dcd",
    r"C:\Users\luizd\OneDrive\Data_Science\Projetos\docking\docking_project\figures_backup\2-SCN-NQ\2026-04-09_19-54-17\complex_hGR_initial.pdb",
    ns=10.0
)
```

### Prioridade BAIXA — Incluir ligante na próxima simulação

Para rodar MD com o ligante parametrizado:
1. Instalar AmberTools: `conda install -c conda-forge ambertools`
2. Ou usar `--mode real` que já tem os fallbacks Gasteiger/SMIRNOFF
3. Re-executar: `python md_module.py --target TcTR hGR --ns 10 --mode real`

---

## 5. Status Geral do Projeto

| Módulo | Status |
|--------|--------|
| md_module.py | PARCIAL — truncado, precisa restaurar ~350 linhas |
| dashboard_v2.py | OK (sintaxe), mas aba MD quebrada por falta de funções em md_module |
| admet_module.py | OK |
| docking_v2.py | OK |
| toxicity_module.py | OK |
| cyclodex_module.py | OK |
| drug_design_module.py | OK |
| monitor_md.py | OK |
| notion_module.py | OK |
| setup_gaff2.py | OK (novo) |

### Compostos

- **2-SCN-NQ:** Docking completo, MD real executada (mas análise falhou), ADMET completo
- **Lapachol, Nifurtimox, alpha-Lapachone:** Apenas ADMET
- **14 compostos restantes:** Apenas ligantes preparados

---

## 6. Resumo Executivo

O problema mais urgente é que `md_module.py` está **truncado** — faltam ~350 linhas incluindo as funções `run_md_analysis()`, `generate_md_text()`, `fig_md_summary_bar()` e o final de `fig_md_comparison()`. Isso quebra completamente a aba de MD no dashboard. A correção é restaurar essas funções a partir do backup V1, adaptando para as melhorias da V2. Além disso, as trajetórias DCD existentes (TcTR: 4.67 GB, hGR: 2.12 GB) nunca foram analisadas corretamente porque o MDTraj falhou ao abrir os arquivos com paths relativos — basta re-executar `reanalyze_existing_trajectory()` com paths absolutos.
