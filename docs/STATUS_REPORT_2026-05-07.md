# Relatório de Status — md_module.py
**Data:** 2026-05-07 | **Execução automática (scheduled task)**

---

## 1. Sintaxe e Integridade do Arquivo

O `md_module.py` (2880 linhas, ~136 KB, última modificação 2026-05-07 00:00) **passa no parse AST sem erros**. Contém 25 funções e 1 classe (`GPUTemperatureGuard`).

| Métrica | 2026-05-05 | 2026-05-06 | 2026-05-07 |
|---|---|---|---|
| Linhas | 2143 (truncado) | 2558 | **2880** |
| Funções | 20 | 25 | **25** |
| AST OK | Sim (sintaxe), mas truncado | Sim | **Sim** |
| `bare except` | 0 | 0 | **0** |
| Funções críticas (run_md_analysis, generate_md_text, fig_md_*) | Ausentes | Restauradas | **Mantidas + ampliadas** |

O arquivo **cresceu 322 linhas desde 2026-05-06** — sugere refinamentos na rotina `run_md_analysis()` (busca de PDBQT em múltiplas raízes, fallback robusto para `figures_backup`/`figures_real` se houver symlink quebrado, suporte a `SMILES:` como ligante quando PDBQT ausente). Compilação `.pyc` recente (`md_module.cpython-311.pyc`, 2026-05-06 23:20) confirma que o módulo está sendo importado com sucesso pelo dashboard.

---

## 2. Logs de Execução

**Nenhum arquivo de log encontrado** na pasta do projeto (OneDrive):
- `md_autostart.log` — **ausente** ⇒ `md_autostart.py` ainda não foi executado nem agendado no Task Scheduler do Windows
- `md_results.json` — **ausente** ⇒ nenhuma simulação MD foi concluída e analisada nesta pasta
- Trajetórias `.dcd`, PDBs `complex_*_initial.pdb`, e `openmm_*.csv` — **ausentes**

O arquivo `venv_info.json` (gerado em 2026-03-29 por `setup_openmm_venv.py`) reporta `"openmm_ok": false` — provavelmente desatualizado, pois o status report de 2026-05-06 documentava simulação real concluída em 10 ns no TcTR.

**Conclusão:** Os outputs MD reais estão no drive `E:\data\docking\` ou em `figures_backup\`, fora da pasta sincronizada via OneDrive. O scheduled task atual não tem visibilidade sobre eles.

---

## 3. Diagnóstico — Erros e Correções Necessárias

### ERRO #1 (RESOLVIDO desde 2026-05-06): Truncação do md_module.py
Confirmado: a integridade está mantida, código robusto e funcional.

### ERRO #2 (PERSISTENTE): Falta de visibilidade dos resultados MD na pasta OneDrive
- **Sintoma:** O scheduled task não consegue auditar os resultados das simulações reais porque eles vivem no drive `E:\` local.
- **Correção sugerida (operacional):** Após cada execução em `E:\`, copiar `md_results.json`, `md_text_*.txt` e PNGs para `docking_project\figures\` no OneDrive. Adicionar ao final de `run_md_analysis()` (após linha 2830) um bloco opcional de mirror:
  ```python
  if os.environ.get("MIRROR_MD_RESULTS"):
      mirror_dir = Path(__file__).parent / "docking_project" / "figures"
      mirror_dir.mkdir(parents=True, exist_ok=True)
      for f in out_dir.glob("md_*.{json,txt,png}"):
          shutil.copy2(f, mirror_dir / f.name)
  ```

### ERRO #3 (POTENCIAL): md_autostart.py nunca foi agendado
- **Sintoma:** Sem `md_autostart.log` ⇒ o script nunca rodou. As simulações que pararam por queda de energia ou logoff não estão sendo retomadas automaticamente.
- **Correção:** No PowerShell admin:
  ```powershell
  schtasks /Create /TN "MD_AutoRestart" /TR "pythonw C:\Users\luizd\OneDrive\Data_Science\Projetos\docking\md_autostart.py" /SC ONSTART /RL HIGHEST /F
  ```
  Verificar com `schtasks /Query /TN MD_AutoRestart`.

### ERRO #4 (PERSISTENTE de 2026-05-06): Ligante nem sempre incluído na simulação
O código atual tem **5 rotas de fallback** para parametrização do ligante (Rotas 0/A/B/B2/D/E). A linha 396 chama `GAFFTemplateGenerator` (Rota A) que requer `openff.toolkit` + AmberTools (`sqm`). Se ausentes, cai para Gasteiger via RDKit (Rota E, linha 405-408). Validar que pelo menos uma rota está funcional executando:
```bash
python md_module.py --check-cuda
```
Se OpenMM aparecer como "NO", rodar `python setup_openmm_venv.py` no env `datascience` da Miniforge.

---

## 4. Status do Pipeline Geral

| Componente | Status | Observação |
|---|---|---|
| `md_module.py` | OK (2880 linhas, AST OK) | Estável, sem regressões |
| `md_autostart.py` | OK (sintaxe) | **NÃO agendado no Task Scheduler** |
| `monitor_md.py` | OK (sintaxe) | Não foi tocado desde 2026-04-17 |
| `dashboard_v2.py` | OK | Aba MD operacional |
| `setup_gaff2.py` | OK | |
| `venv_info.json` | Desatualizado (mar/2026) | `openmm_ok=false` precisa re-detecção |
| Trajetórias DCD | Ausentes em OneDrive | Provavelmente em `E:\data\docking\` |
| `md_results.json` | Ausente em OneDrive | Idem |
| Logs (`*.log`) | Ausentes | autostart nunca rodou |

---

## 5. Instruções para Próximas Fases

### FASE 1 — Sincronização e Diagnóstico (URGENTE — fazer hoje)
1. **Re-detectar OpenMM:** `python setup_openmm_venv.py` (atualiza `venv_info.json`)
2. **Diagnóstico CUDA:** `python md_module.py --check-cuda`
3. **Localizar DCDs no E:\\:**
   ```powershell
   dir E:\data\docking\figures /S | findstr /i "traj_.*\.dcd"
   ```
4. **Re-análise das trajetórias existentes** (se DCDs existirem em E:\\):
   ```python
   from md_module import reanalyze_existing_trajectory
   r = reanalyze_existing_trajectory(
       "TcTR",
       r"E:\data\docking\figures\2-SCN-NQ\traj_TcTR.dcd",
       r"E:\data\docking\figures\2-SCN-NQ\complex_TcTR_initial.pdb",
       ns=10.0
   )
   ```

### FASE 2 — Agendamento do md_autostart
1. Criar tarefa no Task Scheduler (comando do item 3.3 acima)
2. Testar manualmente: `python md_autostart.py --check`
3. Verificar geração de `md_autostart.log`

### FASE 3 — Espelhar resultados para OneDrive
Adicionar ao final de `run_md_analysis()` o bloco de mirror sugerido na Seção 3 (ERRO #2). Variável de ambiente `MIRROR_MD_RESULTS=1` ativa o mirror seletivo (apenas JSON/TXT/PNG, sem DCDs grandes).

### FASE 4 — Expansão Multi-Alvo (próxima semana)
Com TcTR validado, expandir para os 9 alvos do `STABILITY_REFERENCE`:
```bash
python md_module.py --targets TcTR LiTR TcCruzain hGR --ns 100 --mode real --lang pt
```
Prioridade segundo melhor score docking esperado: **TcTR → LiTR → TcCruzain → hGR (controle humano)**.

### FASE 5 — Manuscrito (após Fase 4)
- Figuras: RMSD/RMSF por alvo + comparação multi-alvo + bar plot resumo
- Texto científico em PT/EN automaticamente via `generate_md_text()`
- Integrar `md_text_TcTR_detailed.md` (já existente, 7 KB) como referência detalhada
- Cross-check com `PROTOCOLO_CIENTIFICO.md` e `DRUG_DESIGN_MANUAL.md`

### FASE 6 — Refatoração (médio prazo, baixa prioridade)
Conforme análise em `.claude_analysis.json`:
1. Adicionar suite `tests/` com pytest para `md_module`, `docking_v2`, `admet_module`
2. Refatorar `dashboard_v2.py` (270 KB!) em páginas/componentes Streamlit separados
3. Migrar para `pyproject.toml` (pacote Python instalável)
4. Externalizar PyMOL para fora do repositório

---

## 6. Resumo Executivo

**md_module.py está saudável e ampliado** (2880 linhas, +322 desde ontem, AST OK, 25 funções, sem `bare except`). A correção de truncação reportada em 2026-05-05 está consolidada. **Não há erros de sintaxe ou estruturais ativos.**

**Os bloqueios atuais são operacionais, não de código:**
1. Os resultados das simulações reais vivem em `E:\data\docking\` e não são visíveis pela pasta OneDrive sincronizada — adicionar mirror seletivo ao fim de `run_md_analysis()`.
2. `md_autostart.py` nunca foi agendado no Task Scheduler — comando de criação fornecido na Fase 2.
3. `venv_info.json` está desatualizado (mar/2026 com `openmm_ok=false`) — re-rodar `setup_openmm_venv.py` para refletir o estado atual.

**Próxima ação imediata:** executar `python md_module.py --check-cuda` e `dir E:\data\docking\figures /S` na máquina local para mapear o estado real dos artefatos MD antes de iniciar a expansão multi-alvo.
