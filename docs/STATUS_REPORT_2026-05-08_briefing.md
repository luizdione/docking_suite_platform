# Daily Briefing — Projeto Docking & MD
**Data:** 2026-05-08 (run automático) | **Foco solicitado:** `md_module.py`

---

## 1. TL;DR — o que fazer agora

1. **`md_module.py` está OK** (3091 linhas, AST válido, gate Docking→MD funcionando, suporte a FAD/cofatores presente).
2. **A MD não está bloqueada por bug no `md_module.py`** — está bloqueada porque **NÃO HÁ pose docked nenhum dos 9 alvos**.
3. **Bug crítico real está em `docking_v2.py` rota RDKit fallback** que gerou todos os 9 PDBQTs de receptor com formatação inválida do campo charge — Vina rejeita os 9 com `PDBQT parsing error: Charge "-0.320  " is not valid.`
4. **Boa notícia:** o código fonte atual de `docking_v2.py` (linha 619) já contém o fix correto. Os PDBQTs no disco foram gerados durante uma execução anterior à correção. Bastaria deletar e re-rodar.

**Ação imediata recomendada:**
```powershell
cd C:\Users\luizd\OneDrive\Data_Science\Projetos\docking
del docking_project\proteins\*_receptor.pdbqt
del docking_project\__pycache__\docking_v2.cpython-*.pyc
.\run_docking_batch.bat
```

---

## 2. Auditoria do `md_module.py`

| Métrica | Valor |
|---|---|
| Tamanho | 150.825 bytes (3091 linhas) |
| AST | OK |
| Funções top-level | 26 |
| Classes | 1 (GPUTempReporter) |
| Imports | 91 |
| Funções duplicadas | 0 |
| Funções vazias (`pass` ou só docstring) | 0 |
| TODOs / FIXMEs / HACKs | 0 |
| `bare except:` | 0 |
| Backup `.broken_2026-05-07_backup` | quebrado (parêntese não fechado L2880, 2881 linhas, 139 KB) |

### Funções principais e tamanhos

| Linha | Função | LOC |
|---|---|---|
| 943 | `run_openmm_real` | 989 |
| 2818 | `run_md_analysis` (orchestrator) | 225 |
| 2638 | `_md_text_pt` | 173 |
| 2454 | `_md_text_en` | 181 |
| 2191 | `fig_md_rmsd` | 145 |
| 792 | `_extract_and_parametrize_cofactors` | 149 |
| 627 | `_build_lig_xml_rdkit` | 134 |
| 465 | `_build_lig_xml_smirnoff` | 143 |
| 309 | `_prepare_ligand_for_openmm` | 153 |

A `run_openmm_real` é monolítica (≈1000 linhas) — candidata óbvia a refatoração futura, mas funcional.

### O que está implementado e validado em `md_module.py`

- Gate explícito que aborta MD `real` quando não há `result_<TID>.pdbqt` (L2960-2983), retornando `status=aborted_no_docked_pose`.
- Busca de pose docked em múltiplas raízes: `docking_project/`, `_script_dir`, `E:\data\docking\` (para Dione: cobre o segundo HD).
- Detecção e parametrização de cofatores cristalográficos (FAD, FMN, NAD, NADP+, NADPH, HEM) via SMIRNOFF (rota D) ou RDKit fallback (rota E).
- `GPUTempReporter` com throttling em 88 °C (proteção da RTX 3060) e abort em 92 °C.
- Reanálise de trajetórias DCD existentes (`find_existing_trajectory` + `reanalyze_existing_trajectory`).
- Fallback automático para `simulate_md_data` (sintético) quando `OPENMM_OK=False` ou erro fatal.
- Sanitização JSON de arrays numpy/listas longas antes de serializar.

**Diagnóstico do `md_module.py`: nenhum bug pendente identificável estaticamente. O arquivo está em estado utilizável.**

---

## 3. Bug crítico real: `docking_v2.py` rota RDKit fallback

### Sintoma observado nos logs (`run_20260507_214225.log`)

Todos os 9 alvos falharam no Vina com mensagem idêntica:

```
PDBQT parsing error: Charge "-0.320  " is not valid.
 > ATOM      1  N   MET A   2       1.688  31.073 -11.885  1.00 88.37  -0.320   N
```

### Análise byte-a-byte do PDBQT gerado

Comparando o PDBQT no disco com o formato Vina esperado:

| Posição | Esperado (Vina) | No disco (errado) |
|---|---|---|
| Cols 67-70 | `"    "` (4 espaços) | `"  -0"` (charge invadiu) |
| Cols 71-76 | `"-0.320"` (charge) | `".320  "` (cauda do charge) |
| Col 77 | `" "` (separador) | `" "` |
| Cols 78-79 | `"N "` (atom type) | `" N"` (atom type deslocado) |

### Causa-raiz

O arquivo `docking_v2.py.partial_after_pdbqt_fix` (versão intermediária) tem na linha 613:

```python
line_out = pdb_line[:66] + f"{q:>8.3f}{' ':>3}{at:<2}\n"
```

O `q:>8.3f` formata `-0.320` em **8 chars right-padded** = `"  -0.320"`. Aí concatenando com 3 espaços e tipo, o charge fica nas cols 67-74 em vez de 71-76 — exatamente o padrão do bug observado.

### Estado do fix

`docking_v2.py` atual (mtime 2026-05-08 00:04, linha 619) já tem a forma correta:

```python
charge_str = f"{q:+6.3f}"  # "+0.000", "-0.320" — 6 chars exatos
line_out = pdb_line[:66] + "    " + charge_str + " " + f"{at:<2}\n"
```

Validação local com inputs reais (1BZL_clean.pdb): output produz `cols 71-76 = "-0.320"` e `cols 78-79 = "N "` — formato Vina-compliant. **O fix já está no código.**

### Por que o disco ainda tem o output errado

O batch que rodou em 2026-05-07 21:42 (`run_docking_batch.bat`) provavelmente carregou:

- **OU** o `docking_v2.py` numa versão anterior (antes do fix de 00:04)
- **OU** o `__pycache__/docking_v2.cpython-311.pyc` cacheado (que fica datado em 00:40, suspeitamente após o batch)

**Os 9 `*_receptor.pdbqt` no disco são lixo. Precisam ser regerados com a versão atual.**

---

## 4. Estado da pipeline de docking (snapshot)

### Arquivos receptor PDBQT presentes mas inválidos

| ID | PDB | Tamanho | Status |
|---|---|---|---|
| TcTR | 1BZL | 612 KB | ⚠ Charge format errado |
| LiTR | 2JK6 | – | ⚠ idem |
| TcCruzain | 3KKU | – | ⚠ idem |
| TcFeSOD | 4DVH | – | ⚠ idem |
| LmDHODH | 3MHU | – | ⚠ idem |
| hGR | 3GRS | – | ⚠ idem |
| hCatL | 2XU3 | – | ⚠ idem |
| hDHODH | 1D3H | – | ⚠ idem |
| hMnSOD | 1LUV | – | ⚠ idem |

### `result_*.pdbqt` (poses docked)

**Nenhum existente.** Todos falharam por causa do bug de charge. Como consequência:

- MD em modo `real` continuará abortando para todos os 9 alvos com `aborted_no_docked_pose` (comportamento correto e protegendo a GPU).
- Dashboard:5001 vai mostrar overlay vermelho ao tentar disparar MD real.

---

## 5. Plano de ação (passo a passo)

### Passo 1 — Limpar artifacts contaminados

```powershell
cd C:\Users\luizd\OneDrive\Data_Science\Projetos\docking
Remove-Item docking_project\proteins\*_receptor.pdbqt
Remove-Item __pycache__\docking_v2.cpython-*.pyc -ErrorAction SilentlyContinue
Remove-Item docking_project\docking\log_*.txt -ErrorAction SilentlyContinue
```

### Passo 2 — Validar versão do código antes de rodar

```powershell
python -c "import ast; src=open('docking_v2.py').read(); ast.parse(src); print('AST OK', src.count(chr(10))+1, 'linhas')"
findstr /n "charge_str" docking_v2.py
```

Esperar ver `f\"{q:+6.3f}\"` na saída — se aparecer `q:>8.3f`, o fix foi revertido.

### Passo 3 — Rodar o batch novamente

```powershell
.\run_docking_batch.bat
```

Tempo estimado: ~1h30 (9 alvos × 5-10 min).

### Passo 4 — Validar PDBQTs do receptor com 1 linha

```powershell
python -c "ln=open('docking_project/proteins/1BZL_receptor.pdbqt').readlines()[1]; print('cols71-76:', repr(ln[70:76]))"
```

Esperar ver: `cols71-76: '-0.320'` (ou outro charge válido). **Se mostrar `' .320 '` ou similar, o bug voltou.**

### Passo 5 — Confirmar que `result_*.pdbqt` foi gerado para todos

```powershell
python check_docking_status.py
```

Esperar 9 linhas `[OK] docked_flat` ou `docked_session`.

### Passo 6 — Teste curto de MD para 1 alvo

```powershell
conda activate datascience
python -c "from md_module import run_md_analysis; run_md_analysis(targets=['TcTR'], ns=10, mode='real')"
```

Validar nos logs:
- `[Cofator] Detectados: FAD`
- `[OK] Cofator FAD parametrizado e registrado`
- `[Ligante] Pose docked: ...result_TcTR.pdbqt`
- `RMSD = X.XX +/- Y.YY A | Estavel/Instavel`

---

## 6. Riscos e armadilhas

1. **Cache `.pyc`** — se o batch reusar o pyc antigo, o bug volta. Sempre limpar `__pycache__` antes de re-rodar.
2. **Encoding `\r\n` no Windows** — Python no Windows abre como text mode e converte para `\n`. Se algum dia o arquivo for editado por uma ferramenta que escreve `\r\n` literal e for lido em modo `rb`, o `[:66]` pode capturar `\r` e quebrar o output. **Recomendação preventiva:** acrescentar `.rstrip("\r\n")` em vez de `.rstrip("\n")` na linha 582 de `docking_v2.py`.
3. **OpenBabel BABEL_DATADIR** — a rota A continua falhando (logs confirmam). O ambiente `(base)` no Windows está sem a variável. Solução simples: `set BABEL_DATADIR=C:\path\to\openbabel\share\openbabel\3.x.x` ou continuar dependendo da rota B (que agora está corrigida).
4. **Espaço em disco** — alvos × snapshots DCD podem ocupar ~1-3 GB cada (50-100 ns). Confirme `E:\data\docking\` tem espaço antes de escalar.

---

## 7. Hardening sugerido (opcional, low-risk)

Patch sugerido em `docking_v2.py` (linha 582) para evitar regressão:

```python
# Antes:
pdb_line = pdb_lines[i].rstrip("\n")
# Depois (cobre arquivos com CRLF que poderiam vazar \r):
pdb_line = pdb_lines[i].rstrip("\r\n")
```

E adicionar uma validação self-test ao final de `_convert_to_pdbqt`:

```python
# Sanity check: validar primeira linha ATOM do output
with open(str(pdbqt_path)) as _f:
    for _ln in _f:
        if _ln.startswith(("ATOM", "HETATM")):
            _ch = _ln[70:76]
            try:
                float(_ch)
            except ValueError:
                return False, f"PDBQT charge field corrompido: {_ch!r}"
            break
```

Isso falha rápido em vez de deixar o Vina rejeitar 9× depois.

---

## 8. Resumo executivo

- `md_module.py` está sintaticamente íntegro e funcionalmente completo. **Não é o bloqueador.**
- O verdadeiro gargalo é o **docking falhando** em 9/9 alvos por bug de formatação no PDBQT do receptor.
- O fix já está no código atual de `docking_v2.py`, mas os PDBQTs no disco foram gerados ANTES do fix (provavelmente via cache `.pyc` antigo).
- **Plano:** limpar PDBQTs + pyc, re-rodar `run_docking_batch.bat`, validar formato, então MD destrava sozinha.

---

## 9. Notas de execução do briefing

- Execução autônoma (scheduled task `daily-briefing_docking-e-dm`).
- Sem mudanças escritas ao código (apenas leitura/validação) — convém revisar este relatório antes de aplicar o hardening da seção 7.
- Tempo investido: análise estática + diff de 3 versões de `docking_v2.py` + simulação byte-a-byte do output.

**Próxima execução:** rodar de novo após limpar `*_receptor.pdbqt` e re-executar `run_docking_batch.bat`. Se o problema voltar, consultar a seção 6 (riscos).
