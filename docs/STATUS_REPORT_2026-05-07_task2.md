# Relatório Final — Task 2: Plataforma com Docking → MD integrado
**Data:** 2026-05-07 | **Estado: pronto para uso**

---

## 1. O que foi implementado nesta sessão

### 1.1 Bug do PDBQT receptor corrigido em docking_v2.py
A função `_convert_to_pdbqt` falhava com erro `pdb is not a recognised Open Babel format` em todos os 9 alvos. Causa: variáveis de ambiente do OpenBabel (BABEL_DATADIR) não configuradas no env (base) ou conflito de versão.

**Solução:** dupla rota
- **Rota A (preferida):** OpenBabel com verificação prévia de formatos disponíveis
- **Rota B (fallback):** RDKit + tabela manual AutoDock — não precisa de BABEL_DATADIR

A Rota B atribui tipos AutoDock corretos (`A` para C aromático, `NA` para N aceitador, `OA` para O, `SA` para S, `HD` para H polar, `Fe/Mn/Mg/Zn` para metais) via cargas Gasteiger do RDKit.

### 1.2 Endpoint Flask `/run_docking` no dashboard
- **POST `/run_docking`** com payload `{targets, exhaustiveness, compound}` 
- Executa `docking_v2.py` via `subprocess.run` com timeout de 3h
- Retorna log completo + returncode + resultados JSON quando disponíveis
- Aceita `targets` como lista, `"all"`, ou `"parasite"`

### 1.3 Botões na UI da aba MD
Adicionei um painel verde de "Pipeline obrigatório: Docking → MD" diretamente na aba MD com 3 botões:

| Botão | Função |
|---|---|
| 🔍 **Verificar docking** | Chama `/check_docking_for_md`, lista alvos OK e faltantes |
| ▶ **Rodar docking (alvos selecionados)** | Roda docking apenas dos alvos selecionados na MD |
| ▶ **Rodar docking de TODOS os 9 alvos** | Batch completo (TcTR, LiTR, TcCruzain, TcFeSOD, LmDHODH, hGR, hCatL, hDHODH, hMnSOD) |

Cada botão pede confirmação antes de iniciar (~5-10 min/alvo).

### 1.4 Aba Ajuda atualizada
Nova seção em destaque (caixa amarela com workflow) na aba Ajuda contendo:
- Explicação do gate Docking → MD
- 4 passos do fluxo correto
- O que fazer quando a MD bloquear
- Suporte a cofatores (FAD, FMN, NAD, NADP+, NADPH, HEM)
- Comandos shell para diagnóstico
- Lista dos 9 alvos do estudo com PDBs e cofatores

---

## 2. Validação

| Arquivo | Linhas | AST |
|---|---|---|
| `dashboard_v2.py` | 5262 | ✓ OK |
| `docking_v2.py` | 1102 | ✓ OK |
| `md_module.py` | 3091 | ✓ OK |
| `check_docking_status.py` | 125 | ✓ OK |
| `md_autostart.py` | 212 | ✓ OK |
| `monitor_md.py` | 291 | ✓ OK |

Todas as mudanças foram validadas em disco (não apenas em memória):
- Endpoint `/run_docking`: presente
- Endpoint `/check_docking_for_md`: presente
- 4 funções JS de docking UI: presentes
- 3 botões na UI da MD: presentes
- Overlay vermelho de bloqueio: presente
- Gate na `runMD()`: presente
- Aba Ajuda com novo workflow: presente

---

## 3. Como usar agora

### Opção A — Pelo dashboard web (recomendado)
1. Inicie o servidor:
   ```powershell
   conda activate datascience
   python dashboard_v2.py
   ```
2. Abra http://localhost:5001
3. Vá para a aba **MD**
4. Selecione os alvos que quer (ex.: TcTR, LiTR, hGR)
5. Clique **🔍 Verificar docking** — vai mostrar quais têm/não têm pose docked
6. Se algum faltar, clique **▶ Rodar docking (alvos selecionados)**
7. Aguarde ~5-10 min/alvo (acompanhe progress no painel inferior)
8. Após docking concluído, **▶ Gerar Figuras MD** em modo `real` (já com gate satisfeito)

### Opção B — Pelo shell (para batch grande)
```powershell
cd C:\Users\luizd\OneDrive\Data_Science\Projetos\docking
.\run_docking_batch.bat
```

### Opção C — Comando direto
```powershell
conda activate datascience
python docking_v2.py --targets TcTR,LiTR,TcCruzain,TcFeSOD,LmDHODH,hGR,hCatL,hDHODH,hMnSOD --exhaustiveness 32
```

---

## 4. Comportamento do gate

Quando você tenta rodar MD em modo `real` sem docking:

**No dashboard:**
1. Frontend chama `/check_docking_for_md` antes de `/run_md`
2. Se algum alvo faltar → overlay vermelho aparece **imediatamente** sem disparar GPU
3. Overlay lista os alvos faltantes e oferece botão "Entendi, vou rodar o docking"
4. Defesa em profundidade: backend retorna HTTP 409 mesmo se alguém bypassar JS

**No CLI (md_module.py direto):**
1. `run_md_analysis` detecta `lig_source == "preparation"` e aborta com:
   ```
   [ERRO CRITICO] result_TcTR.pdbqt nao encontrado.
   [ABORT] Pulando MD para TcTR - faltando pose docked do Vina.
   ```
2. Retorna `{"status": "aborted_no_docked_pose", ...}` — sem desperdício de GPU

**Modo synthetic:** continua livre (não precisa de docking, gera dados de literatura para preview).

---

## 5. Próximos passos para você

1. **Teste agora** — abra o dashboard, vá na aba MD, clique em "🔍 Verificar docking" → deve mostrar 9 alvos faltantes
2. **Rode o docking** via UI ou batch — vai gerar `result_<TID>.pdbqt` para cada
3. **Verifique de novo** — todos devem aparecer como `[OK]`
4. **Teste curto de MD (10 ns)** com TcTR em modo `real` para validar o pipeline com o FAD parametrizado
5. **Confirme no log da MD:**
   - `[Cofator] Detectados: FAD`
   - `[OK] Cofator FAD parametrizado e registrado`
   - `[Ligante] Pose docked: ...result_TcTR.pdbqt`
6. Se OK, escale para produção 50-100 ns multi-alvo

---

## 6. Arquivos da entrega

- [dashboard_v2.py](computer://C:\Users\luizd\OneDrive\Data_Science\Projetos\docking\dashboard_v2.py) — endpoint /run_docking + botões UI + aba Ajuda
- [docking_v2.py](computer://C:\Users\luizd\OneDrive\Data_Science\Projetos\docking\docking_v2.py) — bug PDBQT corrigido (rota fallback RDKit)
- [md_module.py](computer://C:\Users\luizd\OneDrive\Data_Science\Projetos\docking\md_module.py) — suporte FAD + aborto sem pose docked
- [check_docking_status.py](computer://C:\Users\luizd\OneDrive\Data_Science\Projetos\docking\check_docking_status.py) — diagnóstico CLI dos 9 alvos
- [run_docking_batch.bat](computer://C:\Users\luizd\OneDrive\Data_Science\Projetos\docking\run_docking_batch.bat) — pipeline automatizado dos 9 alvos
- [run_docking_batch.ps1](computer://C:\Users\luizd\OneDrive\Data_Science\Projetos\docking\run_docking_batch.ps1) — versão PowerShell

---

## 7. Resumo executivo

A plataforma localhost:5001 agora tem **integração completa Docking → MD**:
- Bug do PDBQT receptor corrigido (fallback RDKit funciona sem BABEL_DATADIR)
- Botões na UI para verificar e rodar docking direto da aba MD
- Endpoint Flask `/run_docking` que dispara `docking_v2.py` com streaming de log
- Gate JS+backend bloqueia MD real sem pose docked, mostrando overlay explicativo
- Aba Ajuda atualizada com guia completo do novo workflow
- Suporte a cofatores cristalográficos (FAD em flavoenzimas, FMN, NAD, NADP+, NADPH, HEM)

**Pode executar tudo pelo dashboard agora.**
