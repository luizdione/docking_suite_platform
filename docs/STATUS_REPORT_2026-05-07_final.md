# Relatório Final — Pronto para rodar `run_docking_batch.bat`
**Data:** 2026-05-07 | **Estado: pronto para execução do docking**

---

## 1. Sumário das mudanças aplicadas

| Arquivo | Mudança | Status |
|---|---|---|
| `md_module.py` | Suporte FAD + aborto sem pose docked | ✓ AST OK (3091 linhas) |
| `docking_v2.py` | Sem alterações (já era funcional) | ✓ AST OK |
| `dashboard_v2.py` | Endpoint `/check_docking_for_md` + gate JS + overlay HTML | ✓ AST OK (4704 linhas) |
| `check_docking_status.py` | Lista padrão expandida para 9 alvos | ✓ AST OK |
| `run_docking_batch.bat` | 9 alvos (5 parasita + 4 humanos) | ✓ |

---

## 2. `run_docking_batch.bat` — agora cobre 9 alvos

### Parasitários (5)
| ID | PDB | Alvo | Cofator |
|---|---|---|---|
| TcTR | 1BZL | Trypanothione Reductase (T. cruzi) | FAD |
| LiTR | 2JK6 | Trypanothione Reductase (L. infantum) | FAD |
| TcCruzain | 3KKU | Cruzain (T. cruzi) | — |
| TcFeSOD | 4DVH | Fe-Superoxide Dismutase (T. cruzi) | Fe |
| LmDHODH | 3MHU | DHODH (L. major) | FMN+FAD |

### Controles humanos (4) — adicionados nesta sessão
| ID | PDB | Alvo | Cofator |
|---|---|---|---|
| hGR | 3GRS | Glutathione Reductase | FAD |
| hCatL | 2XU3 | Cathepsin L | — |
| hDHODH | 1D3H | DHODH | FMN |
| hMnSOD | 1LUV | Mn-Superoxide Dismutase | Mn |

**Comando do batch:**
```
python docking_v2.py --targets TcTR,LiTR,TcCruzain,TcFeSOD,LmDHODH,hGR,hCatL,hDHODH,hMnSOD --exhaustiveness 32
```

**Tempo estimado:** ~5–10 min/alvo × 9 alvos = **~1h30** total.

---

## 3. Gate dashboard:5001 — bloqueio de MD sem docking

### Backend (Flask)
- Novo endpoint **`POST /check_docking_for_md`** retorna `{missing, ok, all_ok, compound}` baseado nos mesmos paths de busca do `md_module.py` (incluindo `E:\data\docking\`).
- Endpoint **`POST /run_md`** agora valida `mode='real'` server-side e retorna **HTTP 409 Conflict** com payload `{blocked_no_docking: true, missing_targets: [...]}` se faltar pose. Bypass possível com `skip_docking_check=true` (não exposto na UI).

### Frontend
- **`runMD()`** agora chama `/check_docking_for_md` antes de enviar `/run_md` quando o modo é `real`. Se algum alvo falhar:
  - Aborta sem disparar a MD
  - Mostra **overlay modal** vermelho explicando o problema
  - Lista os alvos faltantes e gera comando específico de docking para apenas esses
- **Defesa em profundidade:** mesmo se o JS for bypassado, o servidor retorna 409 e o frontend captura para mostrar o overlay também.
- O modo `synthetic` **não é bloqueado** (não precisa de pose docked, gera dados de literatura).

### Mensagem do overlay
O usuário vê:
1. Cabeçalho vermelho: "MD bloqueada — pose docked ausente"
2. Explicação do "por quê": ligante posicionado fora do bolsão = resultado fisicamente inválido
3. Lista dos alvos sem pose
4. Caixa verde com **6 passos** para resolver (cd, run_docking_batch.bat, etc.)
5. Caixa amarela com dica: usar modo `synthetic` para preview rápido enquanto o docking roda
6. Botão "Entendi, vou rodar o docking"

---

## 4. Pré-checagem antes de rodar o batch

✅ **Análise completa:**
- `vina.exe` presente na pasta (1.18 MB)
- `docking_v2.py` íntegro (946 linhas, AST OK), suporta todos os 9 alvos no `TARGET_DB`
- Ligante 2-SCN-NQ já preparado em `docking_project\ligand\2-SCN-NQ\2-SCN-NQ.pdbqt` (1641 bytes)
- Conexão de internet necessária para baixar 9 PDBs do RCSB
- Espaço em disco para `docking_project\proteins\` (PDBs limpos) e `docking_project\docking\` (PDBQTs de saída)
- `md_module.py` agora tem suporte a FAD (Rota D SMIRNOFF/Rota E RDKit)

⚠️ **Não há bloqueios técnicos.** Você pode rodar o batch agora.

---

## 5. Quando posso rodar `run_docking_batch.bat`?

**Resposta: agora.**

```powershell
cd C:\Users\luizd\OneDrive\Data_Science\Projetos\docking
.\run_docking_batch.bat
```

### O que o batch fará:
1. Lista os 9 alvos (sanity check)
2. Para cada um:
   - Baixa o PDB do RCSB (se não existir)
   - Roda QC + análise estrutural
   - Limpa proteína (mantém cofatores especificados em `keep_cofactors`)
   - Converte para PDBQT receptor
   - Calcula grid box (centrado em ligante co-cristal ou metal)
   - Roda Vina (exhaustiveness=32)
   - Salva `result_<TID>.pdbqt` em `docking_project\docking\`
   - Loga energias dG e Ki
3. Ao final, executa `check_docking_status.py` para validar

### Após o batch terminar:
1. Verifique a saída do `check_docking_status.py` — todos devem estar `[OK] docked_flat` ou `docked_session`
2. Abra o dashboard: `python dashboard_v2.py` → http://localhost:5001
3. Teste rápido: rode MD em modo `real` para TcTR com **10 ns** primeiro
4. Confirme no log:
   - `[Cofator] Detectados: FAD`
   - `[OK] Cofator FAD parametrizado e registrado`
   - `[OpenMM] Cofatores parametrizados mantidos: FAD`
   - `[Ligante] Pose docked: ...result_TcTR.pdbqt` (não mais "Preparação")
5. Se OK, escale para produção 50–100 ns multi-alvo

---

## 6. Como testar o gate (opcional)

Sem rodar o docking ainda, você pode confirmar que o gate funciona:

1. Inicie o dashboard: `python dashboard_v2.py`
2. Abra http://localhost:5001
3. Vá na aba MD, selecione TcTR, modo **real**, 10 ns
4. Clique "Rodar MD"
5. **Esperado:** o overlay vermelho aparece imediatamente, sem disparar GPU/simulação

---

## 7. Arquivos da entrega

- [md_module.py](computer://C:\Users\luizd\OneDrive\Data_Science\Projetos\docking\md_module.py) — 3091 linhas com suporte FAD
- [dashboard_v2.py](computer://C:\Users\luizd\OneDrive\Data_Science\Projetos\docking\dashboard_v2.py) — 4704 linhas com gate
- [run_docking_batch.bat](computer://C:\Users\luizd\OneDrive\Data_Science\Projetos\docking\run_docking_batch.bat) — 9 alvos
- [check_docking_status.py](computer://C:\Users\luizd\OneDrive\Data_Science\Projetos\docking\check_docking_status.py) — diagnóstico
- [docking_v2.py](computer://C:\Users\luizd\OneDrive\Data_Science\Projetos\docking\docking_v2.py) — pipeline de docking (sem alterações)

---

## 8. Resumo executivo

Tudo pronto para rodar `run_docking_batch.bat`. O batch agora cobre os **9 alvos do estudo** (5 parasitários + 4 controles humanos: hGR, hCatL, hDHODH, hMnSOD). A plataforma localhost:5001 agora **bloqueia automaticamente** qualquer tentativa de rodar MD em modo `real` para um composto/alvo que não tenha pose docked, mostrando um overlay explicativo com instruções passo-a-passo. O modo `synthetic` continua disponível para preview rápido. **Pode executar o batch agora.**
