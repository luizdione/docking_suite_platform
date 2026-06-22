# auto_fix — pacote offline de correcao do MD do TcTR

Roda **localmente, sem LLM**, aplicando correcoes rule-based a partir do diagnostico em [FAD_XML_SCRAMBLE_DIAGNOSIS.md](../FAD_XML_SCRAMBLE_DIAGNOSIS.md). Voce starta e vai embora — quando voltar, le o `auto_md_report.md` com tudo.

## Arquivos

| Script | O que faz |
|---|---|
| `fix_fad_xml.py`     | Corrige `cofactor_FAD_TcTR.xml`: elementos errados + bonds faltando em HBF (rule-based, GAFF-like defaults) |
| `validate_system.py` | Sanity check headless: `md_clean.py` AST, XML estrutura, OpenMM cria System sem rodar steps |
| `auto_md_loop.py`    | Loop principal: valida → fix → roda MD curto (0.1 ns) → classifica log → aplica novo patch → repete |

## Uso (rodar antes de sair, ler quando voltar)

```powershell
conda activate datascience
cd "C:\Users\luizd\OneDrive\Data_Science\Projetos\docking"

# Modo padrao — 8 iteracoes max, MD de 0.1 ns por tentativa
python auto_fix\auto_md_loop.py

# Versao agressiva — mais tempo por tentativa
python auto_fix\auto_md_loop.py --max-iter 12 --ns 0.2 --max-run-sec 1800

# So validar (nao roda MD) — bem rapido
python auto_fix\validate_system.py

# Aplicar so o fix do XML
python auto_fix\fix_fad_xml.py --inplace
```

## Outputs (gerados na raiz `docking/`)

| Arquivo | Conteudo |
|---|---|
| `auto_md_report.md` | Relatorio markdown — uma secao por iteracao com classe do erro e tail do log |
| `auto_md_state.json` | Estado estruturado (JSON) — facil de parsear |
| `auto_md_logs/` | Logs completos por iteracao (`iter_NN_validate.log`, `iter_NN_md.log`, etc.) |

## Classes de erro reconhecidas

| Classe | Patch aplicado |
|---|---|
| `NAN` / `BOND_ANOMALY_THEN_NAN` | (so a 1a iter) fix_fad_xml.py — corrige elementos + adiciona bonds harmonicas |
| `MISSING_TEMPLATE` | (ainda sem fix automatico — registra para inspecao humana) |
| `SYNTAX_IMPORT`    | (ainda sem fix automatico) |
| `TIMEOUT`          | (sem patch; vc pode aumentar --max-run-sec) |
| `RUNTIME_ERROR`    | (sem fix; loga a excecao) |
| `SUCCESS`          | Loop para. PDB de interacao TcTR-FAD-LIG no diretorio do tag |

Stop conditions:
- 3 iteracoes consecutivas com o **mesmo** erro → para (loop detectado)
- `--max-iter` atingido → para
- `SUCCESS` → para com codigo 0

## O que ja foi diagnosticado (nao precisa redescobrir)

1. **33/86 atom types no FAD XML tem element errado** — N3A diz "O", C4A diz "P", etc.
2. **`HarmonicBondForce` esta vazio (0 bonds)** — 91 bonds da Residue sem termo harmonico
3. Sintoma: warm-up explode em NaN imediatamente porque as bonds tem k=0

O `fix_fad_xml.py` resolve ambos. Se isso nao for suficiente para warm-up sobreviver, o proximo bottleneck e:
- Falta de `<HarmonicAngleForce>` (angle bending sem parametros)
- Falta de `<PeriodicTorsionForce>` (torsions sem parametros)
- Cargas nao definidas em `<NonbondedForce>` (chargesn proton)

Esses tres NAO sao patchados automaticamente pelo `fix_fad_xml.py`. Se o relatorio mostrar que apos o fix do XML ainda explode NaN, o proximo passo manual e gerar o XML completo via OpenFF Toolkit (ver Opcao A em `FAD_XML_SCRAMBLE_DIAGNOSIS.md`).

## Politica de tamanho do md_clean.py

Atual: 858 linhas (limite 1500). Se um patch precisar passar de 1500, **NAO inline** — refatorar em modulo separado.

`validate_system.py` ja verifica isso e falha se ultrapassar.
