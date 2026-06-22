# Briefing automatizado — md_module.py v3.2 (2026-05-08)

Tarefa diária retomando do limite de créditos da sessão anterior. Foco em
`md_module.py` da pipeline de docking + dinâmica molecular.

## Diagnóstico recuperado da sessão anterior

A v3 (`2026-05-08.v3.icp-insitu`) tinha duas falhas encadeadas observadas no
último smoke test:

1. **Bug fatal no Embed (Boost.Python ArgumentError)** — `_AC.EmbedMolecule(_mol_h, _AC.ETKDGv3(), randomSeed=42)` mistura a assinatura "params object" (sem kwargs) com a assinatura legada (com kwargs). RDKit não aceita os dois ao mesmo tempo, então a re-adição do FAD foi pulada.
2. **NaN durante NVT fase 1 sub-step 2** — mesmo quando o cofator é pulado, a MD explodiu entre 10 e 20 ps de equilibração @ 1 fs. Causa provável: ligante parametrizado por Rota E (MMFF94 RDKit, marcada PREVIEW) + cargas Gasteiger gerando forças locais ruins no início da dinâmica, sem o FAD que normalmente "preenche" o bolso ativo cristalográfico.

## Fixes aplicados nesta iteração (v3.2)

### 1. Embed do cofator — API correta do RDKit
```python
_params = _AC.ETKDGv3()
_params.randomSeed = 42
_ret = _AC.EmbedMolecule(_mol_h, _params)        # signature 1: params
if _ret == -1:
    _ret = _AC.EmbedMolecule(_mol_h, randomSeed=42, useRandomCoords=True)  # signature 2: kwargs
```
Mais `MMFFOptimizeMolecule` (com fallback UFF) após o Embed para reduzir distorções antes do ICP.

### 2. NVT fase 0 (warm-up conservador) — bloqueia NaN no início
- 5 ps @ dt=0.5 fs e T=100 K (sub-blocos de 1 ps com checagem NaN)
- T-ramp em 4 degraus (150 → 200 → 250 → target K), 1.25 ps cada
- Só depois entra na fase 1 (50 ps @ 1 fs), agora já com sistema acomodado

Resultado esperado: a fase 1 que antes explodia em sub-step 2 agora deve passar limpo, porque o sistema já vem aquecido gradualmente de 0.5 fs.

## Fingerprint de versão atualizada

```python
__MODULE_VERSION__ = "2026-05-08.v3.2.icp-embed-fix+warmup"
```

`smoke_test_fad_fix.py` foi atualizado para validar essa string em runtime — aborta com `sys.exit(8)` se outra versão for importada de cache.

## Validação programática

- `python -c "import ast; ast.parse(open('md_module.py').read())"` — **AST OK**
- `python -c "import ast; ast.parse(open('smoke_test_fad_fix.py').read())"` — **AST OK**

## Próximo passo (ação do usuário)

Em terminal limpo, com TODOS os processos Python fechados (dashboard/jupyter):

```powershell
Get-Process python | Stop-Process -Force -ErrorAction SilentlyContinue
Get-ChildItem -Path C:\Users\luizd -Recurse -Include 'md_module*.pyc' -ErrorAction SilentlyContinue | Remove-Item -Force
conda activate datascience
cd C:\Users\luizd\OneDrive\Data_Science\Projetos\docking
python smoke_test_fad_fix.py
```

Logs esperados se v3.2 funcionar:
```
md_module versao   = 2026-05-08.v3.2.icp-embed-fix+warmup
[OK] Versao correta importada
...
[Cofator->H] FAD #1: ICP RMSD=0.0XX A, 53 heavy + 33 Hs in-situ
[Cofator->H] FAD #2: ICP RMSD=0.0XX A, 53 heavy + 33 Hs in-situ
[OK] 2 cofator(es) re-adicionado(s) ...
[OpenMM] NVT warm-up (fase 0: 5 ps, dt=0.5 fs @ 100K → ramp T) ...
    [OK] warm-up dt=0.5fs estavel @ 300 K
[OpenMM] NVT equilibration (fase 1: 50 ps, dt=1 fs) ...
    [OK] dt=1fs estavel
... (fases 1.5 e 2 sem NaN)
```

Cenários de falha residuais a observar:
- ICP RMSD > 1 Å — bug de mapeamento por elemento+geometria; fixar por outro algoritmo (Hungarian / scipy.optimize.linear_sum_assignment).
- NaN reaparece em fase 1.5 ou 2 — ainda problema de parametrização do ligante (Rota E é PREVIEW). Mitigação real: instalar `openff-toolkit` para Rota D ou `ambertools` para GAFF2 (Rotas A/B/C).

## Arquivos modificados

- `md_module.py` — embed fix + warmup phase + fingerprint
- `smoke_test_fad_fix.py` — fingerprint esperado atualizado
