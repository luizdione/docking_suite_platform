# Adendo ao relatório de desenvolvimento — UI de Subprojetos implementada e verificada

**Complementa:** `RELATORIO_DESENVOLVIMENTO_PLATAFORMAS.md` (mesmo run, 2026-06-19).
**Motivo:** o relatório principal listou como *próximo passo prioritário* a UI de
subprojetos na docking e registrou que "os apps grandes não foram editados". Este adendo
documenta que esse passo foi, de fato, **executado e verificado** nesta mesma execução,
de forma **aditiva e não-destrutiva**.

## O que foi feito em `dashboard_v2.py` (docking, Flask)

Rotas Flask **aditivas** (não alteram nenhuma aba existente), inseridas imediatamente antes
do `ENTRY POINT`, consumindo o layer **autoritativo** `subprojeto/`:

| Rota | Função |
|------|--------|
| `GET /subprojects` | Página-resumo do guarda-chuva: por subprojeto, SVG do pipeline + abstract (`RESUMO_<ID>.md`) + tabela de scripts + etapas com status. |
| `GET /api/subprojects` | JSON: programa (NOME GERAL) + lista de subprojetos com progresso (etapas done/total). |
| `GET /api/subproject/<ID>` | JSON: `config` (de `subprojeto_config.json`) + `status` (de `status.json`). |
| `GET /api/subproject/<ID>/pipeline.svg` | SVG do pipeline; **regenera sob demanda** via `subprojeto/pipeline_svg.py` se o arquivo não existir. |

Além disso, um botão **"☂ Subprojetos"** foi adicionado à barra de abas (link
`target="_blank"` para `/subprojects` — não exige fiação no `switchTab`, minimizando risco).

Helpers adicionados: `_load_subproj_conf`, `_subproj_status`, `_subproj_progress`,
`_ensure_svg`. Nenhuma dependência nova.

## Verificação (demanda 10 — testes em 2º plano)

Servidor Flask iniciado em porta de teste, com `curl` em todas as rotas:

```
GET /api/subprojects                      -> HTTP 200  (SUB01; progress 0.545 = 6/11 etapas)
GET /api/subproject/SUB01/pipeline.svg    -> HTTP 200  (SVG viewBox 0 0 1500 760)
GET /subprojects                          -> HTTP 200
GET /data_status   (rota pré-existente)   -> HTTP 200  (intacta)
GET /              (dashboard principal)  -> HTTP 200  (intacto)
py_compile dashboard_v2.py                -> OK
```

## Por que isto é seguro

- As rotas são **somente-leitura** sobre `subprojeto/` e regeneram SVG via o gerador já
  existente. Não tocam nas abas, no estado de execução, nem nos dados científicos.
- A edição em `dashboard_v2.py` foi cirúrgica (um bloco antes do `__main__` + um link na
  barra de abas) e validada por `py_compile` **e** por servidor vivo.

## Outras entregas deste adendo

- **Ajuda (demanda 8):** `docs/SUBPROJETOS_GUARDA_CHUVA.md` — explica o padrão e os endpoints;
  surge automaticamente na aba **Ajuda** (que lista `docs/*.md`).
- **Pipeline SVG (demanda 2):** regenerados em ambas as plataformas
  (`docking_project/figures/pipeline_sub01.svg`; `epigenoma/www/pipeline_epigenoma.svg`).

## Correção de rumo

No início deste run, a plataforma **LILO/IA.agente** foi erroneamente tratada como 2ª
plataforma e recebeu um modelo divergente de subprojetos; **tudo foi revertido** (LILO
voltou ao git limpo em `arch-v10`) e os artefatos divergentes foram para
`OBSOLETOS/divergent_subprojects_2026-06-19/`. As duas plataformas canônicas são
**docking + epigenoma**, sob a convenção autoritativa `subprojeto/` (subNN maiúsculo).
