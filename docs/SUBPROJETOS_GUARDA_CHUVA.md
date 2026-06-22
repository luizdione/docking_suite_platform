# Padrão guarda-chuva de subprojetos (NQ-Antiparasitários)

Este documento descreve como a plataforma **Docking Suite** organiza vários
manuscritos/estudos sob um mesmo programa ("guarda-chuva"), em harmonia com a
plataforma irmã **epigenoma** (R/Shiny).

## Conceito

- **Programa (NOME GERAL):** o guarda-chuva — aqui, *NQ-Antiparasitários: descoberta
  computacional de naftoquinonas antiparasitárias (docking + MD + ADMET)*.
- **Subprojeto (subNN):** cada estudo recebe um identificador estável `SUB01`, `SUB02`,
  `SUB03`… O **nome específico** de cada subprojeto é o **título atual do paper**.
- Subprojetos podem ser **independentes** ou **dependentes** entre si.

## Onde ficam os dados (fonte de verdade)

```
subprojeto/
  subprojeto_config.json     # registry: programa + cada SUBNN (sigla, título, alvos, scripts, caminhos)
  pipeline_svg.py            # gera o SVG do pipeline por subprojeto (espelha o gerador da epigenoma)
  SUB01/
    RESUMO_SUB01.md          # abstract-âncora (sincronizado com a aba Produção de Texto)
    status.json              # etapas (done/pending/blocked/data) + métricas-chave
    pipeline_sub01.svg       # cópia do diagrama para submissão
```

O diagrama servido ao dashboard fica em `docking_project/figures/pipeline_<id>.svg`.

## Como ver na plataforma

Na barra de abas do dashboard, clique em **☂ Subprojetos** (ou acesse `/subprojects`).
A página mostra, para cada subprojeto:

- o **pipeline em SVG** com o status de cada etapa (verde = concluído, amarelo = a fazer,
  cinza = bloqueado/futuro, azul = dados);
- o **resumo/abstract** (lido de `RESUMO_<ID>.md`);
- os **scripts** usados/a usar (etapa → arquivo → engine/saída);
- as **etapas** com status individual.

### Endpoints (uso programático)

| Rota | Retorno |
|------|---------|
| `GET /subprojects` | página HTML resumo do guarda-chuva |
| `GET /api/subprojects` | JSON: programa + lista de subprojetos com progresso |
| `GET /api/subproject/<ID>` | JSON: config + status do subprojeto |
| `GET /api/subproject/<ID>/pipeline.svg` | SVG do pipeline (regenera sob demanda) |

## Criar um novo subprojeto

1. Crie a pasta `subprojeto/SUB0N/` com `RESUMO_SUB0N.md` e `status.json`.
2. Adicione o bloco `SUB0N` em `subprojeto/subprojeto_config.json` (sigla, título atual do
   paper, alvos, scripts, caminhos).
3. Rode `python subprojeto/pipeline_svg.py SUB0N` para gerar o SVG.
4. O dashboard passa a exibir o novo subprojeto automaticamente.

## Atualizar o status / pipeline

Edite `subprojeto/SUB0N/status.json` (campo `etapas`, cada uma com `status`
`done`/`pending`/`blocked`/`data`) e regenere o SVG com `pipeline_svg.py`. O resumo/abstract
deve refletir sempre a versão mais atual do paper (harmonia com a skill `/build_paper`).
