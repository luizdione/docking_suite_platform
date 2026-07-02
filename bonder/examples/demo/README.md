# Demo - analise completa de exemplo

`proteinGroups_demo.txt` e uma saida de busca **sintetica** (estilo MaxQuant,
~250 proteinas, 3 IP x 3 CTRL) cujas colunas casam com o `config/samples.yaml`
padrao. Serve para ver o pipeline COMPLETO (SAINT + DEP/limma + rede) produzindo
resultados reais, sem precisar dos seus dados ainda.

## Rodar o demo (Docker, com SAINT + R de verdade)

```bash
docker compose -f docker/docker-compose.yml run --rm \
  -v "${PWD}/examples:/examples" worker \
  bonder run --input /examples/demo/proteinGroups_demo.txt --contrast NSs_vs_CTRL
```

Ou pela interface: copie `examples/demo/proteinGroups_demo.txt` para `inputs/`,
suba com `docker compose ... up`, e no dashboard rode o contraste `NSs_vs_CTRL`.

## O que esperar

- ~25 "interatores verdadeiros" (acessos `PINT0001..PINT0025`) mais o `P12345`,
  `Q8R1` no topo da tabela (alto no IP, baixo no controle);
- a isca `NSs_ICOARACI` e o contaminante `P00761` excluidos dos preys;
- volcano (DEP/limma) e rede de interacao gerados em `runs/<run_id>/`.

Os dados sao ficticios: servem para validar o fluxo e treinar a leitura das
saidas, NAO para conclusoes biologicas.
