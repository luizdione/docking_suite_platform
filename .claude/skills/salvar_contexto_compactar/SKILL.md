---
name: salvar_contexto_compactar
description: >-
  Salva a memória da sessão nos arquivos Markdown do repositório (igual à skill
  salvar_contexto) e, em seguida, prepara a REDUÇÃO de contexto mantendo
  continuidade: entrega ao usuário o comando /compact pronto (com foco) para ele
  digitar. Use quando o contexto estiver grande (perto do limite / 1M no Opus) e
  a intenção for "salvar E compactar" numa tacada só, ou quando o usuário pedir
  "salvar e compactar", "reduzir contexto", "compactar mantendo a memória".
---

# /salvar_contexto_compactar — Salvar memória + preparar a compactação

Faz o mesmo que `/salvar_contexto` (grava a memória da sessão nos `.md` do repo) e,
no final, **prepara a redução de contexto**: monta e entrega o comando `/compact`
com foco, para o usuário digitar na sequência — reduzindo o contexto **sem perder
a continuidade**, já que todo o estado ficou salvo em disco.

## ⚠️ Limite importante (ler primeiro)
Uma skill **não consegue executar `/compact` sozinha** — `/compact`, `/clear` e
`/branch` são comandos do próprio Claude Code, digitados pelo usuário; não há
ferramenta que os dispare de dentro de uma skill. Portanto esta skill:
1. **automatiza** o que dá: gravar/atualizar os arquivos de memória; e
2. **instrui** o usuário a digitar o `/compact` (comando já pronto no fim da saída).

Nunca afirme que "o contexto foi compactado" — quem compacta é o usuário ao digitar
`/compact`. A skill deixa tudo pronto para isso.

## Quando usar
- Contexto grande / perto do limite (ex.: perto de 1M no Opus) e quer **salvar e já
  compactar** em um passo só.
- Ao encerrar um bloco de trabalho longo antes de continuar na mesma sessão.

O usuário pode passar um argumento livre com o foco/tema (ex.:
`/salvar_contexto_compactar módulo ADMET`). Use-o como tag do snapshot **e** como
foco sugerido no comando `/compact` final.

## Passo a passo

1. **Gravar a memória** — execute exatamente o mesmo processo da skill
   `salvar_contexto` (mesmos arquivos, mesmo template). A especificação completa e o
   template estão em **`.claude/skills/salvar_contexto/SKILL.md`** — siga-a na
   íntegra. Em resumo, isso atualiza:
   - `docs/historico/CONTEXTO_SESSAO.md` (arquivo vivo, sobrescrito);
   - `docs/historico/BRIEFING_<AAAA-MM-DD>_<tag>.md` (snapshot datado, novo);
   - a seção `## Instruções de Compactação` no `claude.md` da raiz (só se faltar).

2. **Confirmar a gravação** — verifique que os arquivos foram escritos (liste os
   caminhos). Se algo falhar, **pare e avise** — não oriente a compactar antes de a
   memória estar salva (senão o contexto é perdido).

3. **Montar o foco da compactação** — derive 1 linha de foco a partir do objetivo e
   das pendências salvos (ou do argumento do usuário). Ex.: `foque no módulo ADMET e
   nos próximos passos pendentes`.

4. **Entregar o comando ao usuário** — termine a resposta com um bloco claro e
   destacado, por exemplo:

   > ✅ Memória salva em `docs/historico/CONTEXTO_SESSAO.md`.
   > **Agora digite para reduzir o contexto:**
   > ```
   > /compact <foco derivado no passo 3>
   > ```
   > Depois do `/compact`, o `claude.md` é reinjetado automaticamente. Para retomar
   > o detalhe, peça: **"releia `docs/historico/CONTEXTO_SESSAO.md`"**.

## Alternativas a oferecer (se fizer sentido)
- **Começar limpo** mantendo o histórico resumível: `/clear` e depois `/resume`.
- **Ramificar** sem perder o original: `/branch`.
  Em ambos, a memória já salva em `docs/historico/CONTEXTO_SESSAO.md` garante a
  continuidade — basta reler esse arquivo ao retomar.

## Regras importantes
- **NÃO commitar nem fazer push** automaticamente (política de commits do repo — ver
  `claude.md`). Só gravar/atualizar os arquivos locais.
- **Não** copiar segredos (tokens Notion/Anthropic, credenciais) para os `.md`.
- Nunca dizer que compactou: a skill **prepara**; o usuário **executa** o `/compact`.
