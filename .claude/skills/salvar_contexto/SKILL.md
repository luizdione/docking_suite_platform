---
name: salvar_contexto
description: >-
  Gera e atualiza a memória da sessão atual (objetivo, progresso, decisões,
  arquivos tocados e pendências) em arquivos Markdown versionados no repositório,
  para permitir reduzir o contexto com /compact, /clear ou /branch sem perder a
  continuidade da tarefa. Use quando o contexto estiver grande (perto do limite /
  1M no Opus), antes de compactar/limpar/forkar a sessão, ao encerrar um bloco de
  trabalho, ou sempre que o usuário pedir para "salvar contexto", "salvar a
  sessão", "gerar memória" ou "fazer handoff".
---

# /salvar_contexto — Memória de sessão e handoff

Automatiza o fluxo de **preservar a memória da sessão** antes de reduzir o
contexto. O objetivo é permitir que o usuário rode `/compact`, `/clear` (e depois
`/resume`) ou `/branch` **sem perder o estado da tarefa em andamento**, porque
todo o estado relevante fica gravado em arquivos `.md` no próprio repositório —
que sobrevivem a compactação, clear e troca de sessão.

## Quando usar
- O contexto está grande / perto do limite (ex.: perto de 1M no Opus).
- Antes de `/compact`, `/clear` ou `/branch`.
- Ao terminar uma etapa importante e querer um ponto de retomada.
- Quando o usuário pedir explicitamente (`/salvar_contexto`).

O usuário pode passar um argumento livre com o foco/nota da sessão, por exemplo
`/salvar_contexto revisão do módulo ADMET`. Se houver argumento, use-o como
**tag/tema** do snapshot e destaque-o no resumo.

## O que a skill produz
1. **Arquivo vivo** `docs/historico/CONTEXTO_SESSAO.md` — sempre reescrito com o
   estado **mais recente** (é o handoff que você lê ao retomar).
2. **Snapshot datado** `docs/historico/BRIEFING_<AAAA-MM-DD>_<tag>.md` — histórico
   imutável daquela sessão (segue a convenção já existente na pasta). Não
   sobrescreve snapshots anteriores.
3. **Ponteiro de compactação** no `claude.md` da raiz — uma seção
   `## Instruções de Compactação` (criada só se ainda não existir) que instrui a
   preservar o conteúdo do arquivo vivo ao compactar.

> Localização: por padrão em `docs/historico/`. Se o repositório não tiver essa
> pasta, criar. Sempre resolver a **raiz do projeto** com
> `git rev-parse --show-toplevel` e gravar relativo a ela.

## Passo a passo (executar nesta ordem)

1. **Descobrir a raiz e coletar metadados** com um único bloco de shell:
   - Raiz: `git rev-parse --show-toplevel`
   - Data/hora: `date +'%Y-%m-%d %H:%M %z'` e `date +%Y-%m-%d`
   - Branch: `git rev-parse --abbrev-ref HEAD`
   - Estado do working tree: `git status --short`
   - Resumo das mudanças não commitadas: `git diff --stat` e `git diff --stat --cached`
   - Últimos commits: `git log --oneline -10`
   - Garantir a pasta: `mkdir -p <raiz>/docs/historico`

2. **Montar o conteúdo do documento** a partir de DUAS fontes combinadas:
   - **O histórico desta conversa** (a fonte principal do "porquê"): objetivo do
     usuário, decisões tomadas, caminhos descartados, pendências, próximos passos.
   - **Os metadados do git** coletados acima (o "o quê" concreto: branch,
     arquivos, commits).

   Use o **template** abaixo. Seja específico e conciso — sem preenchimento
   artificial. Se alguma seção não se aplicar, escreva `- (nada nesta sessão)`.

3. **Escrever o arquivo vivo** `docs/historico/CONTEXTO_SESSAO.md`
   (sobrescrever por completo — ele reflete só o estado atual).

4. **Escrever o snapshot datado** `docs/historico/BRIEFING_<AAAA-MM-DD>_<tag>.md`
   com o **mesmo conteúdo**. `<tag>` = slug do argumento do usuário, ou um resumo
   curto do tema (ex.: `contexto_sessao`). Se o arquivo do dia+tag já existir,
   acrescente um sufixo `_UPDATE`, `_UPDATE2`, ... (não sobrescrever histórico).

5. **Garantir a seção de compactação no `claude.md` da raiz** (idempotente):
   - Ler o `claude.md` (na raiz; em Windows o nome é case-insensitive).
   - Se **não** houver uma seção `## Instruções de Compactação`, **acrescentá-la
     ao final** com o bloco padrão abaixo. Se já existir, **não** reescrever (para
     manter o diff do arquivo curado limpo).

     ```markdown
     ## Instruções de Compactação
     Ao compactar (`/compact`) esta sessão, **preserve**:
     - o objetivo/tarefa atual e as decisões já tomadas;
     - os arquivos criados/modificados e os próximos passos pendentes.
     O estado detalhado da sessão fica em `docs/historico/CONTEXTO_SESSAO.md` —
     releia esse arquivo ao retomar após `/compact`, `/clear` ou `/branch`.
     ```

6. **Reportar ao usuário** de forma curta:
   - Quais arquivos foram criados/atualizados (caminhos relativos).
   - Um resumo de 2–4 linhas do que foi salvo.
   - O próximo comando sugerido conforme a intenção:
     - continuar na mesma sessão → **`/compact`** (com foco, ex.: `/compact foque em <tema>`);
     - começar limpo mantendo o histórico resumível → **`/clear`** (depois `/resume`);
     - ramificar sem perder o original → **`/branch`**.
   - Lembrar que, após reduzir o contexto, basta pedir para **reler
     `docs/historico/CONTEXTO_SESSAO.md`** para retomar de onde parou.

## Regras importantes
- **NÃO commitar nem fazer push** automaticamente. Esta skill só grava/atualiza os
  arquivos locais. Deixe o commit a critério do usuário (o repositório tem política
  própria de commits — ver `claude.md`, "Commits SEM trailer de co-autoria de IA").
- **Não** versionar segredos: nunca copiar tokens (Notion/Anthropic), caminhos com
  credenciais ou conteúdo de arquivos gitignored para dentro dos `.md` de memória.
- Nunca inventar estado: se algo não foi feito ou está incerto, marque como
  **pendente** ou **a confirmar**.
- Manter os arquivos enxutos: priorize metadados e decisões, não cole saídas longas
  de ferramentas nem trechos gigantes de arquivos (referencie por caminho:linha).

## Template do documento

```markdown
# Contexto da Sessão — <tema/tag>

- **Data/hora:** <AAAA-MM-DD HH:MM ±ZZZZ>
- **Branch:** <branch>
- **Repositório:** <nome/remote>
- **Foco desta sessão:** <argumento do usuário ou resumo de 1 linha>

## 1. Objetivo / tarefa atual
<O que estamos tentando alcançar, na visão do usuário.>

## 2. Estado / progresso
- <O que já foi feito, em bullets.>
- <Resultados/verificações relevantes.>

## 3. Arquivos criados / modificados
- `caminho/arquivo.py` — <o que mudou e por quê>
- <use `git status`/`git diff --stat`; inclua também não-versionados relevantes>

## 4. Decisões tomadas (e justificativa)
- <Decisão> — <por quê / alternativas descartadas>

## 5. Pendências / próximos passos
- [ ] <próximo passo concreto>
- [ ] <bloqueios ou dúvidas a confirmar com o usuário>

## 6. Como retomar
- Comandos/rotas úteis: <ex.: `python dashboard_v2.py` → http://localhost:5001>
- Ao reabrir a sessão (após /compact, /clear ou /branch), **releia este arquivo**
  e continue pelos itens da seção 5.

## 7. Instruções de compactação (o que preservar)
- Preserve: objetivo (seção 1), decisões (seção 4) e pendências (seção 5).
- Pode descartar: saídas de ferramentas antigas e trechos longos já resumidos aqui.
```
