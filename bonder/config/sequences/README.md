# Sequencias das iscas (FASTA)

As iscas do piloto sao **exogenas** ao proteoma murino e precisam entrar no FASTA
de busca (`search_db.fasta`) para serem quantificadas como controle positivo.
Coloque aqui **um FASTA por isca**, com o nome apontado em `config/baits.yaml`.

## 1. Galectina-3 humana (LGALS3, `LGALS3_HUMAN`) — UniProt P17931

Baixe a sequencia canonica exata (nao digite a mao):

```bash
curl -L "https://rest.uniprot.org/uniprotkb/P17931.fasta" -o LGALS3_HUMAN.fasta
```

## 2. NSs do ortobunyavirus Icoaraci (`NSs_ICOARACI`)

Sequencia especifica do seu isolado. Salve em `NSs_Icoaraci.fasta` com um
cabecalho que combine com `search_db_header` em `baits.yaml`, por exemplo:

```
>bait|NSs_ICOARACI|NSs Icoaraci orthobunyavirus
MSEQ...   (a sequencia de aminoacidos do seu isolado)
```

Se houver acesso no GenBank/UniProt para o seu isolado, baixe de la; caso
contrario, use a traducao do ORF de NSs do seu sequenciamento.

## Verificacao

Depois de adicionar os dois arquivos:

```bash
bonder baits check
```

Isso valida as iscas e gera `data/search_db.fasta` (proteoma murino +
contaminantes + iscas) para a busca. Sem os FASTAs, o comando apenas avisa
(nao quebra) e o restante do pipeline continua funcionando a partir de uma
saida de busca ja pronta.
