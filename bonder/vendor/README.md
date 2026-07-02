# vendor/ - binarios externos (opcional)

Escotilha para fornecer o **SAINTexpress-int** sem reconstruir a imagem Docker.

O `docker-compose.yml` monta `./vendor/bin` em `/opt/bonder/bin` (que esta no
PATH do conteiner de analise). Se o SAINTexpress nao compilou no build
(`saint: pulado` no autoteste), basta colocar um binario Linux
`SAINTexpress-int` em `vendor/bin/` e reiniciar - o pipeline passa a usa-lo.

Como obter o binario (Linux):

```bash
curl -fL -o saint.tgz "https://sourceforge.net/projects/saint-apms/files/SAINTexpress_v3.6.3__2018-03-09.tar.gz/download"
tar xzf saint.tgz
cd SAINTexpress_v3.6.3__2018-03-09/SAINTexpress-int && bash compile   # requer g++ + libgsl-dev
cp SAINTexpress-int  <repo>/vendor/bin/
```

Binarios nao sao versionados (ver .gitignore).
