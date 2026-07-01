@echo off
REM Cria o atalho (.lnk) da skill /salvar_contexto na pasta Projetos.
REM Duplo-clique para rodar. Chama o criar_atalho.ps1 que esta ao lado.
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0criar_atalho.ps1" %*
pause
