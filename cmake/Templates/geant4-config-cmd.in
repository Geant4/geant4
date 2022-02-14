@echo off
rem
rem Self locate script when sourced
setlocal
for /F %%i in ("%~dp0.") do set BINDIR=%%~fi
pushd %BINDIR%
sh geant4-config %*
popd
