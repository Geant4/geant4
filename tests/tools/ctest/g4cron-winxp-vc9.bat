@echo off  
echo Execution started: %date% %time%
rem ---Compiler------------------------------------------------------
call "c:\Program Files\Microsoft Visual Studio 9.0\Common7\Tools\vsvars32.bat"
rem ---Xerces-C------------------------------------------------------
set XERCESC_ROOT_DIR=E:/local/lib/lcg/external/XercesC/3.1.1p1/i686-winxp-vc9-opt
set PATH=%XERCESC_ROOT_DIR%\bin;%PATH%

set THIS=%~d0%~p0
set MODE=nightly

set WORKDIR=F:/CDash/G4
set VERSION=g4tags-dev
set CONFIG=winxp-vc9

set SOURCE=%WORKDIR%/%CONFIG%-src
set BINARY=%WORKDIR%/%CONFIG%

rem ---Load the RSA key ---------------------------------------------
start /b "C:\Program Files\PuTTY\pageant.exe" \\cern.ch\dfs\Users\m\mato\Documents\.ssh\rsakey.ppk

rmdir /q/s "%SOURCE%"
rem ---Checkout Geant4 sources ---------------------------------------
python %THIS%g4checkout.py -c %VERSION% -d %SOURCE% --quiet


rem ---Run the CTest script-------------------------------------------
ctest -V -S %THIS%g4nightly.cmake 

echo Execution ended: %date% %time%

