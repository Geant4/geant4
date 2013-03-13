@echo off  
echo Execution started: %date% %time%

rem ---Compiler------------------------------------------------------
call "%VS100COMNTOOLS%vsvars32.bat"
rem ---Xerces-C------------------------------------------------------
set XERCESC_ROOT_DIR=D:\sw\lcg\external\XercesC\3.1.0\i686-winxp-vc9-opt
set PATH=%XERCESC_ROOT_DIR%\bin;%PATH%

rem ---Define basic config parameters--------------------------------
set THIS=%~d0%~p0
set MODE=nightly
set WORKDIR=D:/build/cdash/%MODE%
set VERSION=g4tags-dev
set CONFIG=win7-vc10
set SOURCE=%WORKDIR%/%VERSION%
set BINARY=%WORKDIR%/%CONFIG%

rem rmdir /q/s "%SOURCE%"
rem ---Checkout Geant4 sources ---------------------------------------
python %THIS%g4tagsvn.py update -c %VERSION% -d %SOURCE% --quiet

rem ---Run the CTest script-------------------------------------------
ctest -V -S %THIS%g4%MODE%.cmake 

echo Execution ended: %date% %time%

