@echo off  
rem ---Compiler------------------------------------------------------
call "%VS90COMNTOOLS%vsvars32.bat"
rem ---Xerces-C------------------------------------------------------
set XERCESC_ROOT_DIR=D:\sw\lcg\external\XercesC\3.1.0\i686-winxp-vc9-opt
set PATH=%XERCESC_ROOT_DIR%\bin;%PATH%

rem ---Define basic config parameters--------------------------------
set THIS=%~d0%~p0
if NOT DEFINED MODE    set MODE=nightly
if NOT DEFINED WORKDIR set WORKDIR=D:/build/cdash
if NOT DEFINED VERSION set VERSION=g4tags-dev
set CONFIG=winxp-vc9
set SOURCE=%WORKDIR%/%MODE%/%CONFIG%/%VERSION%
set BINARY=%WORKDIR%/%MODE%/%CONFIG%/build

rem ---Checkout Geant4 sources ---------------------------------------
if NOT EXIST %SOURCE% python %THIS%g4tagsvn.py update -c %VERSION% -d %SOURCE% --quiet

rem ---Optionally remove source code tree
if "%RESET_SOURCE%"=="yes" (
   echo "Removing source code tree %SOURCE%"
   if EXIST %SOURCE% (
      echo "Removing source code tree %SOURCE%"
	   rmdir /q/s "%SOURCE%"
	)
)

rem ---Checkout Geant4 sources ---------------------------------------
if NOT EXIST %SOURCE% python %THIS%g4tagsvn.py update -c %VERSION% -d %SOURCE% --quiet

rem ---Delete libraries and executables (important for continuous)----
if EXIST %BINARY%/outputs rmdir /s/q %BINARY%/outputs

rem ---Run the CTest script-------------------------------------------
ctest -V -S %THIS%g4%MODE%.cmake 


