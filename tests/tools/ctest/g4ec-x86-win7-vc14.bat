@echo off  
rem ---Compiler------------------------------------------------------
call "%VS140COMNTOOLS%\..\..\vc\vcvarsall"
rem ---Xerces-C------------------------------------------------------

set xerc_dirs=D:\sw\lcg\external\XercesC\3.1.2\x86-windows-vc14
set xerc_dirs=D:\sw\lcg\external\XercesC\3.1.1\x86-windows-vc10,%xerc_dirs%
set xerc_dirs=C:\build\sw\lcg\external\XercesC\3.1.1\x86-windows-vc10,%xerc_dirs%

for %%a IN (%xerc_dirs%) DO (
	if exist %%a   set XERCESC_ROOT_DIR=%%a
)         
echo Using XercesC from %XERCESC_ROOT_DIR%

set PATH=%XERCESC_ROOT_DIR%\bin;%PATH%

rem ---------------------------------------

set CONFIG=x86-win7-vc14

rem - location of this .bat file
set THIS=%~d0%~p0

rem ---Define basic config parameters--------------------------------
call %THIS%g4-win-common.bat

echo Geant4 CMake options - 2 : %G4_XOPTS%

rem ---Optionally remove source code tree
if "%RESET_SOURCE%"=="yes" (
   echo "Removing source code tree %SOURCE%"
   if EXIST %SOURCE% (
      echo "Removing source code tree %SOURCE%"
	   rmdir /q/s "%SOURCE%"
	)
)

rem - echo %Path%
rem - python --version 
rem - set
set Path=%Path%;C:\Python27
rem ---Checkout Geant4 sources ---------------------------------------
if NOT EXIST %SOURCE% (
  python %THIS%g4tagsvn.py checkout -c %VERSION% -d %SOURCE% -q
)
rem ---Run the CTest script-------------------------------------------
ctest -V -S %THIS%g4%MODE%.cmake 
