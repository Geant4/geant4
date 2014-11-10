@echo off  
rem ---Compiler------------------------------------------------------
call "%VS100COMNTOOLS%\..\..\vc\vcvarsall" amd64
rem ---Xerces-C------------------------------------------------------
set XERCESC_ROOT_DIR=D:\sw\lcg\external\XercesC\3.1.1\x86_64-windows-vc10
set PATH=%XERCESC_ROOT_DIR%\bin;%PATH%

rem ---Define basic config parameters--------------------------------
set THIS=%~d0%~p0
if NOT DEFINED MODE    set MODE=nightly
if NOT DEFINED WORKDIR set WORKDIR=D:/build/cdash
if NOT DEFINED VERSION set VERSION=g4tags-dev
set CONFIG=x64-win7-vc10
set SOURCE=%WORKDIR%/%MODE%/%CONFIG%/%VERSION%
set BINARY=%WORKDIR%/%MODE%/%CONFIG%/build

rem ---Extra Options--------------------------------------------------
set G4_XOPTS="-DGEANT4_USE_OPENGL_WIN32=ON"

rem ---Optionally remove source code tree
if "%RESET_SOURCE%"=="yes" (
   echo "Removing source code tree %SOURCE%"
   if EXIST %SOURCE% (
      echo "Removing source code tree %SOURCE%"
	   rmdir /q/s "%SOURCE%"
	)
)

rem ---Checkout Geant4 sources ---------------------------------------
if NOT EXIST %SOURCE% (
  python %THIS%g4tagsvn.py checkout -c %VERSION% -d %SOURCE% -q
)
rem ---Run the CTest script-------------------------------------------
ctest -V -S %THIS%g4%MODE%.cmake 
