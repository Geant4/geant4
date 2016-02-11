@echo off

rem -- Common setting 
if DEFINED VERSION_BENCHMARKS  echo "Using benchmarks tag: %VERSION_BENCHMARKS%"

if NOT DEFINED MODE    set MODE=nightly
if NOT DEFINED WORKDIR set WORKDIR=D:/build/cdash
if NOT DEFINED VERSION set VERSION=g4tags-dev
echo "Mode=%MODE%"

set BUILDDIR=%WORKDIR%/%MODE%/%CONFIG%

rem -- Jenkins uses workspace
if DEFINED WORKSPACE   set WORKDIR=%cd%
if DEFINED WORKSPACE   set BUILDDIR=%WORKDIR:\=/%
if DEFINED WORKSPACE   echo "Using jenkings workspace .%WORKSPACE%. and workdir: .%WORKDIR%."
)

if NOT EXIST %BUILDDIR%\ (
   echo "Error: No Directory to build in, using BUILDDIR=%BUILDDIR%"
)


set SOURCE=%BUILDDIR%/%VERSION%
set BINARY=%BUILDDIR%/build

rem ---default build option, to be modified for non default buildtypes
set G4_XOPTS=-DGEANT4_USE_OPENGL_WIN32=ON

rem -- only need to check for options from jenkins
IF NOT DEFINED BUILDOPTIONS goto NOBUILDOPT

rem - turn on delayed exansion to delay expansion of variable in IF
SETLOCAL ENABLEDELAYEDEXPANSION

rem -- for loop requires setlocal, but 
rem -- in setlocal changes are local, so will be lost after endlocal, or end of .bat
rem --  attempt tp set G4_XOPTS global, will only keep last expansion doe within loop,
rem --   i.e. when two options modify G4_XOPTS, first is lost.
rem --  workaround: set a (global) variable for each option, and concatenate after endlocal

rem --- echo Buildoptions are: %BUILDOPTIONS%

FOR %%a IN (%BUILDOPTIONS%) DO ( 
	rem --- echo option %%a with opts !G4_XOPTS!
	IF        "%%a"=="staticlibs" (
	     ENDLOCAL
		set OPT_STATIC=-DBUILD_STATIC_LIBS=ON;-DBUILD_SHARED_LIBS=OFF
		  SETLOCAL ENABLEDELAYEDEXPANSION
	) ELSE IF "%%a"=="EpCheck" (
	     ENDLOCAL
		set G4Hadronic_epReportLevel=-3
		  SETLOCAL ENABLEDELAYEDEXPANSION
	) ELSE IF "%%a"=="BoundsCheck" (
		echo Error: %%a not -yet- supported on windows
	) ELSE IF "%%a"=="UseUsolids" (
	     ENDLOCAL
		set OPT_USOLID=-DGEANT4_USE_USOLIDS=ON
		  SETLOCAL ENABLEDELAYEDEXPANSION
	) ELSE IF "%%a"=="UseInternalCLHEP" (
		echo Error: %%a is still default on windows
	) ELSE IF "%%a"=="UseGranularCLHEP" (
		echo Error: %%a not -yet- supported on windows 
	) ELSE IF 1==1 (
		echo Error: unknown option passed in BUILDOPTIONS: %%a
		REM --exit 1
	)
)

	
ENDLOCAL

:NOBUILDOPT


IF "%THREAD%" == "MT" (
   set OPT_MT=-DGEANT4_BUILD_MULTITHREADED=ON
) ELSE IF "%THREAD" == "MTmax" (
   set OPT_MT=-DGEANT4_BUILD_MULTITHREADED=ON
	set G4FORCENUMBEROFTHREADS="max"
)	


IF DEFINED OPT_STATIC  set G4_XOPTS=%G4_XOPTS%;%OPT_STATIC%
IF DEFINED OPT_USOLID  set G4_XOPTS=%G4_XOPTS%;%OPT_USOLID%
IF DEFINED OPT_MT      set G4_XOPTS=%G4_XOPTS%;%OPT_MT%

set Path=%Path%;C:\Python27

rem - python --version 

rem -  show environment:
rem - SET

echo Geant4 CMake options: %G4_XOPTS%
