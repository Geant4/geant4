@ECHO off
SETLOCAL

REM //////////////////////////////////////////////////////////////
REM /// C++ //////////////////////////////////////////////////////
REM //////////////////////////////////////////////////////////////
CL.exe 1> NUL 2> NUL
IF ERRORLEVEL 1 (
  ECHO CL.exe program not found. You have to setup VisualC++ environment.
  GOTO build_return
)

SET flags=/nologo /DWIN32 /MD /O2 /W3  /GX /GR 

SET head=..\..\..

SET cppflags=%flags%
SET cppflags=%cppflags% /I%head%\tools /I%head%\tools

REM //////////////////////////////////////////////////////////////
REM /// F77 //////////////////////////////////////////////////////
REM //////////////////////////////////////////////////////////////
F77.exe 1> NUL 2> NUL
IF ERRORLEVEL 1 (
  ECHO F77.exe program not found. You have to setup your fortran.
  GOTO build_return
)

SET f77flags=/nologo /MD

SET CERNLIB_home=C:\cern\pro

IF NOT EXIST %CERNLIB_home% (
  ECHO %CERNLIB_home% not found.
  GOTO build_return
)

REM //////////////////////////////////////////////////////////////
REM /// hello_f77.f //////////////////////////////////////////////
REM //////////////////////////////////////////////////////////////
F77.exe /nologo /compile_only /object:.\hello_f77.obj .\hello_f77.f
LINK.exe /nologo /entry:mainCRTStartup /out:.\tools_test_hello_f77.exe .\hello_f77.obj
DEL hello_f77.obj

REM //////////////////////////////////////////////////////////////
REM /// hbook.f //////////////////////////////////////////////////
REM //////////////////////////////////////////////////////////////
F77.exe /nologo /compile_only /object:.\hbook.obj .\hbook.f
LINK.exe /nologo /entry:mainCRTStartup /out:.\tools_test_hbook.exe .\hbook.obj /libpath:%CERNLIB_home%\lib packlib.lib advapi32.lib ws2_32.lib 
DEL hbook.obj

REM //////////////////////////////////////////////////////////////
REM /// chbook.cpp ///////////////////////////////////////////////
REM //////////////////////////////////////////////////////////////
F77.exe %f77flags% /compile_only /object:.\close.obj ..\..\tools\hbook\close.f
F77.exe %f77flags% /compile_only /object:.\setpawc.obj ..\..\tools\hbook\setpawc.f

CL.exe %cppflags% /c /Fo.\chbook.obj /Tp.\chbook.cpp

LINK.exe /nologo /OPT:NOREF /out:.\tools_test_chbook.exe .\chbook.obj .\setpawc.obj .\close.obj /libpath:%CERNLIB_home%\lib packmd.lib shiftmd.lib advapi32.lib ws2_32.lib

DEL chbook.obj

REM //////////////////////////////////////////////////////////////
REM /// whbook.cpp ///////////////////////////////////////////////
REM //////////////////////////////////////////////////////////////

CL.exe %cppflags% /c /Fo.\whbook.obj /Tp.\whbook.cpp
F77.exe %f77flags% /compile_only /object:.\setntuc.obj ..\..\tools\hbook\setntuc.f

LINK.exe /nologo /OPT:NOREF /out:.\tools_test_whbook.exe .\whbook.obj .\setntuc.obj .\setpawc.obj .\close.obj /libpath:%CERNLIB_home%\lib packmd.lib shiftmd.lib advapi32.lib ws2_32.lib

DEL whbook.obj
DEL setntuc.obj setpawc.obj close.obj 

REM //////////////////////////////////////////////////////////////
REM //////////////////////////////////////////////////////////////
REM //////////////////////////////////////////////////////////////

:build_return
ENDLOCAL
@ECHO on
