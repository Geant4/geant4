@echo off
set this=%~d0%~p0
set opts=%1
:next
shift
if '%1'=='' goto end 
set opts=%opts% %1
goto next
:end
python %this%g4tagsvn.py %opts%

