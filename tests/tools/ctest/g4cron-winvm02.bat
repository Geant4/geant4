@echo off  
svn update %~d0%~p0
cmd /C "%~d0%~p0g4cron-winvm02-vc9.bat"
cmd /C "%~d0%~p0g4cron-winvm02-vc10.bat"
