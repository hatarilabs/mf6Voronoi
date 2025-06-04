Copy /Y /A "HillslopeModel_v1.mfsim.nam" "mfsim.nam"
call "C:\Program Files\USGS\ModelMuse5\bin\ModelMonitor.exe" -m C:\WRDAPP\mf6.6.0_win64\bin\mf6.exe -n mfsim.nam -mv 6
Start Notepad.exe HillslopeModel_v1.lst
Start Notepad.exe mfsim.lst
pause
