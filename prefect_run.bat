REM Set environment variables
REM  tb14  --project=C:\Users\lovca\.julia\environments\diffeqpy


set filename=%1
set arg2=%2
set arg3=%3
set arg4=%4
set arg5=%5
set time=%6
set runname=%7
set bf=%8
set rt=%9
for /L %%i in (1,1,9) do shift
set at=%1
set ma=%2
set mi=%3
set sn=%4


julia.exe %JULIAENV% %filename% %arg2% %arg3% %arg4% %arg5% %time% %runname% %bf% %rt% %at% %ma% %mi% %sn%
exit