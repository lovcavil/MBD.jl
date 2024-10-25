@echo off
REM %1 is the JSON string passed by the Python script
set filename=%1
set folder=%2
set subname=%3
set json=%4
REM Call Julia script and pass the JSON string as an argument
julia.exe --project=%JULIAENV% %filename% %folder% %subname% %json%
