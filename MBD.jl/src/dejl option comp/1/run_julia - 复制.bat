@echo off
REM %1 is the JSON string passed by the Python script

REM Call Julia script and pass the JSON string as an argument
julia.exe --project=%JULIAENV% "D:\OneDrive\Articles\10.Working\[D21][20211009]ContactMechanics\MBD.jl\MBD.jl\src\dejl option comp\main agent call.jl"
