import subprocess

# Path to your BAT file
# Arguments to pass to the BAT file
base='C:/OneDrive/Articles/10.Working/[D21][20211009]ContactMechanics/MBD.jl/'
bat_file = base+'prefect_run.bat'

file=base+'MBD.jl/src/door341_350_entrysolvercomp_fix.jl'

arg2='1e-2'
arg3="aaa"
arg4="contact2"
arg5="350"
command = [bat_file,file,arg2,arg3,arg4,arg5]

# Run the command
result = subprocess.run(command, capture_output=True, text=True, shell=True)

# Print the output and errors (if any)
print("STDOUT:")
print(result.stdout)
print("STDERR:")
print(result.stderr)
