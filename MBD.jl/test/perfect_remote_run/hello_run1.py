import subprocess

# Path to your BAT file
base = 'C:/OneDrive/Articles/10.Working/[D21][20211009]ContactMechanics/MBD.jl/'
bat_file = base + 'prefect_run.bat'

# Arguments to pass to the BAT file
file = base + 'MBD.jl/src/door341_350_entrysolvercomp_fix.jl'
arg2 = '1e-2'
arg3 = 'aaa'
arg4 = 'contact2'
arg5 = '350'
command = [bat_file, file, arg2, arg3, arg4, arg5]

# Run the command and capture output in real-time
process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True)

# Read the output line by line
while True:
    output = process.stdout.readline()
    if output == '' and process.poll() is not None:
        break
    if output:
        print(output.strip())

# Capture any remaining output
stdout, stderr = process.communicate()

print("Remaining STDOUT:")
print(stdout)
print("STDERR:")
print(stderr)

# Check the return code for any errors
if process.returncode != 0:
    print(f"Batch file exited with an error. Return code: {process.returncode}")
    print("STDERR:")
    print(stderr)
