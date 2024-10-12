import argparse
from prefect import task, flow, get_run_logger
import subprocess
from datetime import datetime
import os

@task
def run_bat_file(command):
    logger = get_run_logger()
    logger.info("Starting batch file execution")

    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        shell=True,
        encoding='utf-8'
    )

    # Read the output line by line
    for output in process.stdout:
        if output.strip():
            #logger.info(output.strip())
            print(output.strip())

    # Capture any remaining output
    stdout, stderr = process.communicate()

    logger.info("Remaining STDOUT:")
    if stdout:
        #logger.info(stdout)
        print(stdout)
    
    if stderr:
        logger.error("STDERR:")
        logger.error(stderr)
        print(stderr)

    # Check the return code for any errors
    if process.returncode != 0:
        logger.error(f"Batch file exited with an error. Return code: {process.returncode}")
        print(f"Batch file exited with an error. Return code: {process.returncode}")
        print("STDERR:")
        print(stderr)
        return process.returncode

    return 0

@flow(flow_run_name="{name}-app{app}-te{t}-dt{dt}-{cfg2}-{ct}")
def CAKD_SIM(dt: str, cfg1: str, cfg2: str, app: str, t: str,name: str,basefolder: str,ct,
             rt: str, at: str, ma: str, mi,sn):
    base_path = basefolder
    bat_file = f'{base_path}prefect_run.bat'
    script_file = f'{base_path}MBD.jl/src/door341_350_entrysolvercomp_fix.jl'
    command = [bat_file, script_file, dt, cfg1, cfg2, app, t, name, base_path, rt, at, ma, mi,sn]
    
    print(f"Script: {script_file} => Run: {name}")
    print(f"dt: {dt} | Config: {cfg1} | Contact: {cfg2} | t: {t}")
    
    result = run_bat_file(command)
    return result

if __name__ == "__main__":
    onedrivefolder = os.getenv('OneDrive', 'F:/OneDrive')

    parser = argparse.ArgumentParser(description="Run CAKD Simulation")
    parser.add_argument('--name', type=str, default="hello_run_prefect", help="Name Run")
    parser.add_argument('--dt', type=str, default="1e-2", help="Time step value (default: 1e-2)")
    parser.add_argument('--cfg1', type=str, default="aa", help="Configuration 1 (default: aa)")
    parser.add_argument('--cfg2', type=str, default="contact_off", help="Configuration 2 (default: contact_off)")
    parser.add_argument('--app', type=str, default="341", help="Application identifier (default: 341)")
    parser.add_argument('--t', type=str, default="0.5", help="Time step value (default: 1.0)")
    parser.add_argument('--bf', type=str,
                        default=f"{onedrivefolder}/Articles/10.Working/[D21][20211009]ContactMechanics/MBD.jl/",
                        help="")
    parser.add_argument('--rt', type=str, default="1e-4", help="")
    parser.add_argument('--at', type=str, default="1e-4", help="")
    parser.add_argument('--ma', type=str, default="1e-2", help="Time step value (default: 1e-2)")
    parser.add_argument('--mi', type=str, default="1e-4", help="Time step value (default: 1e-2)")
    parser.add_argument('--sn', type=str, default="Euler", help="")
    args = parser.parse_args()
    ct= datetime.now().strftime("%Y/%m/%d-%H:%M:%S")
    CAKD_SIM(dt=args.dt, cfg1=args.cfg1, cfg2=args.cfg2, app=args.app,t=args.t,name=args.name,
             basefolder=args.bf,ct=ct,rt=args.rt,at=args.at,ma=args.ma,mi=args.mi,sn=args.sn)
