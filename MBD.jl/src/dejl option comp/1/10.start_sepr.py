import json
import subprocess
import os
from datetime import datetime

def main():
    # Read JSON file
    file_path = os.path.join(os.path.dirname(__file__), "my_data.json")
    
    # Load JSON data into Python dictionary
    with open(file_path, "r") as json_file:
        json_content = json.load(json_file)

    # Prepare folder name based on the current time
    now = datetime.now()
    folder = now.strftime("%Y-%m-%d-%H-%M-%S")
    os.makedirs(os.path.join(os.path.dirname(__file__),folder), exist_ok=True)  # Create the folder if it doesn't exist

    # Call the batch file
    bat_file = os.path.join(os.path.dirname(__file__), "run_julia.bat")
    jl_file = os.path.join(os.path.dirname(__file__), "main agent call.jl")

    # Dispatch each key-value pair as a separate task
    for i, (key, value) in enumerate(json_content.items(), start=1):
        # Create a separate JSON object for this task
        single_task = {key: value}

        # Convert the task to a JSON string
        json_task_string = json.dumps(single_task)

        # Call the batch file and pass the JSON string directly to Julia
        result = subprocess.run([bat_file, jl_file, folder, str(i), json_task_string], capture_output=True, text=True, encoding="utf-8")

        # Print the captured stdout and stderr
        print(result.stdout)  # Print standard output from the subprocess
        print(result.stderr)  # Print any errors that occurred

        # Check if Julia returned any errors
        if result.returncode != 0:
            print(f"Error executing Julia for task {i}: {result.stderr}")
        else:
            # Read the result from the output JSON file (produced by Julia)
            output_json_file = os.path.join(os.path.dirname(__file__),folder, f"{str(i)}.json")
            if os.path.exists(output_json_file):
                with open(output_json_file, "r") as output_file:
                    julia_output = json.load(output_file)
                    #print(f"Task {i} - Energies:", julia_output.get("energies"))
                    print(f"Task {i} - Times:", julia_output.get("times"))

if __name__ == "__main__":
    main()
