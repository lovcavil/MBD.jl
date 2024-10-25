# ThreadPoolExecutor this is phase ok version
import json
import subprocess
import os
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
from prefect import task, flow, get_run_logger
# Function to run each task
@task(name="run_task", tags=[f"test"])
def run_task(task_id, key, value, folder, bat_file, jl_file):
    # Print before starting the task
    print(f"Task {task_id} started.")
    
    # Create a separate JSON object for this task
    single_task = {key: value}

    # Convert the task to a JSON string
    json_task_string = json.dumps(single_task)

    try:
        # Start the process
        process = subprocess.Popen(
            [bat_file, jl_file, folder, str(task_id), json_task_string],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, encoding="utf-8",
            creationflags=subprocess.CREATE_NEW_PROCESS_GROUP  # For Windows
        )

        # Wait for the process to complete
        stdout, stderr = process.communicate()

        # Print the captured stdout and stderr
        #print(stdout)  # Print standard output from the subprocess
        #print(stderr)  # Print any errors that occurred

        # Check if the process returned any errors
        if process.returncode != 0:
            print(f"Error executing Julia for task {task_id}: {stderr}")
            return False  # Mark the task as failed
        else:
            # Ensure that the output JSON file exists
            output_json_file = os.path.join(os.path.dirname(__file__), folder, f"{str(task_id)}.json")
            if os.path.exists(output_json_file):
                with open(output_json_file, "r") as output_file:
                    julia_output = json.load(output_file)
                    print(f"Task {task_id} - Times:", julia_output.get("times"))
                print(f"Task {task_id} completed successfully.")
                return True  # Task completed successfully
            else:
                print(f"Task {task_id} output file not found.")
                print(f"Task {task_id} failed.")
                return False  # Task failed because output was not generated

    except Exception as e:
        print(f"Task {task_id} encountered an error: {e}")
        return False  # Mark the task as failed due to an exception
    
@flow(name="dejlcomp"
     )
def main():
    # Read JSON file
    file_path = os.path.join(os.path.dirname(__file__), "my_data.json")
    
    # Load JSON data into Python dictionary
    with open(file_path, "r") as json_file:
        json_content = json.load(json_file)

    # Prepare folder name based on the current time
    now = datetime.now()
    folder = now.strftime("%Y-%m-%d-%H-%M-%S")
    os.makedirs(os.path.join(os.path.dirname(__file__), folder), exist_ok=True)  # Create the folder if it doesn't exist

    # Call the batch file
    bat_file = os.path.join(os.path.dirname(__file__), "20.run_julia.bat")
    jl_file = os.path.join(os.path.dirname(__file__), "30.main agent call.jl")

    # Use ThreadPoolExecutor to handle multiple tasks concurrently
    num_threads = 2  # Number of concurrent threads
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        # Submit tasks to the pool
        future_to_task = {
            executor.submit(run_task, i, key, value, folder, bat_file, jl_file): (i, key)
            for i, (key, value) in enumerate(json_content.items(), start=1)
        }

        # Process the results as they are completed
        for future in as_completed(future_to_task):
            task_id, key = future_to_task[future]
            try:
                result = future.result()
                if result:
                    print(f"Task {task_id} ended successfully.")
                else:
                    print(f"Task {task_id} ended with failure.")
            except Exception as e:
                print(f"Task {task_id} generated an exception: {e}")

if __name__ == "__main__":
    main()