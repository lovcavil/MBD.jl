import json
import subprocess
import os
from datetime import datetime
from prefect import task, flow, get_run_logger
from prefect.states import Completed, Failed
from prefect.task_engine import TaskRunTimeoutError
from prefect.client.orchestration import get_client
import asyncio
import time

# Function to set concurrency limit
async def set_concurrency_limit():
    async with get_client() as client:
        # Set a concurrency limit of 28 on the 'test' tag
        limit_id = await client.create_concurrency_limit(
            tag="test",
            concurrency_limit=18  # Adjust limit as needed
        )
        print(f"Concurrency limit set for 'test': {limit_id}")
async def check_concurrency_limit():
    async with get_client() as client:
        concurrency_limit = await client.read_concurrency_limit_by_tag("test")
        
        logger = get_run_logger()
        active_slots_count = len(concurrency_limit.active_slots)
        limit = concurrency_limit.concurrency_limit
        
        logger.info(f"Active slots: {active_slots_count}/{limit}")

        # Return True if there is space to submit a new task, otherwise wait a short period
        if active_slots_count < limit:
            return True
        else:
            time.sleep(1)  # Wait briefly before the next check to avoid overloading checks
            return False


# Function to run each task
@task(name="run_task", tags=["test"])
def run_task(task_id, key, value, folder, bat_file, jl_file):
    logger = get_run_logger()
    logger.info(f"Task {task_id} started.")

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
        print(stdout)
        # Check if the process returned any errors
        if process.returncode != 0:
            logger.error(f"Error executing Julia for task {task_id}: {stderr}")
            return Failed(message="Error executing Julia for task")
        else:
            # Ensure that the output JSON file exists
            output_state = os.path.join(os.path.dirname(__file__), folder, f"{str(task_id)}.csv.state")
            if os.path.exists(output_state):
                with open(output_state, 'r') as file:
                    # Read the content of the file
                    file_content = file.read()
                    if file_content=="timeout":
                        raise TaskRunTimeoutError
            output_json_file = os.path.join(os.path.dirname(__file__), folder, f"{str(task_id)}.json")
            if os.path.exists(output_json_file):
                with open(output_json_file, "r") as output_file:
                    julia_output = json.load(output_file)
                    logger.info(f"Task {task_id} - Times: {julia_output.get('times')}")
                logger.info(f"Task {task_id} completed successfully.")
                return True  # Task completed successfully
            else:
                logger.error(f"Task {task_id} output file not found.")
                return Failed(message=f"Task {task_id} output file not found.")  # Task failed because output was not generated

    except Exception as e:
        logger.error(f"Task {task_id} encountered an error: {e}")
        return False  # Mark the task as failed due to an exception

@flow(name="dejlcomp", flow_run_name="{folder}")
async def main_flow(folder):
    await set_concurrency_limit()
    logger = get_run_logger()

    # Read JSON file
    file_path = os.path.join(os.path.dirname(__file__), "my_data.json")
    with open(file_path, "r") as json_file:
        json_content = json.load(json_file)

    # Prepare folder name based on the current time
    os.makedirs(os.path.join(os.path.dirname(__file__), folder), exist_ok=True)  # Create the folder if it doesn't exist

    # Call the batch file and Julia script
    bat_file = os.path.join(os.path.dirname(__file__), "20.run_julia.bat")
    jl_file = os.path.join(os.path.dirname(__file__), "30.main agent call.jl")
    
    
    # Submit tasks one by one based on concurrency availability
    results = []
    for i, (key, value) in enumerate(json_content.items(), start=70):
        # Check concurrency limit availability before submitting the task
        while True:
            concurrency_available = await check_concurrency_limit()
            if concurrency_available:
                logger.info(f"Concurrency available, submitting task {i}")
                logger.info(f"Progress: {i}/{len(json_content)}")
                result = run_task.submit(i, key, value, folder, bat_file, jl_file)
                results.append(result)
                time.sleep(5)
                break  # Task submitted, move to the next one
            else:
                logger.info(f"No concurrency slot available for task {i}, waiting...")
                time.sleep(5)  # Wait for 5 seconds before checking again

    # Collect all task results
    for result in results:
        if result.result():
            logger.info("A task ended successfully.")
        else:
            logger.error("A task ended with failure.")


if __name__ == "__main__":
    now = datetime.now()
    folder = now.strftime("%Y-%m-%d-%H-%M-%S")
    asyncio.run(main_flow(folder))
