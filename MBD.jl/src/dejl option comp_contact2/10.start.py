import json
import subprocess
import os
from datetime import datetime
# def main():
#     # Read JSON file
#     file_path = os.path.join(os.path.dirname(__file__), "my_data.json")
    
#     # Load JSON data into Python dictionary
#     with open(file_path, "r") as json_file:
#         json_content = json.load(json_file)

#     # Convert the dictionary to a JSON string
#     json_string = json.dumps(json_content)

#     # Call the batch file and pass the JSON string as an argument
#     bat_file = os.path.join(os.path.dirname(__file__), "run_julia.bat")

#     # Note: You may need to handle quotes or escape characters depending on the complexity of the JSON.
#     subprocess.call([bat_file,"D:/OneDrive/Articles/10.Working/[D21][20211009]ContactMechanics/MBD.jl/MBD.jl/src/dejl option comp/main agent call.jl",json_string])

def main():
    # Read JSON file
    file_path = os.path.join(os.path.dirname(__file__), "my_data.json")
    
    # Load JSON data into Python dictionary
    with open(file_path, "r") as json_file:
        json_content = json.load(json_file)

    # Convert the dictionary to a JSON string
    json_string = json.dumps(json_content)
    # Call the batch file and pass the JSON string as an argument
    bat_file = os.path.join(os.path.dirname(__file__), "20.run_julia.bat")
    jl_file=os.path.join(os.path.dirname(__file__), "30.main agent call.jl")
    now = datetime.now()
    folder=now.strftime("%Y-%m-%d-%H-%M-%S")
    num=1
    print(folder)
    result = subprocess.run([bat_file,jl_file,folder,str(num), json_string],
                             capture_output=True, text=True, encoding="utf-8")

    # Check if Julia returned any errors
    if result.returncode != 0:
        print(f"Error executing Julia: {result.stderr}")
    else:
        # Read the result from the output JSON file (produced by Julia)
        output_json_file = os.path.join(os.path.dirname(__file__),folder, f"{str(num)}.json")  # Assuming the batch prints the output file path
        with open(output_json_file, "r") as output_file:
            julia_output = json.load(output_file)
            #print("Energies:", julia_output["energies"])
            print("Times:", julia_output["times"])

if __name__ == "__main__":
    main()
