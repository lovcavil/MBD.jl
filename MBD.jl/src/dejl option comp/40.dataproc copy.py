import h5py
import pandas as pd

# Open the file
fn="D:\OneDrive\Articles/10.Working\[D21][20211009]ContactMechanics\MBD.jl\MBD.jl\src\dejl option comp/2024-10-24-20-42-43/3.csv.jld2"




def load_results_from_jld(file_path):
    # Open the file using h5py
    with h5py.File(file_path, "r") as f:
        data = {}
        
        # Iterate through all the items in the file (could be groups or datasets)
        def visit_fn(name, obj):
            if isinstance(obj, h5py.Dataset):
                # If it's a scalar dataset
                if obj.shape == ():
                    data[name] = obj[()]  # Store the scalar value
                # If it's an array or table, convert to a pandas DataFrame
                else:
                    data[name] = pd.DataFrame(obj[:])
            elif isinstance(obj, h5py.Group):
                print(f"Skipping group: {name}")

        # Recursively visit all items in the file
        f.visititems(visit_fn)

    return data


loaded_data = load_results_from_jld(fn)

# Access and print loaded data
for key, value in loaded_data.items():
    if isinstance(value, pd.DataFrame):
        print(f"DataFrame {key}:")
        print(value)
    else:
        print(f"Scalar {key}: {value}")
