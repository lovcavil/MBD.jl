import shutil

# Source file path
source_file_path = r"C:\OneDrive\Articles\10.Working\[D21][20211009]ContactMechanics\31.AdamsCompare\sph_pla.tab"

# Destination file path
destination_file_path = r"C:\OneDrive\Articles\10.Working\[D21][20211009]ContactMechanics\MBD.jl\demo205\sph_pla.tab"

# Copy the file
shutil.copyfile(source_file_path, destination_file_path)

print(f"File copied successfully from {source_file_path} to {destination_file_path}")
