# Function to process the file
def process_adams_rsp_file(file_path):
    structured_data = {}
    current_section = None

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith("FILE TYPE AND SOURCE"):
                current_section = "file_info"
                structured_data[current_section] = {}
            elif line.startswith("DATA SET TITLE") or line.startswith("RESULTS TITLE"):
                current_section = line.lower().split()[0] + "_" + line.lower().split()[1]
                structured_data[current_section] = {}
            elif line.startswith("!"):
                current_section = None
            elif current_section:
                # Split line into parts and process according to the current section
                parts = line.split()
                if len(parts) >= 2:  # Basic check to ignore empty lines or lines that don't match expected pattern
                    # Assuming the last element is the value and the rest form the key
                    key = " ".join(parts[:-1])
                    value = parts[-1]
                    structured_data[current_section][key] = value

    return structured_data

# Replace 'file_path' with the actual path to your file
file_path = 'Post\Adams.res'
structured_data = process_adams_rsp_file(file_path)

# Example to print structured data
for section, data in structured_data.items():
    print(f"Section: {section}")
    for key, value in data.items():
        print(f"  {key}: {value}")
