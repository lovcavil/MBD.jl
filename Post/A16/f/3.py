# List of filenames
filenames = [
    '213.asc', '214.asc', '215.asc', '216.asc', '217.asc', 
    '218.asc', '219.asc', '313.asc', '314.asc', '315.asc', 
    '316.asc', '317.asc', '318.asc', '319.asc', '413.asc', 
    '414.asc', '415.asc', '416.asc', '417.asc', '418.asc', 
    '419.asc', '513.asc', '514.asc', '515.asc', '516.asc', 
    '517.asc', '518.asc', '519.asc'
]

# Group by first character
grouped = {}
for filename in filenames:
    first_char = filename[0]  # Get the first character
    if first_char not in grouped:
        grouped[first_char] = []  # Create a new list if it doesn't exist
    grouped[first_char].append(filename)  # Append the filename to the appropriate group

# Print the grouped result
for key, group in grouped.items():
    print(f"{key}: {group}")