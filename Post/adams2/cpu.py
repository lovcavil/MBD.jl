import pandas as pd
import matplotlib.pyplot as plt

# Data from the table
data = {
    'co': [2.71484375, 2.87890625, 2.40625, 2.44140625, 2.48046875, 2.66015625, 2.62890625],
    'linear': [2.00390625, 2.75, 2.37890625, 2.06640625, 2.12890625, 2.28125, 2.43359375],
    'smooth': [2.24609375, 2.23828125, 2.71875, 2.125, 2.2890625, 2.14453125, 2.3984375],
    'modified': [2.42578125, 2.59765625, 4.390625, 3.12890625, 2.0546875, 2.58203125, 3.41015625]
}

df = pd.DataFrame(data, index=['F', 'Hc', 'Ln', 'Hm', 'Go', 'Hg', 'Gh'])

# Create Clustered Bar Chart
plt.figure(figsize=(12, 6))
df.plot(kind='bar', width=0.8, figsize=(12, 6))
plt.title('Clustered Bar Chart of Values')
plt.xlabel('Factor r (7 levels)')
plt.ylabel('Values')
plt.xticks(rotation=0)
plt.legend(title='Factor c (4 levels)')
plt.grid(axis='y')

plt.tight_layout()
plt.show()
