import json
import pandas as pd
import matplotlib.pyplot as plt
import os

# 确定脚本所在目录
script_dir = os.path.dirname(__file__)

# 加载 CSV 数据
csv_file = os.path.join(script_dir, 'adams_output.csv')  # 根据需要调整 CSV 文件名
data = pd.read_csv(csv_file)

# 从 JSON 中加载要绘制的列
json_file = os.path.join(script_dir, 'plot_columns.json')  # 根据需要调整 JSON 文件名
with open(json_file, 'r') as file:
    plot_info = json.load(file)
subplots_info = plot_info['subplots']

# 为每个子图创建绘图
num_subplots = len(subplots_info)
fig, axs = plt.subplots(num_subplots, 1, figsize=(8, 4 * num_subplots))

# 确保 axs 是列表
if num_subplots == 1:
    axs = [axs]

for i, (subplot_name, y_cols) in enumerate(subplots_info.items()):
    for y_col in y_cols:
        axs[i].plot(data['t'], data[y_col], label=f'{y_col} vs t')
    axs[i].set_xlabel('Time (t)')
    axs[i].set_ylabel('Values')
    axs[i].legend()
    axs[i].set_title(subplot_name)

plt.tight_layout()
plt.show()
