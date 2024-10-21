import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import gaussian_kde
from scipy.interpolate import interp1d  # Added import

def create_custom_colormap():
    """创建自定义的颜色映射"""
    rgb_colors = [
        (38, 70, 83),
        (40, 114, 113),
        (41, 157, 143),
        (138, 176, 125),
        (232, 197, 107),
        (243, 162, 97),
        (230, 111, 81)
    ]
    normalized_colors = [(r / 255, g / 255, b / 255) for r, g, b in rgb_colors]
    cmap_name = 'custom_cmap'
    return LinearSegmentedColormap.from_list(cmap_name, normalized_colors)

def load_test(file):
    """加载测试数据"""
    column_names = [
        "StrainB1_1", "StrainB2_1", "StrainB2_2", "StrainB2_3", "StrainB2_4", "StrainB3_1", "StrainB3_2",
        "StrainB3_3", "StrainB3_4", "accA1_X", "accA1_Y", "accA1_Z", "accA2_X", "accA2_Y", "accA2_Z",
        "accA3_X", "accA3_Y", "accA3_Z", "accA4_X", "accA4_Y", "accA4_Z", "accA5_X", "accA5_Y", "accA5_Z",
        "accA6_X", "accA6_Y", "accA6_Z", "accA7_X", "accA7_Y", "accA7_Z", "canpress"
    ]
    column_types = {name: 'float64' for name in column_names}
    try:
        df = pd.read_csv(
            file,
            delim_whitespace=True,
            skiprows=26,
            names=column_names,
            dtype=column_types,
            engine='python'
        )
        return df
    except Exception as e:
        print(f"加载文件 {file} 时出错: {e}")
        return None

def post_test(file):
    """处理测试数据，添加时间列"""
    test = load_test(file)
    if test is not None:
        delta = 1.9531250e-3
        t = np.arange(delta, delta * (len(test) + 1), delta)
        test['t'] = t
    return test

def normalize_label_data(label_data):
    """归一化标签数据"""
    min_val = min(label_data.values())
    max_val = max(label_data.values())
    normalized_label_data = {
        key: (val - min_val) / (max_val - min_val) for key, val in label_data.items()
    }
    return normalized_label_data, min_val, max_val

def plot_smoothed_pressure(label_data, normalized_label_data, cmap, min_val, max_val):
    """绘制平滑的压力曲线"""
    dforcesme = {}
    dforcesme0 = []
    plt.figure(figsize=(10,6.5),facecolor='white')
    ax = plt.subplot(1, 1, 1)
    ax.grid(True)

    for i1 in range(24, 31):
        for i2 in range(1, 9):
            if (i1, i2) in [(24, 6), (31, 1), (32, 4)]:
                continue
            file_path = f"D:/D11.Proj/A009.BYD_R2/20220710/DATA_PROC/CUT/CLOSE_NORMAL/MR_door_{i1}{i2}.asc"
            if (i1, i2) in normalized_label_data and os.path.isfile(file_path):
                norm_val = normalized_label_data[(i1, i2)]
                color = cmap(norm_val)
                test = post_test(file_path)
                if test is None:
                    continue

                x = test['t'].to_numpy()
                y = test['canpress'].to_numpy()
                aaa = y > 10
                V1 = np.argmax(aaa) if np.any(aaa) else None
                M = np.max(y)
                V2 = np.argmax(y)
                V3 = len(y) - np.argmax(aaa[::-1]) - 1 if np.any(aaa) else None

                if V1 is not None and V3 is not None:
                    dforcesme[(i1, i2)] = [x[V1], x[V2], x[V3], M]
                    dforcesme0.append([x[V1], x[V2], x[V3], M])

                    # 平滑数据
                    window_length = min(17, len(y) if len(y) % 2 != 0 else len(y) - 1)
                    window_length = max(window_length, 5)
                    if window_length % 2 == 0:
                        window_length += 1

                    y_smoothed = savgol_filter(y, window_length=window_length, polyorder=3)

                    # 绘制曲线和填充区域
                    st = 1
                    en = min(V3 + 40, len(x))
                    ax.plot(x[st:en], y_smoothed[st:en], color=color, alpha=0.5)
                    ax.fill_between(x[st:en], y_smoothed[st:en], color=color, alpha=0.3)

    # 添加颜色条，传入 min_val 和 max_val
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=min_val, vmax=max_val))
    sm.set_array([])
    # plt.colorbar(sm, ax=ax, label='Velocity')
    textstyle = {'family': 'Times New Roman', 'size': 28}
    tick_fontsize = 20
    tick_font = 'Times New Roman'
    ax.set_xlabel('t (s)', fontdict=textstyle)
    ax.set_ylabel('Driving Forces (N)', fontdict=textstyle)
    # ax.set_title('Smoothed Pressure Over Time with Velocity Color Mapping')
    plt.yticks(fontsize=tick_fontsize, fontfamily=tick_font)
    plt.xticks(fontsize=tick_fontsize, fontfamily=tick_font)
    #plt.show()
    output_folder = r"D:\OneDrive\Articles\10.Working\[D21][20211009]ContactMechanics\MBD.jl\Post\A16"
    output_path = os.path.join(output_folder, 'drives2.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    return dforcesme, dforcesme0

def plot_violin(label_data, cmap):
    """绘制带有平滑颜色映射的小提琴图，并添加数据点"""
    data = list(label_data.values())
    kde = gaussian_kde(data)
    y = np.linspace(min(data), max(data)+0.1, 500)
    density = kde(y)
    density = density / density.max()

    x = np.linspace(-1, 1, 200)
    X, Y = np.meshgrid(x, y)
    violin_mask = np.abs(X) <= density[:, np.newaxis]

    norm = plt.Normalize(vmin=min(y), vmax=max(y)+0.1)
    colors = cmap(norm(Y))
    colors[..., -1] = violin_mask.astype(float)

    fig, ax = plt.subplots(figsize=(2, 8))
    ax.imshow(
        colors,
        extent=[x.min(), x.max(), y.min(), y.max()],
        aspect='auto',
        origin='lower'
    )
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.set_ylabel('Value')

    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(y.min(), y.max())

    # 添加数据点
    density_interp = interp1d(y, density, kind='linear', bounds_error=False, fill_value=0)
    xi_list = []
    yi_list = []
    for yi in data:
        hi = density_interp(yi)
        xi = 0
        xi_list.append(xi)
        yi_list.append(yi+np.random.uniform(-0.01,0.01))
    ax.scatter(xi_list, yi_list, color='black', s=10, zorder=3)

    #plt.show()
    output_folder = r"D:\OneDrive\Articles\10.Working\[D21][20211009]ContactMechanics\MBD.jl\Post\A16"
    output_path = os.path.join(output_folder, 'drives.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')

def main():
    """主函数"""
    # 创建自定义颜色映射
    cmap = create_custom_colormap()





    # 标签数据
    label_data = {
        (24, 1): 0.86, (24, 2): 0.80, (24, 3): 0.82, (24, 4): 1.00, (24, 5): 1.04, (24, 6): 1.20, (24, 7): 1.21,
        (25, 1): 1.21, (25, 2): 1.17, (25, 3): 1.06, (25, 4): 0.94, (25, 5): 0.89, (25, 6): 0.68, (25, 7): 0.95, (25, 8): 0.73,
        (26, 1): 0.80, (26, 2): 0.87, (26, 3): 0.87, (26, 4): 0.92, (26, 5): 0.99, (26, 6): 1.17, (26, 7): 0.80,
        (27, 1): 0.69, (27, 2): 0.74, (27, 3): 0.92, (27, 4): 1.05, (27, 5): 1.00,
        (28, 1): 1.16, (28, 2): 0.82, (28, 3): 0.90, (28, 4): 1.04,
        (29, 1): 0.56, (29, 2): 0.67, (29, 3): 0.80, (29, 4): 0.78, (29, 5): 0.77, (29, 6): 0.95, (29, 7): 0.95, (29, 8): 1.03,
        (30, 1): 0.95, (30, 2): 1.14, (30, 3): 0.73, (30, 4): 0.92, (30, 5): 0.97,
    }

    # 归一化标签数据
    normalized_label_data, min_val, max_val = normalize_label_data(label_data)

    # 绘制平滑的压力曲线，传入 min_val 和 max_val
    dforcesme, dforcesme0 = plot_smoothed_pressure(label_data, normalized_label_data, cmap, min_val, max_val)

    # 绘制小提琴图并添加数据点
    plot_violin(label_data, cmap)

if __name__ == "__main__":
    main()
