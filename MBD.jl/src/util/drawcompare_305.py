import pandas as pd
import matplotlib.pyplot as plt


project_name="demo305"

data_jl = pd.read_csv(f"{project_name}/data.csv")
data_jl = data_jl.drop(data_jl.index[0])

data_matlab = pd.read_csv(f"{project_name}/matlab.csv")
data_matlab = data_matlab.drop(data_matlab.index[0])

data_adams = pd.read_csv(f"{project_name}/posveldoor.tab", sep='\t')
data_adams = data_adams.rename(columns={'.2.Last_Run.request_POS_doorcm.TIME': 't'})


def map_series_to_dataframes(series_list, column_info):
    """
    Maps each series to its corresponding data in different dataframes.
    
    :param series_list: A list of names for the data series (e.g., ['x_series', 'y_series']).
    :param column_info: A dictionary mapping each series to its column names in different data sources.
                        Example: {'x_series': {'adams': 'adams_x_column', 'jl': 'jl_x_column'}, ...}
    :return: A dictionary mapping each series to its relevant data for plotting.
    """
    mapped_series = {}
    markers = ['*','o', 'p']  # Markers for plotting 
    sources = ['adams', 'matlab','jl']  # Data sources 

    for series in series_list:
        series_info = []
        for source, marker in zip(sources, markers):
            series_info.append({
                'dataframe': globals()[f"data_{source}"],  # Accessing DataFrame by its name
                'column_name': column_info[series][source],
                'series_label': f"{series}_{source}",
                'marker_style': marker
            })
        mapped_series[series] = series_info

    return mapped_series

# List of series names
series_names = ['x_series','vx_series', 'y_series', 'z_series']

# Mapping of series to their respective column names in different data sources
column_mapping = {
    'x_series': {'adams': '.MODEL_CLOSE_ST2.Last_Run.request_POS_doorcm.U2','matlab':'x1', 'jl': 'x1'},#
    'vx_series': {'adams': '.MODEL_CLOSE_ST2.Last_Run.request_POS_doorcm.U6','matlab':'x1d', 'jl': 'xd1'},#
    'y_series': {'adams': '.MODEL_CLOSE_ST2.Last_Run.request_POS_doorcm.U3','matlab':'y1', 'jl': 'y1'},#
    'z_series': {'adams': '.MODEL_CLOSE_ST2.Last_Run.request_POS_doorcm.U4','matlab':'z1', 'jl': 'z1'} #
}

# Creating the series mapping
series_mapping = map_series_to_dataframes(series_names, column_mapping)

def plot_data_series(series_mapping, project_name):
    """
    Plots a comparison of different data series.
    
    :param series_mapping: Dictionary containing the series and their data for plotting.
    :param project_name: Project name, used for saving the plot.
    """
    for series, data_info in series_mapping.items():
        plt.figure(figsize=(16, 10))
        for data in data_info:
            y_values = data['dataframe'][data['column_name']]
            y_min = y_values.min() - 10
            y_max = y_values.max() + 10
            plt.plot(data['dataframe']['t'], data['dataframe'][data['column_name']], 
                     label=data['series_label'], marker=data['marker_style'])
        plt.ylim([y_min, y_max])
        plt.legend()
        plt.title(f"Comparison of {series}")
        plt.savefig(f"{project_name}/compare-{series}.png")
        plt.close()

# Plotting each series
plot_data_series(series_mapping, project_name)