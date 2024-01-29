import pandas as pd
import matplotlib.pyplot as plt


project_name="demo205"

data_jl = pd.read_csv(f"{project_name}/data.csv")

dfmt = pd.read_csv(f"{project_name}/t.csv")
dfmx = pd.read_csv(f"{project_name}/x1.csv")
dfmy = pd.read_csv(f"{project_name}/y1.csv")
dfmz = pd.read_csv(f"{project_name}/z1.csv")
data_matlab = pd.concat([dfmt, dfmx,dfmy,dfmz], axis=1)

data_adams = pd.read_csv(f"{project_name}/sph_pla.tab", sep='\t')
data_adams = data_adams.rename(columns={'.sphpla.Last_Run.A.TIME': 't'})

def map_series_to_dataframes(series_list, column_info):
    """
    Maps each series to its corresponding data in different dataframes.
    
    :param series_list: A list of names for the data series (e.g., ['x_series', 'y_series']).
    :param column_info: A dictionary mapping each series to its column names in different data sources.
                        Example: {'x_series': {'adams': 'adams_x_column', 'jl': 'jl_x_column'}, ...}
    :return: A dictionary mapping each series to its relevant data for plotting.
    """
    mapped_series = {}
    markers = ['*', 'o','p']  # Markers for plotting
    sources = ['adams', 'jl','matlab']  # Data sources

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
series_names = ['x_series', 'y_series', 'z_series']

# Mapping of series to their respective column names in different data sources
column_mapping = {
    'x_series': {'adams': '.sphpla.Last_Run.A.X', 'jl': 'x1','matlab':'x1'},
    'y_series': {'adams': '.sphpla.Last_Run.A.Y', 'jl': 'y1','matlab':'y1'},
    'z_series': {'adams': '.sphpla.Last_Run.A.Z', 'jl': 'z1','matlab':'z1'}
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
            plt.plot(data['dataframe']['t'], data['dataframe'][data['column_name']], 
                     label=data['series_label'], marker=data['marker_style'])
        plt.legend()
        plt.title(f"Comparison of {series}")
        plt.savefig(f"{project_name}/compare-{series}.png")
        plt.close()

# Plotting each series
plot_data_series(series_mapping, project_name)