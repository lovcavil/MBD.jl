import pandas as pd
import numpy as np
from scipy.fft import fft, ifft
import plotly.graph_objs as go

scaling_dict = {
    'Strain@B1_1.RN_6': 0.339,
    'Strain@B2_1.RN_6': 0.3505,
    'Strain@B2_2.RN_6': 0.5762,
    'Strain@B2_3.RN_6': -0.5072,
    'Strain@B2_4.RN_6': 0.0763,
    'Strain@B3_1.RN_6': 0.8693,
    'Strain@B3_2.RN_6': 0.8576,
    'Strain@B3_3.RN_6': -0.2255,
    'Strain@B3_4.RN_6': 0.2455,
    
    # Add other columns and their scaling factors here
}

def load_data(path1, path2):
    df1 = pd.read_csv(path1)
    df2 = pd.read_csv(path2)
    return df1, df2

def filter_data(df1, df2, start_time1, end_time1, start_time2, end_time2, selected_columns1, selected_columns2, data_amount):
    filtered_df1 = df1.loc[(df1['t '] >= start_time1) & (df1['t '] <= end_time1)]
    filtered_df2 = df2.loc[(df2['Time'] >= start_time2) & (df2['Time'] <= end_time2)]
    selected_columns1 = selected_columns1[:data_amount]
    selected_columns2 = selected_columns2[:data_amount]
    return filtered_df1, filtered_df2, selected_columns1, selected_columns2

def apply_scaling(filtered_df1, selected_columns1):
    for column in selected_columns1:
        scale = scaling_dict.get(column, 1)
        filtered_df1.loc[:, column] = filtered_df1.loc[:, column] / scale
        print(scale)
    return filtered_df1



def apply_smoothing(filtered_df1, filtered_df2, selected_columns1, selected_columns2, smoothing_method, window_sizes1, window_sizes2):
    if smoothing_method == 'Rolling Mean':
        for column, window_size in zip(selected_columns1, window_sizes1):
            filtered_df1.loc[:, column] = filtered_df1.loc[:, column].rolling(window=window_size, min_periods=1).mean()
        for column, window_size in zip(selected_columns2, window_sizes2):
            filtered_df2.loc[:, column] = filtered_df2.loc[:, column].rolling(window=window_size, min_periods=1).mean()
    elif smoothing_method == 'Fourier Filter':
        for column, window_size in zip(selected_columns1, window_sizes1):
            signal = filtered_df1[column].values
            fft_vals = fft(signal)
            fft_vals[window_size:] = 0  # Zero out high frequencies
            filtered_signal = np.real(ifft(fft_vals))
            filtered_df1.loc[:, column] = filtered_signal
        for column, window_size in zip(selected_columns2, window_sizes2):
            signal = filtered_df2[column].values
            fft_vals = fft(signal)
            fft_vals[window_size:] = 0  # Zero out high frequencies
            filtered_signal = np.real(ifft(fft_vals))
            filtered_df2.loc[:, column] = filtered_signal
    return filtered_df1, filtered_df2

def create_figure(filtered_df1, filtered_df2, selected_columns1, selected_columns2, start_time1):
    fig = go.Figure()
    for column in selected_columns1:
        fig.add_trace(go.Scatter(x=filtered_df1['t '] - start_time1, y=filtered_df1[column], mode='lines', name=f'{column} (Data Source 1)'))
    for column in selected_columns2:
        fig.add_trace(go.Scatter(x=filtered_df2['Time'], y=filtered_df2[column], mode='lines', name=f'{column} (Data Source 2)'))
    fig.update_layout(title="Plot of Selected Columns vs Time", xaxis_title="Time", yaxis_title="Value")
    return fig

def create_figureacc(filtered_df1, filtered_df2, selected_columns1, selected_columns2, start_time1):
    fig = go.Figure()
    for column in selected_columns1:
        fig.add_trace(go.Scatter(x=filtered_df1['t '] - start_time1, y=filtered_df1[column], mode='lines', name=f'{column} (Data Source 1)'))
    for column in selected_columns2:
        fig.add_trace(go.Scatter(x=filtered_df2['Time'], y=filtered_df2[column]/ 1000 /9.8 , mode='lines', name=f'{column} (Data Source 2)'))
    fig.update_layout(title="Plot of Selected Columns vs Time", xaxis_title="Time", yaxis_title="Value")
    return fig
