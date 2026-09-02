# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 12:40:09 2024

@author: PierreA.DERRIEN
"""

import pandas as pd

def zscore(data:pd.DataFrame, data_column:str, time_column:str, window:list) -> tuple:
    """
    Calculate the Z-score for a given column in a DataFrame over a specified time window.
    
    Parameters:
    data (pd.DataFrame): The input DataFrame.
    data_column (str): The name of the column to calculate the Z-score for.
    time_column (str): The name of the column containing the time values.
    window (list): A list containing the start and end times of the window.
    
    Returns:
    pd.DataFrame: The DataFrame with the Z-score column added.
    float: The mean of the data in the specified window.
    float: The standard deviation of the data in the specified window.
    """
    data = data.copy()
    mean= data[(data[time_column] > window[0])&(data[time_column] < window[1])][data_column].mean()
    std= data[(data[time_column] > window[0])&(data[time_column] < window[1])][data_column].std()
    data["Zscore"] = (data[data_column]-mean)/std
    
    return data , mean , std

def zscorewOoutlier(data:pd.DataFrame, data_column:str, outlier_threshold:float=95) -> tuple:
    """
    Calculate the Z-score for a given column in a DataFrame over a specified time window, excluding outliers.
    
    Parameters:
    data (pd.DataFrame): The input DataFrame.
    data_column (str): The name of the column to calculate the Z-score for.
    outlier_threshold (float): The percentage of outliers to exclude (default is 95).
    
    Returns:
    pd.DataFrame: The DataFrame with the Z-score column added.
    float: The mean of the data in the specified window, excluding outliers.
    float: The standard deviation of the data in the specified window, excluding outliers.
    """
    data = data.copy()
    upthreshold_value = pd.Series(data[data_column]).quantile(outlier_threshold/100)
    bottomthreshold_value = pd.Series(data[data_column]).quantile((100 - outlier_threshold)/100)
    meandata = data[(data[data_column] > bottomthreshold_value)&(data[data_column] < upthreshold_value)][data_column]
    mean= meandata.mean()
    std= meandata.std()
    zscores = (data[data_column]-mean)/std
    data["Zscore"] = zscores
    return data , mean , std



def zscore_trials(data:pd.DataFrame ,psth_window:list) -> pd.DataFrame:
    """
    Calculate the Z-score for each trial in a DataFrame over a specified time window.

    Parameters:
    data (pd.DataFrame): The input DataFrame containing trial data.
    psth_window (list): A list containing the start and end times of the window for Z-score calculation.
    
    Returns:
    pd.DataFrame: The DataFrame with the Z-score columns added.
    """

    data = data.copy()
    data["Time(s)"] = data.index
    all_zscores = pd.DataFrame()
    for trial in data.columns:
        if trial == "Time(s)":
            continue
        zscored_trial, mean, std = zscore(data, trial, "Time(s)", psth_window)
        all_zscores[trial] = zscored_trial["Zscore"]
    
    return all_zscores