from pathlib import Path
#Setup custom library path
import sys
#IF custom library path is not in the sys.path, add it/them 
lib_paths = []
sys.path.extend(lib_paths)

#Generals imports
import pandas as pd
import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter, butter, filtfilt

#Fiber specific imports
from custom_Zscore import zscore
from trials_extraction import extract_trials
from photobleaching_correctors import local_min_photobleaching_correction as photobleach_corrector




def ProcessData(data:pd.DataFrame, SignalsNames:list[str], DictEvents:dict[str, list[float]], plots=True, # Extract PSTH around stim events
    PSTH_window: list = [1, 2], #seconds
    ZscoreWind: list = [-1,0], #seconds
    Savgol: bool = True,
    savgolLong: bool = True, 
    LowPassFilter: bool = False,
    FilterCutoff: int = 100,
    eventColor: str = None, 
    signals_colors: list = None,
    MedianFilter: bool = False,
    NormQuantile: bool = False, 
    returnUnprocessed: bool = False, 
    SizeSavgol: str = "10sec",
    SavgolTime: int = 5
    )-> tuple[pd.DataFrame, dict, dict, pd.DataFrame, pd.DataFrame]:
    """
    Process fluorescence signals for a dataset by applying photobleaching correction,
    optional smoothing/normalization steps, and event-aligned trial extraction.

    Parameters
    ----------
    data : pandas.DataFrame
        Time-indexed signal dataframe. Each column contains a raw fluorescence trace.
    SignalsNames : list[str]
        Names of the raw signal columns to process.
    DictEvents : dict[str, list[float]]
        Event times grouped by event type. Each list contains timestamps used to
        build event-aligned trials.
    plots : bool, default=True
        If True, display intermediate plots of corrected signals and PSTHs.
    PSTH_window : list[float], default=[1, 2]
        Window (in seconds) used around each event for trial extraction.
    ZscoreWind : list[float], default=[-1, 0]
        Time window used to compute z-score baselines for each trial.
    Savgol : bool, default=True
        If True, apply a Savitzky-Golay filter to the raw signal before
        photobleaching correction with a 5ms window to avoid oscillation induced by amplitude modulation.
    savgolLong : bool, default=True
        If True, compute a long-timescale Savitzky-Golay trend and subtract it from
        the corrected signal. This is used to separate slow trends from event-related responses.
    LowPassFilter : bool, default=False
        If True, apply a low-pass Butterworth filter after bleaching correction.
    FilterCutoff : int, default=100
        Cutoff frequency (Hz) for the low-pass filter.
    eventColor : str or None, default=None
        Color map used to draw event markers in plots.
    signals_colors : list[str] or None, default=None
        Colors used when plotting each signal.
    MedianFilter : bool, default=False
        If True, subtract a rolling median from the corrected signal.
    NormQuantile : bool, default=False
        If True, normalize signals by a 1st/99th percentile scaling.
    returnUnprocessed : bool, default=False
        If True, also return the unprocessed median-corrected signal stack.
    SizeSavgol : str, default="10sec"
        Window size label used for the long Savitzky-Golay detrending step.
    SavgolTime : int, default=5
        Time window (in milliseconds) used for the short Savitzky-Golay filter applied
        

    Returns
    -------
    tuple
        (Corrected_signals, trialsDataset, zscoredDataset, rolling_median,
         Uncorrected_signals)
        where:
        - Corrected_signals : processed signal dataframe
        - trialsDataset : event-aligned trial data for each event type and signal
        - zscoredDataset : z-scored version of each trial dataset
        - rolling_median : rolling median used for detrending (or None)
        - Uncorrected_signals : unprocessed copy of the corrected signals when
          returnUnprocessed is True, otherwise an empty DataFrame
    """
    if data.isnull().values.any():
        raise ValueError("Input data contains NaN values. Please clean the data before processing.")
    Window5ms = data.index[(data.index > 1.000)&(data.index < (1.000+0.001*SavgolTime))].shape[0] # 1 milliseconds window

    length10sec = data.index[(data.index >1.)&(data.index < 11.)].shape[0]
    length30sec = data.index[(data.index >1.)&(data.index < 31.)].shape[0]
    for wind in [Window5ms, length10sec, length30sec]:
        if wind % 2 == 0:
            wind += 1  # Ensure the window length is odd for Savitzky-Golay filter
    #Preprocesss signals with photobleaching correction
    Corrected_signals = pd.DataFrame(index=data.index)
    Uncorrected_signals = pd.DataFrame(index=data.index)
    for signal in SignalsNames:
        #print(f"Photobleaching correction for signal: {signal}")
        if Savgol:
            data[signal] = savgol_filter(data[signal], window_length=Window5ms ,polyorder=1)
        dF, dFF, (Fit_type, control_arti) = photobleach_corrector(data.reset_index(), signal, "Time(s)")
        Corrected_signals[signal.replace("raw", "corrected")] = dFF["df_over_f"].values
    #Frequency over time analysis + plot
    # downsample to 1ms => 

    if LowPassFilter:
        # #Butterworth filter (10 Hz)
        fs = 1/np.mean(np.diff(data.index))  # Sampling frequency
        nyq = 0.5 * fs
        normal_cutoff = FilterCutoff / nyq
        b, a = butter(N=4, Wn=normal_cutoff, btype='low', analog=False)
        for signal in SignalsNames:
             signal = signal.replace("raw", "corrected")
             Corrected_signals[signal] = filtfilt(b, a, Corrected_signals[signal])

    if MedianFilter:
        window_size = int(5 / np.mean(np.diff(data.index)))  # 20 seconds window
        for signal in SignalsNames:
            signal = signal.replace("raw", "corrected")
            rolling_median = Corrected_signals[signal].rolling(window=window_size, center=True).median()
            rolling_median.fillna(method='bfill', inplace=True)
            rolling_median.fillna(method='ffill', inplace=True)
            Uncorrected_signals[signal] = Corrected_signals[signal]
            Corrected_signals[signal] = Corrected_signals[signal] - rolling_median
            #Fill NaN values at the beginning and end after rolling

    if savgolLong:
        SavgolData = pd.DataFrame(index=data.index)
        for signal in SignalsNames:
            signal = signal.replace("raw", "corrected")
            if SizeSavgol == "10sec":
                SavgolData[signal] = savgol_filter(Corrected_signals[signal], window_length=length10sec ,polyorder=1)
            elif SizeSavgol == "30sec":
                SavgolData[signal] = savgol_filter(Corrected_signals[signal], window_length=length30sec ,polyorder=1)
       
    if NormQuantile:
        #99th percentile normalize signals
        for signal in SignalsNames:
            signal = signal.replace("raw", "corrected")
            perc99 = np.percentile(Corrected_signals[signal], 99)
            perc1 = np.percentile(Corrected_signals[signal], 1)
            Corrected_signals[signal] = (Corrected_signals[signal] - perc1) / (perc99 - perc1)
    

    if plots:           
        #plot corrected signals
        for i, signal in enumerate(SignalsNames):
            plt.figure(figsize=(30,10))
            plt.plot(data[signal], label=signal, color="grey", linewidth=0.25, linestyle='--')
            sb.lineplot( x= data.index, y= control_arti, label = "Local Minima Fit", color = "black", linewidth=0.25, linestyle='--')
            signal = signal.replace("raw", "corrected")
            plt.plot(Corrected_signals[signal], label=signal, color=signals_colors[i], linewidth=0.25)
            if savgolLong:
                plt.plot(SavgolData[signal], label=signal + " Savgol", color="red", linewidth=0.25)
                plt.plot(Corrected_signals[signal] - SavgolData[signal], label=signal + " after Savgol", color="orange", linewidth=0.25)
            for EventsName, CleanedEvents in DictEvents.items():
                for x in CleanedEvents:
                    plt.axvline(x=x, color=eventColor[EventsName], linestyle='--', label=EventsName if x == CleanedEvents[0] else "", alpha = 0.5)   
            sb.despine()
            plt.xlabel("Time (s)")
            plt.title("Corrected Signals")
            plt.legend()
            plt.show()
    if savgolLong:
       Corrected_signals = Corrected_signals - SavgolData
    trialsDataset = {}
    zscoredDataset = {}
    for EventsName, CleanedEvents in DictEvents.items():
        # Extract PSTH around stim events
        Corrected_signals["CleanedEvents"] = 0
        for event_time in CleanedEvents:
            # Find the nearest index position to each event time
            nearest_idx = Corrected_signals[Corrected_signals.index >= event_time].index.min()
            Corrected_signals.at[nearest_idx, "CleanedEvents"] = 1
        trials_data = {}
        zscored_trials_data = {}
        for signal in SignalsNames:
            signal = signal.replace("raw", "corrected")
            trials_data[signal] = extract_trials(Corrected_signals, signal, "CleanedEvents", winsize =PSTH_window, compact= True)
            zscored_trials_temps = pd.DataFrame()
            for trial in trials_data[signal].columns:
                zscoredvalues, _, _ = zscore(trials_data[signal].reset_index(), trial, "Time", window=ZscoreWind)
                zscoredvalues.index = trials_data[signal].index
                zscored_trials_temps[trial] = zscoredvalues["Zscore"]
            zscored_trials_data[signal.replace("corrected", "zscored")] = zscored_trials_temps

        if plots:
            print(f"Plotting PSTH for event: {EventsName}")
            #plot PSTH for df/f and zscored
            fig, axs = plt.subplots(2,len(SignalsNames), figsize=(12,8), sharex=True)
            for i, signal in enumerate(SignalsNames):
                signal = signal.replace("raw", "corrected")
                if len(SignalsNames) > 1:
                    ax = axs[0, i]
                else:
                    ax = axs[0]
                for trial in trials_data[signal].columns:
                    sb.lineplot(data=trials_data[signal][trial], ax=ax, color="black", legend=False, alpha=0.5)
                sb.lineplot(data=trials_data[signal].mean(axis=1), ax=ax, color=signals_colors[i], label="Mean dF/F")
                ax.set_title(f"PSTH of {signal} (df/f)")
                ax.set_ylabel("df/f")
                ax.axvline(x=0, color='black', linestyle='--', label='Stimulation', alpha = 0.5)
                if len(SignalsNames) > 1:
                    ax2 = axs[1, i]
                else:
                    ax2 = axs[1]
                for trial in zscored_trials_data[signal.replace("corrected", "zscored")].columns:
                    sb.lineplot(data=zscored_trials_data[signal.replace("corrected", "zscored")][trial], ax=ax2, color="black", alpha=0.5, legend=False)
                sb.lineplot(data=zscored_trials_data[signal.replace("corrected", "zscored")].mean(axis=1), ax=ax2, color=signals_colors[i], label="Mean Zscore")
                ax2.set_title(f"PSTH of {signal.replace('corrected', 'zscored')} (Zscore)")
                ax2.set_ylabel("Zscore")
            plt.show()

        trialsDataset[EventsName] = trials_data
        zscoredDataset[EventsName] = zscored_trials_data
    
    if not MedianFilter:
        rolling_median = None
    if returnUnprocessed:
        return( Corrected_signals, trialsDataset, zscoredDataset ,rolling_median, Uncorrected_signals)
    else:
        return (Corrected_signals, trialsDataset, zscoredDataset , rolling_median, pd.DataFrame())


def SavePreprocess(output_directory: str, AnimalID: str, Corrected_signals: pd.DataFrame, trialsDataset: dict, zscoredDataset: dict, rolling_median: pd.DataFrame = None)-> None:
    output_dir = Path(output_directory)
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    print(f"Saving Corrected signals",)
    # Save Corrected Signals
    corrected_file = output_dir / f"{AnimalID}_CorrectedSignals.npy"
    np.save(corrected_file, Corrected_signals)
    
    print(f"Saving Trials Dataset")
    # Save Trials Dataset
    for event_name, trials_data in trialsDataset.items():
        for region, trials in trials_data.items():
            print(f"Saving Trials Dataset - Event: {event_name}, Region: {region}")
            trials_region_file = output_dir / f"{AnimalID}_{region}_{event_name}.csv"
            trials.to_csv(trials_region_file)
        
    print(f"Saving Zscored Dataset",end="\r")
    # Save Zscored Dataset
    for event_name, zscored_data in zscoredDataset.items():
        for region, zscored_trials in zscored_data.items():
            print(f"Saving Zscored Dataset - Event: {event_name}, Region: {region}")
            zscored_region_file = output_dir / f"{AnimalID}_{region}_{event_name}.csv"
            zscored_trials.to_csv(zscored_region_file)
        
    # Save Rolling Median if available
    if rolling_median is not None:
        print(f"Saving Rolling Median")
        rolling_median_file = output_dir / f"{AnimalID}_RollingMedian.csv"
        rolling_median.to_csv(rolling_median_file)
    print("Finished saving all preprocessed data.")