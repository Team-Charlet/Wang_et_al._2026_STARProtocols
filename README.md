# Wang_et_al._2026_STARProtocols

Code for fiber photometry analysis of DORIC recordings, with preprocessing, photobleaching correction, event alignment, and z-score analysis.

## Overview

This repository contains a lightweight analysis pipeline for fiber photometry recordings. It is designed around DORIC export files and supports:

- reading DORIC : .csv  and .doric (HDMF5 based) files
- photobleaching correction using exponential or local-minimum fits
- optional smoothing and normalization
- event-triggered trial extraction around behavioral events
- z-scoring relative to a baseline window
- saving processed signal and trial outputs

The main workflow is implemented in [Preprocess.py](Preprocess.py), with support functions in [doric_system_file.py](doric_system_file.py), [photobleaching_correctors.py](photobleaching_correctors.py), [trials_extraction.py](trials_extraction.py), and [custom_Zscore.py](custom_Zscore.py).

## Repository structure

- [Preprocess.py](Preprocess.py): end-to-end preprocessing pipeline
- [photobleaching_correctors.py](photobleaching_correctors.py): photobleaching correction functions
- [trials_extraction.py](trials_extraction.py): event-aligned trial extraction
- [custom_Zscore.py](custom_Zscore.py): z-score calculation utilities
- [doric_system_file.py](doric_system_file.py): DORIC file reader utilities
- [ExampleProcessing.ipynb](ExampleProcessing.ipynb): example notebook showing the intended workflow
- [ExampleProcessing.md](ExampleProcessing.md): fixed, online-readable version of the example workflow
- [requirements.txt](requirements.txt): Python package requirements

## Requirements
Create the environment from the shared Conda definition:

```bash
conda env create -f environment.yml
conda activate STAR_Env
```

To update an existing environment after changing the dependencies, run:

```bash
conda env update -f environment.yml --prune
```

The project expects common scientific packages such as pandas, numpy, scipy, matplotlib, seaborn, and scikit-learn.

## Typical input format

The pipeline is intended for DORIC CSV exports where the first row is metadata and the time column is named `Time(s)`.

Your DataFrame should look conceptually like this:

- index: time in seconds (`Time(s)`)
- signal columns: raw fluorescence channels (for example `AIn-1 - Dem (AOut-1)` or `CeA_raw`)
- event columns: digital inputs like `DI/O-1`, `DI/O-2`, etc.

A common pattern is to rename your raw fluorescence input and map relevant DI/O channels to event names such as:

- `StartVideo`
- `FootShock`
- `SoundCue`

## Quick start

The notebook in [ExampleProcessing.ipynb](ExampleProcessing.ipynb) gives the canonical workflow. The key steps are reproduced below.

```python
from Preprocess import ProcessData
import pandas as pd

filepath = "Path/to/your/file.csv"  # replace with your actual file path
Data = pd.read_csv(filepath, skiprows=1).set_index("Time(s)")
```

### 1) Identify event channels

```python
availableDIO = ["DI/O-1", "DI/O-2", "DI/O-3"]
DIO_names = {
    "DIO01": "StartVideo",
    "DIO02": "FootShock",
    "DIO03": "SoundCue",
}

EventStart = {}
dios = []
for i, dio in enumerate(availableDIO):
    event_name = DIO_names.get(f"DIO{str(i+1).zfill(2)}", dio)
    dio_events = Data[Data[dio].diff() == 1].index
    EventStart[event_name] = dio_events.tolist()
    dios.append(event_name)

Data.rename(columns={dio: DIO_names.get(f"DIO{str(i+1).zfill(2)}", dio) for i, dio in enumerate(availableDIO)}, inplace=True)
```

### 2) Define the signal to process

```python
Signals = ["CeA_raw"]
Data.rename(columns={"AIn-1 - Dem (AOut-1)": "CeA_raw"}, inplace=True)
Data = Data.dropna(subset=Signals)[Signals + dios]
```

### 3) Run preprocessing

```python
Preprocessed = ProcessData(
    Data,
    Signals,
    EventStart,
    plots=True,
    PSTH_window=[5, 30],      # seconds around each event
    ZscoreWind=[-2.5, 0],     # baseline window for z-scoring
    savgolLong=False,
    LowPassFilter=False,
    eventColor={
        "StartVideo": "blue",
        "FootShock": "red",
        "SoundCue": "orange",
    },
    signals_colors=["green"],
    NormQuantile=True,
)
```

This returns a tuple containing:

- corrected signals
- event-aligned trial datasets
- z-scored datasets
- rolling median information (if used)
- unprocessed data (optional)

## Main functions

### ProcessData

The main entry point is `ProcessData` in [Preprocess.py](Preprocess.py). It performs the full analysis workflow:

- validates the input data
- applies photobleaching correction
- optionally smooths or filters the signal
- computes a long-timescale detrending step
- normalizes values if requested
- extracts trials around each event time
- computes z-scored PSTHs relative to a baseline window
- optionally plots signal traces and PSTHs

Key arguments:

- `data`: pandas DataFrame indexed by time
- `SignalsNames`: list of fluorescence channels to process
- `DictEvents`: dictionary of event names and event times
- `plots`: whether to display plots
- `PSTH_window`: window used to extract trial segments around events
- `ZscoreWind`: time window used for baseline z-scoring
- `Savgol`: apply Savitzky-Golay smoothing
- `LowPassFilter`: apply low-pass filtering
- `MedianFilter`: subtract rolling median
- `NormQuantile`: percentile normalization
- `returnUnprocessed`: optionally return unprocessed traces

### Photobleaching correction

The correction utilities are in [photobleaching_correctors.py](photobleaching_correctors.py). They support:

- exponential fitting of bleaching trend
- local-minimum-based bleaching correction
- calculation of `df/fFrozen`/`df_over_f` style outputs

This step is essential before event analysis because photobleaching can create slow drifts that obscure real signals.

### Trial extraction

The helper in [trials_extraction.py](trials_extraction.py) extracts signal traces centered on each event. Each trial is stored as a column in a DataFrame, with time relative to the event onset.

### Z-score computation

The custom z-scoring utilities in [custom_Zscore.py](custom_Zscore.py) compute baseline-normalized values using a time window defined by the user. This is useful for comparing responses across trials and events.

## Saving outputs

The project includes a helper function to save processed data:

```python
from Preprocess import SavePreprocess

SavePreprocess(
    output_directory="output_folder",
    AnimalID="Animal_01",
    Corrected_signals=Preprocessed[0],
    trialsDataset=Preprocessed[1],
    zscoredDataset=Preprocessed[2],
    rolling_median=Preprocessed[3],
)
```

This saves:

- corrected signals as `.npy`
- per-event trial tables as `.csv`
- z-scored tables as `.csv`
- rolling median data if available

## Recommended workflow

1. Import your DORIC CSV file.
2. Set the time index to `Time(s)`.
3. Identify your fluorescence channel(s).
4. Identify event timestamps from DI/O channels.
5. Rename channels to clear names.
6. Run `ProcessData`.
7. Inspect plots and adjust parameters such as `PSTH_window`, `ZscoreWind`, and filtering settings.
8. Save the processed outputs.

## Notes

- The example notebook is the best reference for actual usage.
- Event times should be aligned to the precise onset of the behavior or stimulus being analyzed.
- If signals are noisy or unstable, start by disabling aggressive filtering and checking the plots before tightening the preprocessing pipeline.
- This repository is research-oriented and is best used as a flexible analysis template rather than as a fully packaged GUI pipeline.

## License

This project is distributed under the repository license included in the root directory.
