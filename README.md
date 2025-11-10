# BurstingDynamicsInference

Contains the code used for the model inference, analysis and fitting in the manuscript "Enhancer control of promoter activity and variability via frequency modulation of clustered transcriptional bursts" by
Jana Tünnermann, Gregory Roth, Julie Cramard, and Luca Giorgetti.

## Usage

- **data_analysis**: functions to estimate from the experimental burst tracks the burst and interburst duration survival curves, correlation, etc
- **StoThyLiveCell**: functions to calculate analystically the observatbles from a multi-state model.
- **model_fit_2s2r**: functions to run the double regime, two state model to the experimental data.
- **model_fit_3s1r**: functions to run the three state model to the experimental data.
- **model_simulation**: functions to simulate the the double regime, two state model.
- **model_bootstrap**: functions to run a boostrap parametric bootstrap analysis.
- **modelassumption_robustness**: functions to run a robustness analysis of the interburst timescales and interburst correlation against missed detection of weak bursts.
