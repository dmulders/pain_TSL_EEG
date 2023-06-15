# Confidence of probabilistic predictions modulates the cortical response to pain

Dounia Mulders, Ben Seymour, André, Mouraux, Flavia Mancini

Article published in [PNAS](https://www.pnas.org/doi/10.1073/pnas.2212252120). 

Associated [preprint](https://doi.org/10.1101/2022.08.11.503296)



![Image](https://user-images.githubusercontent.com/15798671/197595357-e06a29cf-adf8-48a3-83ce-a242a3203ab9.png)



## Description

This project provides the codes written to analyze the data of the above article and produce all figures shown in the manuscript. 

The folder /TSL_experiments/ contains the codes to generate and deliver the sequences of stimuli. 

Associated data are available on an [OSF repository](https://osf.io/8xvtg/). 

## Running

All codes are written in Matlab and were ran using Matlab R2019b.

The codes for the Bayesian models were written by [Florent Meyniel](https://github.com/florentmeyniel/MinimalTransitionProbsModel). They are provided in this repository in /IdealObserversCode/ with some updates to test variants of the initial models (with different priors, learning AF, ...). 

* To collect data and run experiments in the lab, you need 
    - a stimulator and DAQ device 
    - Matlab with the DAQ and [Psychtoolbox](http://psychtoolbox.org/)
    - all codes can be ran from run_all_stim_TCS2.m (check sessions, training, test sessions)

* To analyze the behavioral data 
    - behavioral data from the [OSF repository](https://osf.io/8xvtg/)  
    - always start by running add_all_paths_TSL.m that will add the required sub-folders to the Matlab path 
    - TSL_anayze_ratings.m: loads and analyzes the behavioral data of all the subjects, for one model and one parameter set. The path fn_dir should correspond to the behavioral data folder. 
    - TSL_fit_on_ratings.m: computes the fit of different models (with different parameters) and does the model comparison. 

* To analyze the EEG data 
    - EEG data from the [OSF repository](https://osf.io/8xvtg/).
    - always start by running add_all_paths_TSL.m that will add the required sub-folders to the Matlab path 
    - TSL_analyze_EEG.m: loads and analyzes the EEG recordings. The path fn_dir_EEG should correspond to the EEG data folder. Data are saved as specified in the function and can be reloaded and plotted using other functions. 
    - TSL_plot_avg_EEG.m: reloads useful data and displays the average EEG responses. Data must have been saved by running TSL_analyze_EEG.m with save_avg_eeg = 1 beforehand.  
    - TSL_plot_IO_fit.m: reloads useful data and displays the model fitting. Data must have been saved by running TSL_analyze_EEG.m with IO_fit_opt = 1 beforehand.  

* To perform the parameter recovery analysis, using codes from the folder /param_recovery/
    - start by running add_paths_recov.m to add the required folders to the Matlab path
    - simulate_behavior.m: simulates behavior using a range of parameters consistent with the ones observed in the original data set. 
    - fit_simulated_data.m: computes the quality of fit on data simulated in simulate_behavior.m.
    - disp_param_recovery.m: plots the outcomes of the parameter recovery analysis. The data saved in /data_simu/ enables producing the figures without re-computing the simulations. 

## Contact

You can contact me at dounia **dot** mulders **at** uclouvain.be for any question. :-)
