% Add paths to the different toolboxes at the beginning of a session

% VBA toolbox for Bayesian model selection. 
% https://mbb-team.github.io/VBA-toolbox/
addpath(genpath('./VBA-toolbox'));

% EEGLAB functions useful to do the scalp plots
addpath(genpath('./topoplots')) ;

% Perceptually uniform colormaps
addpath(genpath('./Colormaps')) ; 

% minimal TP model, from https://github.com/florentmeyniel/MinimalTransitionProbsModel, 
% with some updates to test variants of the initial models (with different priors, learning AF, ...)
addpath(genpath('./IdealObserversCode')) ; 

% solve low resolution on plot
set(0, 'DefaultFigureRenderer', 'painters') ;










