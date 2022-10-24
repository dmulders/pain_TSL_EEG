function [s, trials_resp] = load_extract_sequence_info(fn_test)
% Load the useful information about the sequence of stimuli. 
% Similar to load_extract_ratings.m, but without loading the ratings and
% without fitting any model. 
% 
% s: sequence with entries in {1, 2}.
% trials_resp: vector of the indices of the behavioral questions. 


% ===== * ===== Load useful data
test_data = load(fn_test) ; 

s = test_data.s ; 
trials_resp = test_data.trials_resp ; 


end
