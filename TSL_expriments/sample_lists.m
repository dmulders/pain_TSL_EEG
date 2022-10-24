function [trials_resp,tot_trial_time,timings,tot_run_time] = sample_lists(...
    trial,isi,resptime, IRI)
% Define when the behavioral questions are asked to the subjects and
% compute the timeline of the experiment. 
%
% Inputs
%   trial   total nb of trials in a sequence.
%   isi     (sec) vector of possible ISI between each pair of stimuli
%           within the sequence.
%   resptime    (sec). Additional time assigned to the subjects for the 
%               behavioral questions. 
%
% Outputs
%   trials_resp     indices of the trials when the behavioral questions are
%                   asked. 
%   tot_trial_time  (trials x 1) vector of time per trial. 
%   timings     (trials x 3) matrix with, for each trial: 
%               [total trial time, trial isi, question time].%               
%   tot_run_time    (trials x 1) vector of total time up to trial i. 

% ======================================================================= %
% (1) Define the trials when the behavioral questions should be asked.
% ======================================================================= %

if nargin<4
    % Questions asked to the subjects every 15\pm 3 trials.
    IRI = 12:18 ; % inter-response-interval
end
n_IRI = length(IRI) ; 

max_nb_qu = ceil(trial/min(IRI)) ;    % max nb of questions within the sequence
trials_resp = zeros(1,max_nb_qu) ; 

trial_idx_qu = 0 ;  % idx of the trial when a question is asked
idx_qu = 0 ;        % idx of the asked question

while trial_idx_qu<trial
    curr_IRI = IRI(randi(n_IRI,1)) ; 
    trial_idx_qu = trial_idx_qu + curr_IRI ; 
    idx_qu = idx_qu + 1 ;
    trials_resp(idx_qu) = trial_idx_qu ;      
end
trials_resp = trials_resp(1:idx_qu) ; 
trials_resp(trials_resp>trial) = [] ; 
% ======================================================================= %
% (2) Compute the timing of the trials.
% ======================================================================= %
% === Initializations === 
% Duration attributed to the responses for all the trials (0 when no
% question, resptime otherwise).
resp_dur_list = zeros(trial,1);
resp_dur_list(trials_resp)= resptime;

tot_trial_time = zeros(trial,1) ; 
timings = zeros(trial,3) ; 
tot_run_time = zeros(trial,1) ; 

n_isi = length(isi) ; 

for t = 1:trial    
    %isi=Shuffle(isi);
    curr_isi = isi(randi(n_isi)) ; 
    tot_trial_time(t) = (curr_isi + resp_dur_list(t) );
    
    timings(t,:)= [ tot_trial_time(t), curr_isi, resp_dur_list(t) ]; 
    tot_run_time(t) = sum( tot_trial_time(1:t) ) ;    
end

end