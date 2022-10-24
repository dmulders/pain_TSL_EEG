function [n_data, n_data_conf, est_conf_tmp, est_p1_tmp, RW_p1_tmp, ...
    out_conf, trials_resp_TP, trials_resp_conf, s, ...
    all_RT, all_RT_conf, p1_mean, out_conf_full] = ...
    load_extract_delta(fn_test, remove_missed_resp, learned_param, a_LR)
% Load the rating data and extract useful information + fit using delta
% rule (Rescorla-Wagner model). 
%
% When behavioral questions are asked to the subject, we ask them to
% estimate:
% - proba of next stimulus
% - confidence in this estimation.


% ===== * ===== Load useful subjective data
test_data = load(fn_test) ; 
% save(fn_session, '-v7.3', 's', 'trials_resp', 'tot_trial_time', ...
% 'timings', 'tot_run_time','subj_name','foreperiod','postperiod',...
% 'laser_data','proba_same','all_RT','missed_resp','prob_type',...
% 'vas','all_start_times','yCenter','xCenter','rect', 'tstart',...
% 'conf_proba','missed_conf', 'all_RT_conf') ; 
% prob_type(idx_stim) = 0 ; if curr_I == 2, 1 otherwise

s = test_data.s ; 
trials_resp = test_data.trials_resp ; 
missed_resp = test_data.missed_resp ; 
missed_conf = test_data.missed_conf ; 

proba_same = test_data.proba_same ; 
prob_type = test_data.prob_type ; 
conf_proba = test_data.conf_proba ; 

% to check corr btw RT and conf
all_RT = test_data.all_RT ; 
all_RT_conf = test_data.all_RT_conf ; 

trials_resp_TP = trials_resp ; trials_resp_conf = trials_resp ; 
if remove_missed_resp
    % Remove the behavioral data when it was not
    % confirmed by the subject. 
    missed_asked_qu = missed_resp(trials_resp) ; 
    trials_resp_TP = trials_resp(~missed_asked_qu) ; 
    % only keep un-missed responses

    missed_asked_conf = missed_conf(trials_resp) ; 
    trials_resp_conf = trials_resp(~missed_asked_conf) ; 
end
n_data = length(trials_resp_TP) ; 
n_data_conf = length(trials_resp_conf) ; 

% Only keep data of interest
proba_same = proba_same(trials_resp_TP) ; 
conf_proba = conf_proba(trials_resp_conf) ;    

all_RT = all_RT(trials_resp_TP) ; 
all_RT_conf = all_RT_conf(trials_resp_conf) ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ===== * ===== Model inference (delta rule - RW)
[p1_mean, p1_update, confidence] = RescorlaWagner(s, learned_param, a_LR) ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RW_p1_tmp       = make_row_vector(p1_mean(trials_resp_TP)) ;
RW_update_tmp   = make_row_vector(p1_update(trials_resp_conf)) ; 
RW_conf_tmp     = make_row_vector(confidence(trials_resp_conf)) ; 

est_p1_tmp      = zeros(1,n_data) ; 
est_conf_tmp    = conf_proba ; 

for i_data = 1:n_data
    % select optimal quantities of interest
    i_trial = trials_resp_TP(i_data) ; 
    if prob_type(i_trial)==0 % last stim received: I2
        est_p1_tmp(i_data) = 100 - proba_same(i_data) ;  
    else
        est_p1_tmp(i_data) = proba_same(i_data) ;
    end
end

out_conf_full = confidence ; 
out_conf = RW_conf_tmp ;
%out_conf_full = p1_update ;
%out_conf = RW_update_tmp ; 

end
