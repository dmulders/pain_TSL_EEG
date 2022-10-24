% test a call to exp_one_condition.m with dummy parameters

% test_no_LSD = 1 in exp_one_cond
% to check that the whole script can be ran without error (eg when saving data at the end)

fn_dir = ['./all_data_test/'] ; 
if ~exist(fn_dir, 'dir') % 0 or 7 if it exists
    mkdir(fn_dir)
end

cond_idx    = 0 ;   
subj_idx    = 0 ;
subj_name   = 'test' ; 
p_train     = [0.3,0.3] ; 
str_TP      = get_str_TP(p_train) ; 
IRI_train   = 2:3 ; % 12:18
n_stim_test = 10 ; 
pulse_dur   = 0.05 ; 
temp_base   = 20 ; 
temp_I1     = 50 ; 
temp_I2     = 60 ; 
V0_temp     = 20 ; 
V10_temp    = 70 ; 
foreperiod  = 0.2 ; 
postperiod  = 0.4 ; 

fn_train = [fn_dir, 'data_test', num2str(subj_idx), ...
    num2str(cond_idx), str_TP, '_',num2str(n_stim_test),'stim.mat'] ; 

%%%%%%%%%%%%%%%% TO DO

exp_one_condition(NaN, p_train, n_stim_test, pulse_dur, temp_base, ...
    temp_I1, temp_I2, V0_temp, V10_temp, fn_train, subj_name, foreperiod,...
    postperiod,'IRI',IRI_train) ; 