%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to employ for the sequence experiments: 
%       TCS2 with 2 NI devices.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Paremeters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Params indep. of stimuli %%%%%%
subj_idx    = 36 ;              %
subj_name   = 'XXX' ;           %
plus_sessions = 0 ;             % 0: cross, 1: plus, 2: levels. 
p_small = 0.3 ;  % 1/3          % TO TEST
all_TP_test = get_TPs_config(p_small,plus_sessions) 
n_conds     = size(all_TP_test,1) ; %
TP_order    = randperm(n_conds) ; %
TP_order2   = randperm(n_conds) ; %
TP_order_full = [TP_order, TP_order2] ; % do twice each condition
n_conds_full = length(TP_order_full) ; 
 
n_stim_test = 100 ;             % per session
n_stim_train= 50 ;              %
str_stim = '_TCS' ;             %
stim_type   = 2 ;               % type of stimuli
                                % 1: LSD, 2: TCS
TCS_USB     = 0 ;               % TCS: controlled through USB cable, otherwise 2 NI DAQ devices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def_temp_I1 = 15 ;              % temperature in [°C] for stim 1
def_temp_I2 = 58 ;              % temperature in [°C] for stim 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =*= All sessions (checks, staircase, train, test)
pulse_dur   = 0.25 ;             % duration of laser pulses [s]
temp_base   = 31 ;              % temperature in [°C] before/after the stim
SR          = 1000 ;            % samplingrate
active_zones= [1,2,3,4,5] ;     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =*= To record RT (staircase, check sessions)
time_RT = 2 ;                   % time to record RT from the stimulus trigger
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =*= Train and test sessions
foreperiod  = 0.2 ;             % no temperature imposed during for and postperiod (no stimulation trigger, cfr get_waveform_TCS2.m)
                                % before stimulation. 
postperiod  = 0.4 ;             % to be able to record inputs after the stim (& leave time before the beep if foreground stimulation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Only to call exp_one_condition (LSD stim)
V0_temp     = 20 ; 
V10_temp    = 70 ;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param_str = [str_stim, '_subj',num2str(subj_idx)] ; 

fn_dir = ['./all_data/'] ; 
if ~exist(fn_dir, 'dir')
    mkdir(fn_dir)
end

fn_cond_order = [fn_dir, 'Cond',param_str,'.mat'] ; 
save(fn_cond_order,'-v7.3', 'subj_name','all_TP_test','TP_order','TP_order_full') ; 


%% (0) Initialization
NI_TCS2 = initialize_TCS2_NI(SR) ;



%% (0b) test send one pulse -- can send a few pulses before starting the xp here
temp_test = 20 ;
pulse_dur = 0.25 ; 
fp_test = 0.2 ; % fore25period5
pp_test = 2 ; % postperiod
use_foreground = 0 ; 

record_RT = 1 ; 

wave_test = get_waveform_TCS2(NI_TCS2, temp_test, temp_base, pulse_dur,...
    fp_test, pp_test, active_zones) ; 

disp(['-- Press enter to deliver one stimulus of temp ',...
        num2str(temp_test)]) ; pause ; 

[data, time] = stimulate_TCS2(NI_TCS2, wave_test, use_foreground) ;    

if record_RT
    % in1: should be YAP response and in2: TCS2 trigger
    RT = get_RT_from_YAP_button(data(:,1), data(:,2), SR) ; 
end


%% (1) Staircase for the pain threshold, using reaction times
pulse_dur = 0.25 ;

fn_dir = './all_data/' ; % necessary eg if test_call_exp_one_cond has been called
fn_thresh = create_new_fn([fn_dir, 'Staircase',param_str]) ; 

results_stair = staircase_TCS2(NI_TCS2, pulse_dur, temp_base, active_zones, ...
    record_RT, fn_thresh) ;
save(fn_thresh, '-v7.3', 'subj_name','results_stair') ; 


temp_I1 = input('Enter the low temperature used: ') ; 
if isempty(temp_I1)
    temp_I1     = def_temp_I1 ;
end

temp_I2 = input('Enter the high temperature used: ') ; 
if isempty(temp_I2)
    temp_I2     = def_temp_I2 ;
end
save(fn_thresh, '-v7.3', 'subj_name','results_stair', 'temp_I2','temp_I1') ;

%% (2) check session
% //!\\ Can change temp_I1 and/or temp_I2 here!
temp_I1 = 20 ; 
temp_I2 = 57 ; 
pulse_dur = 0.25 ;

fn_dir = ['./all_data/'] ; 
fn_pre_check = create_new_fn([fn_dir, 'Pre_check',param_str]) ; 
% ATTENTION: used for supervised

disp(['--- Starting the pre-check session with I_1 = ',num2str(temp_I1),...
    ' and I_2 = ',num2str(temp_I2),' ---'])
check_session_TCS2(NI_TCS2, pulse_dur, temp_base, ...
    temp_I1, temp_I2, active_zones, fn_pre_check, subj_name, record_RT) ; 


%% (3) Training session
%p_train     = [0.7, 0.7] ; % alternation biased ;  % p_1from2, p_2from1
p_train     = [0.7, 0.4] ; % I_1 biased, avoid too much heat  % p_1from2, p_2from1 
%p_train     = [0.5, 0.5] ; 
str_TP      = get_str_TP(p_train) ; 
IRI_train   = 8:12 ;

temp_I1 = 20 ; 
temp_I2 = 57 ; % 52
pulse_dur = 0.25 ; % 0.5

fn_dir = './all_data/' ; 
fn_train = create_new_fn([fn_dir, 'Train',param_str, str_TP, '_',...
    num2str(n_stim_train),'stim']) ; 

exp_one_condition(NI_TCS2, p_train, n_stim_train, pulse_dur, temp_base, ...
    temp_I1, temp_I2, V0_temp, V10_temp, fn_train, subj_name, foreperiod,...
    postperiod,'IRI',IRI_train,'stim_type',stim_type,'TCS_USB',TCS_USB, ...
    'active_zones', active_zones) ; 
sca ; 
beep ; % end of sequence

%% (4) Test sessions

% Can record and note skin temperature, using Tempett: perpendicular to
% skin, not too far

IRI_test = 12:18 ;  
fn_dir = './all_data/' ; 

temp_I1 = 20 ; 
temp_I2 = 57 ; % 52
pulse_dur = 0.25 ;

for idx_cond = 1:n_conds_full %5:n_conds %1:n_conds   
    % Test session numbered idx_cond -> get the TPs
    ind_TP = TP_order_full(idx_cond) ; 
            
    curr_TP = all_TP_test(ind_TP,:) ; 
    str_TP = get_str_TP(curr_TP) ; 
    
    fn_test = create_new_fn([fn_dir, 'Test_cond',num2str(idx_cond), param_str, ...
        str_TP, '_',num2str(n_stim_test),'stim']) ; 
        
    disp(['[Cond ',num2str(idx_cond),'/',num2str(n_conds_full),']',...
        ' TPs: (', num2str(round(curr_TP(1),3)),', ',...
        num2str(round(curr_TP(2),3)),') == Press enter']) ; 
    pause ; 
    
    exp_one_condition(NI_TCS2, curr_TP, n_stim_test, pulse_dur, temp_base, ...
        temp_I1, temp_I2, V0_temp, V10_temp, fn_test, subj_name, foreperiod,...
        postperiod,'IRI',IRI_test,'stim_type',stim_type,'TCS_USB',TCS_USB, ...
        'active_zones', active_zones,'VAS_min_max',0) ;
    sca ; 
    beep ; % end of sequence
end

% Can record and note skin temperature, using Tempett: perpendicular too
% skin, not too far

%% (5) Post-check session
fn_dir = './all_data/' ; 
fn_post_check = create_new_fn([fn_dir, 'Post_check', param_str]) ;

disp(['--- Starting the post-check session with I_1 = ',num2str(temp_I1),...
    ' and I_2 = ',num2str(temp_I2),' ---'])
check_session_TCS2(NI_TCS2, pulse_dur, temp_base, ...
    temp_I1, temp_I2, active_zones, fn_post_check, subj_name,record_RT) ; 


