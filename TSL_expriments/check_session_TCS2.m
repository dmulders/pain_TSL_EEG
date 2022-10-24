function check_session_TCS2(NI_session, pulse_dur, temp_base,...
    temp_I1, temp_I2, active_zones, fn_check, subj_name, record_RT)
% Check session: sequence of stimuli of 2 intensities, subject has to
% evaluate the intensity received and RT is recorded. 

% Parameters %%%%%%%%%%%%%%%%%%%%
p           = [0.5,0.5] ;       % 
trial       = 15 ;              %
use_foreground = 0 ;            % NEVER foreground, can bug
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<9
    record_RT = 0 ; 
end
if record_RT
    foreperiod  = 0.2 ;         % to see the laser temperature feedback 
                                % before stimulation. 
    postperiod  = 2 ;           % to be able to record inputs (eg RT) after the stim
else
    foreperiod  = 0.2 ;         %
    postperiod  = 0.2 ;         %
end
SR = NI_session.Rate ; 

% ======================================================================= %
% == (1) Generate the sequence of stimuli.
% ======================================================================= %
% ==*== Stimuli intensities
s = GenRandSeq(trial, p) ;
% s:    sequence with entries in {1,2} of length = # trials = trial.

save(fn_check, '-v7.3', 's', 'subj_name','temp_I1','temp_I2') ; 
% '-v7.3': to be able to save .mat file > 2 GB

% ======================================================================= %
% == (2) Get output waveforms
% ======================================================================= %

% Intensity (V) waveforms
wave_I1 = get_waveform_TCS2(NI_session, temp_I1, temp_base, pulse_dur,...
    foreperiod, postperiod, active_zones) ; 
wave_I2 = get_waveform_TCS2(NI_session, temp_I2, temp_base, pulse_dur,...
    foreperiod, postperiod, active_zones) ; 

% ======================================================================= %
% == (3) Start the sequence
% ======================================================================= %
disp(['-- Press enter to start the sequence'])
pause ; 
result = repmat(struct(),trial,1) ; 
% to store the temperature feedback, RT, ratings of I
rated_I = zeros(1, trial) ; 
all_RT = -1.*ones(1, trial) ; 

for idx_stim=1:trial
    curr_I = s(idx_stim) ; 
    
    % =================================================================== %
    % Laser stimulation
    % =================================================================== %
    if curr_I==1
        outputData = wave_I1 ;
    else %curr_I==2
        outputData = wave_I2 ;
    end
    disp(['-- Press enter to deliver stimulus ',num2str(idx_stim),...
        '/',num2str(trial),': ',num2str(curr_I)])
    pause ; 
    
    [data,time] = stimulate_TCS2(NI_session, outputData, use_foreground) ;
    if record_RT
        RT = get_RT_from_YAP_button(data(:,1), data(:,2), SR) ; 
        all_RT(idx_stim) = RT ; 
    end
      
    result(idx_stim).data = data;
    result(idx_stim).time = time;
    
    disp(['  (H: high; L: low)'])
    disp(['RT = ', num2str(all_RT(idx_stim))]) ; 
    rated_str = input('Intensity of the received stimulus: H or L? ','s') ; 
    if strcmpi(rated_str,'l')
        rated_I(idx_stim) = 1 ; 
    end
    if strcmpi(rated_str,'h')
        rated_I(idx_stim) = 2 ; 
    end
    
    % figure
    % plot(time,data) ; legend('temp','RT')    
    
    save(fn_check, '-v7.3', 's', 'subj_name','result','rated_I','all_RT',...
        'idx_stim','trial', 'temp_I1','temp_I2') ; 
end

% nb of false reports
n_miss = sum(rated_I~=s) ; 
disp(['* Nb of mistakes to evaluate the intensity: ',num2str(n_miss)])


disp(['[s; rated_I; all_RT]: '])
I_estI = [reshape(s,1,trial); rated_I]
disp([' ------------ End check ------------ '])

end

