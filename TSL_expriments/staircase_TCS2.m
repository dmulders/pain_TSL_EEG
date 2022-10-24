function results = staircase_TCS2(NI_session, pulse_dur, temp_base, ...
    active_zones, record_RT, fn_thresh)
% Staircase procedure to find the Adelta threshold with the TCS2 (NI). 
%
% The subject has to use the response button to indicate the detection of a
% stimulus. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
use_foreground = 0 ;            % usually best = 0 
% can have bug "Timeout expired before operation could complete." in foreground
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<5
    record_RT   = 0 ;        % use RT for the staircase (require a response button)
end
if record_RT
    foreperiod  = 0.2 ;         %
    postperiod  = 2 ;           % to be able to record inputs (eg RT) after the stim
else
    foreperiod  = 0.2 ;         %
    postperiod  = 0.2 ;         %
end

SR = NI_session.Rate ; 

init_temp = 45 ; % starts at this temperature
all_delta = [5, 2, 1, 0.1] ; 
n_delta = length(all_delta) ; 
idx_delta = 1 ;


all_detections = [] ;   % 0 = not detected
                        % 1 = detected
all_RT = []; 
threshold_found = 0 ;   
all_temp = init_temp ; % vector of tested temperatures
 
while threshold_found~=1
    curr_temp = all_temp(end) ; 
    
    % == (1) == Stimulate with all_temp(end)   
    curr_wave = get_waveform_TCS2(NI_session, curr_temp, temp_base, pulse_dur,...
        foreperiod, postperiod, active_zones) ; 
    disp(['-- Press enter to stimulate at ', num2str(curr_temp), ' °C'])
    pause ;
    % Random pause of a duration between 0 and 4 seconds
    duration_pause = 4*rand(1,1) ; % 1 + 4*rand(1,1)
    pause(duration_pause) ; 

    [data, time] = stimulate_TCS2(NI_session, curr_wave, use_foreground) ;
    if record_RT
        % in1: should be YAP response and in2: TCS2 trigger
        RT = get_RT_from_YAP_button(data(:,1), data(:,2), SR) ; 
        all_RT = [all_RT, RT] ; 
    end 
        
    % == (2) == Painful or not
    if record_RT
        if RT == -1
            % stimulus not detected 
            curr_detected = 0 ;  
            disp(['RT = -1, stimulus not detected'])
        elseif RT>0.65
            % RT slower than 650ms ==> no Adelta resp, cfr Mouraux2012
            curr_detected = 0 ; 
            disp(['RT = ',num2str(RT),': stimulus detected but temp not high enough ( RT',...
                ' > 0.650 s (no Adelta response)'])
        else        
            curr_detected = 1 ; 
            disp(['RT = ',num2str(RT),': stimulus detected '])
        end
        
    else    
        curr_detected = input(['Stimulus painful/detected ? type 0 for no; 1 for yes ']) ;    
    end
    
    all_detections = [all_detections, curr_detected] ;
    
    % == (3) == Threshold found?
    if length(all_detections)>=4
        if sum(abs(diff(all_detections((end-3):end))))==3
            % 3 consecutive changes
            threshold_found = 1  ; 
        end
    end
    
    % == (4) == Change the delta if search direction changed
    if length(all_detections)>1 && abs(diff(all_detections((end-1):end))) 
        % change detected <-> not detected 
        idx_delta = min(idx_delta + 1, n_delta) ; 
        % otherwise: keep the same delta
    end
    
    % == (5) == Compute the next temperature
    if all_detections(end)
        % stimulus painful --> decrease temperature
        next_temp = curr_temp - all_delta(idx_delta) ; 
    else
        % stimulus not painful --> increase temperature
        next_temp = curr_temp + all_delta(idx_delta) ; 
    end 
    % Don't go above 60deg
    next_temp = min(next_temp, 60) ; 
    
    all_temp = [all_temp, next_temp] ;
    
    % Save tmp data
    save(fn_thresh, '-v7.3', 'all_temp','all_RT', 'all_detections') ; 
end

all_temp = all_temp(1:(end-1)) ; 

% Threshold = mean of 4 temp which lead to 3 consecutive changes
threshold = mean(all_temp((end-3):end)) ; 
disp(['=== Threshold = ', num2str(threshold), '°C'])

if all_detections(end)
    disp(['=== Last painful temp = ', num2str(curr_temp), '°C'])
else
    % Last tested temp was not painful
    disp(['=== Last painful temp = ', num2str(all_temp(end-1)), '°C'])
end
disp([' ------------ End staircase ------------ '])

results.all_temp = all_temp ; 
results.all_detections = all_detections ; 
results.all_RT = all_RT ; 
results.threshold = threshold ; 

end






