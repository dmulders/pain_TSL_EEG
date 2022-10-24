function V_wave = get_waveform_TCS2(NI_params, target_temp, base_temp, pulse_dur, ...
    foreperiod, postperiod, active_zones)
% Returns the intensity (V) for one TCS2 pulse. 
% 
% Inputs
%   target_temp     targeted temperature.
%   pulse_dur       duration of the laser pulse.
%   SR              sampling rate
%   foreperiod      recording period before the stim starts
%   postperiod      recording period after the end of the stim (eg to
%                   record an input such as RT). 
%                   ATTENTION: with the NI, the time window of the
%                   recordings is the same for the inputs and outputs.
%
% Outputs
%   V_wave      (n_tot_bins, 8): voltage waveform for the NI. 
%               col 1: neutral temperature.
%               col 2:6: temperature of the 5 zones.
%               col 7: neutral temperature trigger. 
%               col 8: temperature trigger.
%
trigger_EEG = 0 ; % use the neutral temp trigger to define a single onset trigger
% KO on the EEG ANT systems, needs a continuous trigger

if nargin<7
   active_zones = [1,2,3,4,5] ;  
end
if nargin<6
    postperiod = 0 ; 
end
if nargin<5
    foreperiod = 0 ; 
end

if target_temp>60
    disp(['Setting the temp to 60° instead of ', num2str(target_temp)])
    target_temp = 60 ; 
end

% ==== (0) ==== Parameters, prepare voltage matrix
SR = NI_params.Rate ; 
n_bins = SR*pulse_dur ;
fore_bins = foreperiod*SR ;
post_bins = postperiod*SR ;
n_tot_bins = n_bins+fore_bins+post_bins ; 

V_wave = zeros(n_tot_bins, 8) ;

% voltage for neutral tempearture
neutral_volt = NI_params.V_intercept + (base_temp*NI_params.V_slope) ;
active_volt = NI_params.V_intercept + (target_temp*NI_params.V_slope) ;

% ==== (1) ==== Temperature waveforms
% temperature per zone: either base_temp or target_temp
V_pulses = neutral_volt*ones(1,5) ; 
V_pulses(active_zones) = active_volt ; 
V_pulses = repmat(V_pulses,n_bins,1) ;

V_wave(:,1:6) = neutral_volt ; 
% Pulse
V_wave((fore_bins+1):(end-post_bins),2:6) = V_pulses;

% ==== (2) ==== Trigger waveforms
% neutral_temperature_trigger stays to 0
V_wave(:,7) = 0 ;
if trigger_EEG
    % ==*== Use the neutral temperature trigger port (7) to send single 
    %       triggers to the EEG system.
    n_bins_trigger = 2 ; % if SR = 512 Hz, need at least 2 bins for trigger
    volt_mat((fore_bins+1):(fore_bins+1+n_bins_trigger),7) = 5 ;
    volt_mat(1:n_bins_trigger,7) = 5 ;
end

% active_temperature_trigger w/o stimulating during for and postperiod
V_wave((fore_bins+1):(end-post_bins),8) = NI_params.trigger_voltage ;
% Set the active temperature trigger, also during for and
% postperiod to stimulate at neutral temperature
%V_wave(:,8) = NI_params.trigger_voltage ;


% add 0 and the end to bring the voltage back to 0 (on the NI pins)
% (otherwise, voltage stays the same as imposed in the last bin)
V_wave(end,:) = 0 ;

end
