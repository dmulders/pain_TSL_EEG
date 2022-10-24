function RT = get_RT_from_YAP_button(data, data_trigger, SR)
% Fetch the reaction time (RT) from the NI input data from the YAP response
% button. 

% Inputs
%   data:      YAP RT trigger (data is a column vector). 
%   data_trigger:   trigger from the TCS2 (red plug). 


% data: is always O and then 5V when button is pressed.  --> TO CHECK?!
rt_trig = find(data>4,1) ;
rt_stim = find(data_trigger>3,1) ; 

if isempty(rt_trig)
    RT = -1 ;
else
    if isempty(rt_stim)
        disp('No trigger received for the stimulus onset') ; 
        RT = -1 ; 
    else
        RT = (rt_trig-rt_stim)/SR ;
    end
end

end
