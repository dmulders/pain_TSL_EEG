function GFP = compute_GFP(scalp_data)
% Compute the Global Field Power (GFP) for the input (real) vector. 
%
% Cfr Skrandies1990_GFP or Khanna2015_microstates. 
% GFP: is helpful to identify overall ERP component from multichannel data.
% GFP = spatial standard deviation; quantifies the amount of activity at
% each time point in the field considering the data from all recording
% electrodes simultaneously resulting in a reference-independent descriptor
% of the potential field. 
% 
%
% Inputs
%   scalp_data:     input scalp potentials, as a (n_chan, n_time, ...) matrix. 

GFP = std(scalp_data, 0, 1) ; 

end
