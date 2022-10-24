function [chan_order, order_desired] = define_chan_order_colors(chanlocs, desired_chans)
% Return the indices to use to obtain the channels in the order specified
% by desired_chans.
% The channels that are specified in desired_chans but cannot be found in
% chanlocs are ignored. 
%
% Out
%   chan_order  idx of the channels in chanlocs to consider to match
%               desired_chans. All channels from chanlocs are considered.
%   order_desired   same as chan_order but without including the indices of
%                   channels NOT mentioned in desired_chans. 

if nargin<2
    desired_chans = {'Fpz', 'Fz', 'FCZ', 'CZ', 'CPZ', 'PZ', 'POz', 'Oz', ...
        'O1', 'PO3', 'PO7', 'P7', 'P5', 'P3', 'P1', 'CP1', 'CP3', 'CP5', 'TP7', ...
        'T7', 'C5', 'C3', 'C1', 'FC1', 'FC3', 'FC5', 'FT7', 'F7', 'F5', 'F3', ...
        'F1', 'AF3', 'AF7', 'Fp1', 'Fp2', 'AF8', 'AF4', 'F2', 'F4', 'F6', 'F8',...
        'FT8', 'FC6', 'FC4', 'FC2', 'C2', 'C4', 'C6', 'T8', 'TP8', 'CP6', ...
        'CP4', 'CP2', 'P2', 'P4', 'P6', 'P8', 'PO8', 'PO4', 'O2'} ;
end

all_names = {chanlocs.labels} ; 
n_chans = length(chanlocs) ; 
n_desired = length(desired_chans) ; 
order_desired = NaN(1,n_desired) ;


for idx_desired = 1:n_desired
    curr_label = desired_chans{idx_desired} ; 
    curr_idx = find(strcmpi(all_names, curr_label)) ;
    
    if ~isempty(curr_idx)
        order_desired(idx_desired) = curr_idx ; 
    end
end
order_desired = order_desired(~isnan(order_desired)) ; 

idx_remain = setdiff(1:n_chans, order_desired) ; 
chan_order = [order_desired, idx_remain] ; 


end

