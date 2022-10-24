function time_idx = compute_time_CC(all_h) 
% Based on a h-mask, compute connected components of 1's accounting for
% channel topography. 
%
% Inputs
%   all_h: (n_chan, nT_data) ;
% 
% Outputs
%   time_idx:   time index of the cluster including all input channels. 

[nchan, nT] = size(all_h) ; 
h_tot = sum(all_h, 1) ; 
time_idx = [] ; 
%time_idx = find(h_tot==nchan) ; 

CC = bwconncomp(all_h, conndef(2,'minimal')) ;
n_CC = CC.NumObjects ;
idx_clusters = CC.PixelIdxList ; 

      
for i_cc=1:n_CC
    curr_idx = idx_clusters{i_cc} ; 
    % Find indices along each dimension
    [Is, Js] = ind2sub([nchan, nT], curr_idx) ;
    channel_ko = setdiff(1:nchan,Is) ; 
    
    mat_tmp  =zeros(nchan,nT) ; mat_tmp(curr_idx) = 1 ;  
    if isempty(channel_ko)
        time_idx = [time_idx, unique(reshape(Js, 1, length(Js)))] ; 
    end

end
time_idx = sort(time_idx) ; 

end

