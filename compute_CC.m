function [chan_idx, time_idx] = compute_CC(all_h, nT_data, all_names, all_t) 
% Based on a h-mask, compute connected components of 1's accounting for
% channel topography. 
% Called in TSL_plot_IO_fit.m
%
% Inputs
%   all_h: (n_chan, nT_data) ;
% 
% Outputs
%   chan_idx:   cell containing the channels indices for each cluster found. 
  

% Convert n_chan indices to 2D matrix accounting for topo adjacency
% using <= 4 neighbors for each electrode

%n_chan = size(all_h,1) ; 
%figure(); imagesc(all_h) ; set(gca,'YDir','Normal') ; colormap(jet) ; colorbar ; pause ; 

order_CC = 0 ; % using t-stats instead of number of elements
if nargin==4
    order_CC = 1  ;
end
if order_CC
    [all_t,all_h,~,idx_elec_mat] = convert_chans_to_topo_matrix(...
        all_t, all_h, all_names) ;
else
    [~,all_h,~,idx_elec_mat] = convert_chans_to_topo_matrix(...
        all_h, all_h, all_names) ;
end
n_chan_mat = size(all_h,1) ; 
CC = bwconncomp(all_h, conndef(3,'minimal')) ; %bwconncomp(all_h, conndef(3,'minimal')) ; % maximal
n_CC = CC.NumObjects ;
idx_clusters = CC.PixelIdxList ; 

chan_idx = cell(1, n_CC) ; 
time_idx = cell(1, n_CC) ; 

CC_size = zeros(1, n_CC) ; 


%test_tmp = zeros(n_chan, nT_data) ; 
      
for i_cc=1:n_CC
    curr_idx = idx_clusters{i_cc} ; 
    if order_CC
        CC_size(i_cc) = sum(abs(all_t(curr_idx))) ; 
    else
        CC_size(i_cc) = length(curr_idx) ; % number of elements in cluster
    end

    % Find indices along each dimension
    [Is, Js, Ktmp] = ind2sub([n_chan_mat, n_chan_mat, nT_data],curr_idx) ;
    % Convert the Is and Js (in artificial channel matrix)
    % into linear indices for the channels
    idx_chans = idx_elec_mat(sub2ind([n_chan_mat, n_chan_mat], Is, Js)) ; 
    chan_idx{i_cc} = unique(idx_chans) ; 
    time_idx{i_cc} = unique(Ktmp) ; 
    % take ALL the times included in this cluster, no matter the nb of electrodes
    
    %test_tmp(idx_chans, Ktmp) = i_cc ; 
end
%figure(); imagesc(test_tmp) ; set(gca,'YDir','Normal') ; colormap(jet(n_CC+1)) ; colorbar ; pause ; 

[~, I_sizes] = sort(CC_size, 'descend') ; 
chan_idx = chan_idx(I_sizes) ; 
time_idx = time_idx(I_sizes) ; 

end

