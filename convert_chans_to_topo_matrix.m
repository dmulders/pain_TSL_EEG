function [pval_elec_mat,hmask_elec_mat,names_elec_mat,idx_elec_mat] = ...
    convert_chans_to_topo_matrix(pval_array,hmask_array,all_names)
% Convert the input data, whose row indices correspond to channels from a
% 64 EEG cap, by reshaping this first dim into a 2-way matrix accounting
% for the electrodes adjacency on the head. 
%
% Inputs
%   pval_array:     (n_chan,m1,m2) matrix of pvalues
%   hmask_array:    (n_chan,m1,m2) matrix of corresponding h_mask
%                   (boolean entries, H0 rejected or not).
%   all_names:      n_chan-length cell indicating the channel names ordered
%                   as in the 2 previous inputs.
%
% Outputs
%   pval_elec_mat:  (9,9,m1,m2) matrix of the input pvalues arranged
%                   according to adjacency on the head.
%   hmask_elec_mat: (9,9,m1,m2) matrix of the inpu hmask values arranges
%                   according to 2D adjacency on the head.
%   names_elec_mat: (9,9) cell of strings indicating the corresponding
%                   channel names.
%   idx_elec_mat:   (9,9) matrix containing in each entry the index of the
%                   corresponding electrode in the original all_names 
%                   vector.
%
% Only works for n_chan=64 
% (positions of electrodes in 2D matrix are assigned by hand).


names_elec_mat = get_elec_mat() ; 
% Dimensions of 2D matrix for electrode positions
[n_rows_mat,n_cols_mat] = size(names_elec_mat) ; 
idx_elec_mat = zeros(n_rows_mat,n_cols_mat) ; 

[n1_pval,n2_pval,n3_pval] = size(pval_array) ; 
[n1_h,n2_h,n3_h] = size(hmask_array) ; 

pval_elec_mat = NaN*ones(n_rows_mat,n_cols_mat, n2_pval, n3_pval) ;
hmask_elec_mat = zeros(n_rows_mat,n_cols_mat, n2_h, n3_h) ;
n_chan = length(all_names) ; 
if n_chan~=64
    disp('... [compute_CC.m -- convert_chans_to_topo_matrix.m] Only the 64-chan 10-20 EEG system is considered...')
end
if n_chan~=n1_pval || n_chan~=n1_h
    error(['The channels names should correspond to the rows ',...
        'of pval_array and hmask_array'])
end

for idx_chan=1:n_chan
    idx_chan_in_mat = find(strcmpi(names_elec_mat,all_names{idx_chan})) ; 
    [Is,Js] = ind2sub([n_rows_mat,n_cols_mat], idx_chan_in_mat) ; 
    n_idx = length(idx_chan_in_mat) ; 
    for idx_tmp = 1:n_idx
        pval_elec_mat(Is(idx_tmp),Js(idx_tmp),:,:) =  pval_array(idx_chan,:,:) ; 
        hmask_elec_mat(Is(idx_tmp),Js(idx_tmp),:,:) =  hmask_array(idx_chan,:,:) ; 
        
        idx_elec_mat(Is(idx_tmp),Js(idx_tmp)) = idx_chan ; 
    end
end



end

function names_elec_mat = get_elec_mat()
% Returns a 9x9 matrix containg an electrode name in each entry, or NaN if
% no electrode is assigned to the entry.
names_elec_mat = ...
    {NaN, NaN, 'FP1', 'FPz','FPz','FPz','FP2',NaN,NaN;...
    NaN, 'AF7','AF3','AF3',NaN, 'AF4','AF4','AF8',NaN;...
    'F7','F5','F3','F1','Fz','F2','F4','F6','F8';...
    'FT7','FC5','FC3','FC1','FCz','FC2','FC4','FC6','FT8';...
    'T7','C5','C3','C1','Cz','C2','C4','C6','T8';...
    'TP7','CP5','CP3','CP1','CPz','CP2','CP4','CP6','TP8';...
    'P7','P5','P3','P1','Pz','P2','P4','P6','P8';...
    'PO7','PO5','PO3','PO3','POz','PO4','PO4','PO6','PO8';...
    NaN,'PO5','O1','Oz','Oz','Oz','O2','PO6',NaN;...
    } ;



end
