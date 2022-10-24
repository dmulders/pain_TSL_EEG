function h_mask = significance_FDR_corrected(p_values, alpha_level)
% Indicates which null hypothesis are rejected based on the matrix of p
% values p_values, with the FDR bound alpha_level in [0,1].
% NaN values in the matrix p_values are ignored (they are not considered as 
% hypotheses).
% 
% Reference: Genovese et al., 2002, "Thresholding of statistical maps in functional
% neuroimaging using the False Discovery Rate". 
%
% Out: 
%   h_mask:     a binary matrix of the same size as p_values, with 1
%               indicating that the corresponding null hypothesis is rejected
%               with the maximum FDR of alpha_level on avg.

h_mask = zeros(size(p_values)) ;% h=1 means that H0 is rejected
                                % h=0 means that we cannot reject H0

p_values_vec = p_values(:) ; 
nb_H0 = sum(~isnan(p_values_vec)) ; % nb of null hypotheses

% sort in ascending order, with NaN at the end of the sort
[p_values_sorted, lin_idx] = sort(p_values_vec) ; 
% lin_idx(1:nb_H0) = indices of the non-nan pval in ascending order 

c_V = sum( 1./([1:nb_H0]) ) ; % applies for any joint distribution of the p-values
% if indep pv: can use c_V = 1 cfr Genovese2002
% 4.7439 with 64 chan

all_thresholds = ([1:nb_H0].*alpha_level)./(nb_H0.*c_V) ; 
all_thresholds = reshape(all_thresholds, size(p_values_vec)) ; % ensure col or row, idem as p-values

% find largest p-value which is ok
idx_pv = find(p_values_sorted<=all_thresholds, 1, 'last') ;
h_mask(lin_idx(1:idx_pv)) = 1 ; 


end
