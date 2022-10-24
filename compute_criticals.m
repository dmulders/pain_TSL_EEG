function criticals = compute_criticals(perm_cluster_t, n_betas) 
% Compute the criticals based on the distributions in perm_cluster_t. 
%
% perm_cluster_t:  (n_perm, n_betas, 2, 2)
   
cluster_threshold = 95 ; 

% indices:
%   n_betas
%   2: I1 and I2
%   2: mean cluster, max cluster
criticals = zeros(n_betas, 2, 2) ; 

% =============================================================== %
% ======== Compute criticals
% =============================================================== %
% criticals: (n_betas, 2, 2) ;
for i_beta=1:n_betas
    for i_I=1:2 
        criticals(i_beta,i_I,1) = prctile(...
            squeeze(perm_cluster_t(:,i_beta,i_I,1)), cluster_threshold) ; % mean cluster stat
        criticals(i_beta,i_I,2) = prctile(...
            squeeze(perm_cluster_t(:,i_beta,i_I,2)), cluster_threshold) ; % max cluster stat
    end
end



end

