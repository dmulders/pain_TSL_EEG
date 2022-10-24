function perm_cluster_t = compute_cluster_perm_stats(beta_CP, ...
    alpha_level, all_names, n_perm, n_betas, n_chan, nF_data, nT_data, n_subj) 
% Compute the distribution of the cluster-level statistics based on 
% shuffled data. 
%
% beta_CP:  (n_betas, n_chan, nF_data, nT_data, n_subj, 2, n_perm)

paired_tests = 0 ; mean_H0 = 0 ; always_t_test = 1 ; two_sided = 1 ;    

% indices:
%   n_perm
%   n_betas
%   2: I1 and I2
%   2: mean cluster, max cluster
perm_cluster_t = zeros(n_perm, n_betas, 2, 2) ; 

% =============================================================== %
% ======== Compute significance of beta coefficients (shuffled)
% =============================================================== %
% beta_CP: (n_betas, n_chan, nF_data, nT_data, n_subj, 2, n_perm)
  
for i_perm=1:n_perm
    for i_beta = 1:n_betas
        for i_I=1:2                
            curr_betas = shiftdim(beta_CP(i_beta, :, :, :, :,i_I,i_perm),1) ;
            % (n_chan, nF_data, nT_data, n_subj)
            curr_betas = permute(curr_betas, [4, 1, 2, 3]) ;
            curr_betas = reshape(curr_betas, n_subj, n_chan*nF_data*nT_data) ;

            [~, h_mask, ~, ~, ~, t_stats] = ...
                compute_significance_matrix(curr_betas, alpha_level, 10, ...
                paired_tests, mean_H0, always_t_test, two_sided) ;
            
            t_stats = reshape(t_stats.*h_mask, n_chan, nF_data, nT_data) ;
            if n_chan>1
                % Convert n_chan indices to 2D matrix accounting for topo adjacency
                % using <= 4 neighbors for each electrode
                
                t_stats = convert_chans_to_topo_matrix(...
                    t_stats, zeros(size(t_stats)), all_names) ;
                
                t_stats(isnan(t_stats)) = 0 ; % don't take in clusters
                
                CC = bwconncomp(t_stats, conndef(3,'minimal')) ;
            else
                t_stats = squeeze(t_stats) ; 
                CC = bwconncomp(t_stats,4) ;
            end
            n_CC = CC.NumObjects ; 
            idx_clusters = CC.PixelIdxList ; 
            CC_sizes = zeros(n_CC,1) ; % cluster level statistics
            
            for i_cc=1:n_CC
                v = sum(abs(t_stats(idx_clusters{i_cc}))) ; 
                % take abs here otherwise mean cluster size doesn't make sense
                if v>0
                    CC_sizes(i_cc) = v ; 
                end
            end
            if isempty(CC_sizes)
                CC_sizes = 0 ; 
            end
            % curr_p_val = sum(abs(T_permut)>=abs(T_CC))./n_permut ;

            perm_cluster_t(i_perm, i_beta, i_I, 1) = mean(abs(CC_sizes)) ;% mean cluster stat
            perm_cluster_t(i_perm, i_beta, i_I, 2) = max(abs(CC_sizes)) ;% max cluster stat
        end
    end
end


end

