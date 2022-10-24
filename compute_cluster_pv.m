function [p_beta_CP, p_beta_CP_max] = compute_cluster_pv(all_t_beta, all_p_beta, ...
    perm_cluster_t, criticals, all_names, n_perm, n_betas, n_chan, nF_data, nT_data, alpha) 
% Compute the distribution of the cluster-level statistics based on 
% shuffled data. 
%
% Inputs
%   all_t_beta: (n_betas, n_chan, nF_data, nT_data, 2) ;
%   idx_elec_mat: (9,9)
%   perm_cluster_t: (n_perm, n_betas, 2, 2) ; 
%   criticals:(n_betas, 2, 2) ; last index: mean or max cluster stat
% 

p_beta_CP = ones(n_betas, n_chan, nF_data, nT_data, 2) ; 
p_beta_CP_max = ones(n_betas, n_chan, nF_data, nT_data, 2) ; 
% (n_betas, n_chan, nF_data, nT_data, 2) ; 

h_mask = double(all_p_beta<=alpha) ; 
all_t = all_t_beta.*h_mask ; 

  
for i_beta = 1:n_betas
    for i_I=1:2                
        curr_t = shiftdim(all_t(i_beta,:,:,:,i_I),1) ; % n_chan, nF_data, nT_data
        
        if n_chan>1
            % Convert n_chan indices to 2D matrix accounting for topo adjacency
            % using <= 4 neighbors for each electrode
            [curr_t,~,~,idx_elec_mat] = convert_chans_to_topo_matrix(...
                curr_t, curr_t, all_names) ;
            curr_t(isnan(curr_t)) = 0 ; % don't take in clusters
            n_chan_mat = size(curr_t,1) ; 
            
            CC = bwconncomp(curr_t, conndef(3,'minimal')) ;
        else
            curr_t = squeeze(curr_t) ;
            CC = bwconncomp(curr_t,4) ;
        end
        
        n_CC = CC.NumObjects ;
        idx_clusters = CC.PixelIdxList ; 
        
        tmp_p_mean = ones(n_chan, nF_data, nT_data); %test_tmp = zeros(n_chan, nF_data, nT_data) ; 
        tmp_p_max = ones(n_chan, nF_data, nT_data);        
        for i_cc=1:n_CC
            curr_idx = idx_clusters{i_cc} ; 
            v = sum(abs(curr_t(curr_idx))) ; 
            
            %figure('Name', ['i_cc = ', num2str(i_cc), ' pv = ', num2str(sum( abs(v)<perm_cluster_t(:,i_beta,i_I,1) )/n_perm)]); 
            %histogram(squeeze(perm_cluster_t(:,i_beta,i_I,1))); hold on ; 
            %yL = get(gca,'YLim'); plot([v, v], yL, 'r-') ; hold on ; plot([1, 1].*criticals(i_beta,i_I,1), yL, 'b-') ;
            %pause ; 
            % =========================================================== %
            % Found significant cluster (based on mean shuffled clusters)
            % =========================================================== %
            if abs(v)>=criticals(i_beta,i_I,1) % mean permut cluster
                pval = sum( abs(v)<perm_cluster_t(:,i_beta,i_I,1) )/n_perm ; 
                
                if n_chan>1
                    % Find indices along each dimension
                    [Is, Js, Fs, Ts] = ind2sub([n_chan_mat, n_chan_mat,nF_data, nT_data],curr_idx) ;
                    % Convert the Is and Js (in artificial channel matrix)
                    % into linear indices for the channels
                    idx_chans = idx_elec_mat(sub2ind([n_chan_mat, n_chan_mat], Is, Js)) ; 
                    idx_final = sub2ind([n_chan, nF_data, nT_data], idx_chans,Fs, Ts) ; 
                    tmp_p_mean(idx_final) = pval ;
                    %test_tmp(idx_final) = i_cc ; 
                else
                    tmp_p_mean(curr_idx) = pval ;
                    %test_tmp(curr_idx) = i_cc ; 
                end
                                 
            end
            
            % =========================================================== %
            % Found significant cluster (based on max shuffled clusters)
            % =========================================================== %
            if abs(v)>=criticals(i_beta,i_I,2) % max permut cluster
                pval = sum( abs(v)<perm_cluster_t(:,i_beta,i_I,2) )/n_perm  ;
                if n_chan>1
                    % Find indices along each dimension
                    [Is, Js, Fs, Ts] = ind2sub([n_chan_mat, n_chan_mat,nF_data, nT_data],curr_idx) ;
                    idx_chans = idx_elec_mat(sub2ind([n_chan_mat, n_chan_mat], Is, Js)) ; 
                    idx_final = sub2ind([n_chan, nF_data, nT_data], idx_chans,Fs, Ts) ; 
                    tmp_p_max(idx_final) = pval ;
                else
                    tmp_p_max(curr_idx) = pval ; 
                end
            end
        end
        p_beta_CP(i_beta,:,:,:,i_I) = tmp_p_mean ;                
        % p_beta_CP: (n_betas, n_chan, nF_data, nT_data, 2) ;
        p_beta_CP_max(i_beta,:,:,:,i_I) = tmp_p_max ; 
        % p_beta_CP_max: (n_betas, n_chan, nF_data, nT_data, 2) ; 
        
        %figure(); imagesc(test_tmp) ; set(gca,'YDir','Normal') ; colormap(jet(n_CC+1)) ; colorbar ; 
        %pause ; 
    end
end



end

