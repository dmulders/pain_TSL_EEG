function [chan_idx_pos, all_pts_pos] = keep_positive_clusters(chan_idx, all_pts,time_vec)
% Onle keep cell entries in chan_idx and all_pts when clusters are across
% positive times.

[~,idx_0] = min(abs(time_vec-0)) ;
n_cells = length(all_pts) ; 

chan_idx_pos = {} ; 
all_pts_pos = {} ; 

for i_c=1:n_cells
    curr_t = all_pts{i_c} ; 
    if any(curr_t>idx_0)
        chan_idx_pos{end+1} = chan_idx{i_c} ; 
        all_pts_pos{end+1} = all_pts{i_c} ; 
    end
end


end

