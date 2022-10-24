function [idx_select_mod, n_tot_epoch, vec_n_epochs] = find_idx_conds(...
    lwdata_subj, str_to_find)
% Function used to find the index of conditions with names containing the 
% input string 'str_to_find'. 


n_mod = size(lwdata_subj,2);
vec_n_epochs = [] ;                     % vector containing the nb of epoch 
                                        % for each modality found.
idx_select_mod = [] ;                   % idx for the select_mod 

for i = 1:n_mod
    curr_header = lwdata_subj(1,i).header ;
    curr_header_name = curr_header.name ; 
    
    curr_str_idx_mod = strfind(curr_header_name, str_to_find); 
    
    % only do something when the modality is the one selected
    if ~isempty(curr_str_idx_mod)
        vec_n_epochs = [vec_n_epochs, curr_header.datasize(1,1)] ; 
        idx_select_mod = [idx_select_mod, i] ;
    end    
end
n_tot_epoch = sum(vec_n_epochs) ;       % total nb of epochs for 
                                        % modality select_mod.   

end

