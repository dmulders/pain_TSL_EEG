function count_ep_kept()
% Check the nb of epochs kept for the analysis (after epoch rejections by
% visual inspection, amplitude criterion & possibly removing 10 first ep).

% Nb of ep kept saved in TSL_analyze_EEG.m w/ save_ep_kept = 1.

folder_fn = './results_EEG_lica/data/' ; 
%folder_fn = './results_EEG_lica_wo10/data/' ; 

fn_ep_kept = [folder_fn, 'Ep_kept_31subj_LP10_maxA80.mat'] ; 
%fn_ep_kept = [folder_fn, 'Ep_kept_31subj_LP30_maxA80.mat'] ; 
%fn_ep_kept = [folder_fn, 'Ep_kept_31subj_maxA150.mat'] ; 

data_ep = load(fn_ep_kept) ; 
all_ep_kept = data_ep.all_ep_kept ; % n_stim_test, n_conds_max, n_subj
all_s = data_ep.all_s ; % n_stim_test, n_conds_max, n_subj
n_subj = size(all_ep_kept,3) ; 

n_I1 = NaN(n_subj,1) ; 
n_I2 = NaN(n_subj,1) ; 
for i_subj=1:n_subj
    ep_tmp = squeeze(all_ep_kept(:,:,i_subj)) ; 
    s_tmp = squeeze(all_s(:,:,i_subj)) ;
    n_I1(i_subj) = sum(ep_tmp(s_tmp==1)) ; 
    n_I2(i_subj) = sum(ep_tmp(s_tmp==2)) ; 
end

n_I1_I2 = [n_I1, n_I2]
mean_n_I1_I2 = mean(n_I1_I2)
std_n_I1_I2 = std(n_I1_I2)

nep_rej = size(all_ep_kept,1)*size(all_ep_kept,2)*n_subj - sum(n_I1) - sum(n_I2)

end

