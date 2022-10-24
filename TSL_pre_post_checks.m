function TSL_pre_post_checks()
% Load and save summary info from the behavioral data of all the subjects. 
% 
% Count numbers of pre- and post-check mistakes. 

% ======================================================================= %
% ===== * ===== General parameters
% ======================================================================= %
all_subj_idx    = setdiff(1:36,[1, 11, 15, 28, 33]) ;
all_subj_idx    = 1:36 ; 

str_stim        = '_TCS' ;
n_stim_test     = 100 ; 
fn_dir          = ['./data_ratings/'] ; 
n_subj = length(all_subj_idx) ; 
n_checks = zeros(n_subj,2) ; % # of stim identified during pre & post check sessions
temps_used = zeros(n_subj,2) ; 

% ======================================================================= %
% ===== * ===== Load data
% ======================================================================= %
for i_subj = 1:n_subj
    real_subj_idx = all_subj_idx(i_subj) ; 
    
    subj_idx = num2str(real_subj_idx) ; 
    subj_idx_padd = num2str_padd(real_subj_idx, 2) ; 
    param_str = [str_stim, '_subj',subj_idx] ; 
    fn_dir_subj = [fn_dir, 'subj',subj_idx_padd,'/'] ; 
        
    % =============================================================== %
    % *** Pre_check...
    % =============================================================== %
    fn_pre = [fn_dir_subj,'Pre_check',param_str,'_1.mat'] ; 
    tmp_data = load(fn_pre) ;
    s = tmp_data.s ; 
    rated_I = tmp_data.rated_I ; 
    
    n_checks(i_subj,1) = sum(s==rated_I) ; 
            
    % =============================================================== %
    % *** Post_check...
    % =============================================================== %
    fn_post = [fn_dir_subj,'Post_check',param_str,'_1.mat'] ; 
    tmp_data = load(fn_post) ; 
    s = tmp_data.s ; 
    rated_I = tmp_data.rated_I ; 
    
    n_checks(i_subj,2) = sum(s==rated_I) ; 
    temps_used(i_subj,1) = tmp_data.temp_I1 ; 
    temps_used(i_subj,2) = tmp_data.temp_I2 ; 
    
end

all_ages = [21, 22, 28, 20, 22, 22, ...
    24, 21, 27, 20, 25, 23, ...
    22, 27, 39, 22, 23, 20, ...
    24, 26, 26, 18, 20, 28, ...
    30, 25, 26, 24, 25, 23, ...
    18, 21, 52, 27, 24, 18]' ; 

checks_T1_T2 = [15-n_checks, temps_used, all_ages] ;
IDs_checks_T1_T2 = [all_subj_idx', checks_T1_T2]
figure(); hold on ; 
plot([1,2], n_checks, 'k-o', 'linewidth', 2, 'MarkerSize', 10) ; 


% Save data to a latex table 
fn_res = ['./results_ratings/'] ; 
fn_tab_res = [fn_res, 'Checks_', num2str(n_subj), 'subj.txt'] ;
write_latex_table_latencies(strsplit(num2str(all_subj_idx)), checks_T1_T2, fn_tab_res, ...
    {'Pre-mistakes (\#/15)', 'Post-mistakes (\#/15)', '$I_1$ ($\degrees$C)', '$I_2$ (\degrees C)', 'age'}, ...
    '', '')


end
