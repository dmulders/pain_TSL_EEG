function [] = TSL_plot_IO_fit()
% Reload the data allowing to display the IO fitting on the EEG, for *one*
% given parameter. 
%
% Data must have been computed and saved by running TSL_analyze_EEG.m with
% IO_fit_opt = 1 and the corresponding model parameters. 

warning('off','all') ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
all_subj_idx    = setdiff(1:36,[1, 11, 15, 28, 33]) ; 
LP_filter       = 1 ;       % Low-pass filter the epochs
LP_freq         = 30 ;      % frequency cutoff for the LP filtering
reject_ampl     = 1 ;       % reject epochs based on amplitude
max_ampl_rej    = 80 ;      % amplitude criterion to reject epochs
                            % applied after LP AND after bc...
                            % EEG data of all subjects. 
reject_usual    = 0 ;       % reject usual epochs based on signal with LP30 >80
random_inter    = 1 ;       % random intercepts for the different conditions (TPs).
moving_avg      = 0 ;       % moving avg on each EEG *trial*
movmu_ms        = 80 ;      % msec for the moving avg
p_predictor     = 0 ;       % include proba as predictor
rconf_no_logN   = 0 ;       % Add log(N) as a regressor to compute residual conf
conf_prev       = 1 ;       % use confidence in t-1 as predictor in t!
CB_PT           = 0 ;       % cluster-based permutation testing to correct for MC
n_perm          = 2000 ;    % nb of permutations for CB-PT
cf_alpha        = 0.05 ;    % cluster-forming p-value 
% (1) cluster-forming p-val: alpha_level in compute_cluster_perm_stats.m & compute_cluster_pv
% (2) cluster p: 
% - in TSL_analyze_EEG.m, cluster_thershold in compute_criticals
% - in TSL_plot_IO_fit: in all_h_beta=all_p_beta<=alpha_level ;
reject_first    = 0 ;       % reject the first epochs
reject_perI     = 1 ;       % for I1 and I2 separately
n_first_rej     = 2 ;      % # of initial epochs to reject
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
select_subj     = 0 ;       % only show coeff from one subject
subj_idx        = 1 ;       % 2, 23, 27(!), 30
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FDR_correct     = 1 ;       % for the t-stats (if ~CB_PT)
show_R2 = 0 ;               % explained var plots
show_avg = 0 ;              % regression coefficients averaged over all the channels
show_topos      = 0 ;       % 
% Topoplot options 
% times at which to show topo of t-stats, in ms -- for BPE then rconf
t_topo_I1 ={[200, 720, 820],[220, 330, 564]} ; % check topo I1 at peaks of t-stats
t_topo_I1 = {[],[220, 308]} ; % -- TP/IF leak 8 


% N1-N2: 190-200ms
% P2: 330
t_topo_I2 = {[432, 708, 824, 908], [376, 512]} ; % at peaks of t stats
t_topo_I2 = {[], [388, 500]} ; % -- TP leak 8
%t_topo_I2 = {[736, 936], [388, 500]} ; % -- IF leak 8
%t_topo_I2 = {[736, 936], [180, 208, 300, 388, 500]} ; % -- indiv optimal models

% 420ms: min of GFP --> topo meaningless
% N1-N2: 355
% P2: 508
show_GFP = 0 ; 

chan_selected       = {'FCZ'} ; %{'FCZ','C3', 'C4', 'CPz'} ; 
chan_selected       = {'FCZ', 'C3', 'C4', 'CPz', 'Cz'} ;

topo_colors         = 1 ;   % show color codes for the avg signals on topo


% Define channels to show on plot
show_signif_chans   = 1 ;   % show channels where there are significant effects
selected_chan       = 0 ;   % only plot chan_selected on the butterfly plot
only_keep_topo_chans = 0 ;  % don't show channels out of topoplot in butterfly plot
% Default to reproduce figs in article
if CB_PT
    show_signif_chans = 1 ;selected_chan = 0 ; only_keep_topo_chans = 0 ;  
elseif FDR_correct
    selected_chan = 1 ;
end
    

% for significant clusters, connected components are computed using all 64 chans
if show_signif_chans ; only_keep_topo_chans = 0 ; end
if selected_chan ; show_signif_chans = 0 ; end

% IO model -------------------------------------------------------------- %
use_indiv_mod   = 0 ;   % optimal individual models
use_indiv_leak  = 0 ;   % optimal individual leaks
MemParam = {'Decay', [8]} ; % 6
AboutFirst      = 'WithoutFirst' ;
learned_param   = 'transition' ; % 'frequency', 'transition', 'alternation'

% General plot params --------------------------------------------------- %
eps_fig = 0 ; fig_fig = 0 ; pdf_fig = 0 ; % otherwise: png
% Specific plot params -------------------------------------------------- %
plot_ms = 1 ;               % msec time unit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[window_int, windowF_int, window, windowF, leaky_int, leakyF_int, ...
    decay, decayF, weight, ~] = get_memParams(MemParam) ; 

str_model = get_model_param_str(AboutFirst, learned_param, decay, decayF, ...
    window, windowF, leaky_int, leakyF_int, window_int, windowF_int, weight) ;
if use_indiv_mod
    str_model = ['_opt', str_model(4:end)] ; 
    if use_indiv_leak        
        idx_leak = strfind(str_model, 'leak') ;
        str_model = [str_model(1:(idx_leak-1)), 'opt'] ; 
    end        
end


fn_res = ['./results_EEG/'] ; 

if reject_first
    if reject_perI
        fn_res = [fn_res(1:end-1), '_wo', num2str(n_first_rej), 'I/'] ; 
    else
        fn_res = [fn_res(1:end-1), '_wo', num2str(n_first_rej), '/'] ; 
    end
end

n_subj = length(all_subj_idx) ; 
if n_subj>1
    str_subj = ['_',num2str(n_subj),'subj'] ; 
else
    str_subj = ['_subj',num2str(num2str_padd(all_subj_idx(1), 2))] ;
end
str_LP = '' ; 
if LP_filter
    str_LP = ['_LP', num2str(LP_freq)] ; 
end
str_ampl_max = '' ; 
if reject_ampl
    str_ampl_max = ['_maxA', num2str(max_ampl_rej)] ; 
end
str_movmu = '' ; 
if moving_avg
    str_movmu = ['_mmu', num2str(movmu_ms)] ;
end
str_params = [str_subj, str_LP, str_movmu] ;

n_chans_plot = length(chan_selected) ;
str_plot = '' ; 
if show_signif_chans
    str_plot = ['_signif_chans'] ; 
elseif selected_chan
    str_plot = ['_',num2str(n_chans_plot),'chans'] ; 
elseif only_keep_topo_chans
    str_plot = '_topo_chans' ; 
end


str_inter = '' ; 
if random_inter
   str_inter = '_rinter' ;  
end
str_p = '' ; 
if p_predictor ; str_p = '_pp' ; end
if rconf_no_logN ; str_p = [str_p, '_noN'] ; end % noN = no Nb of observations
if conf_prev ; str_p = [str_p, '_pc'] ; end % pc = previous confidence used
if CB_PT 
    str_p = [str_p, '_',num2str(n_perm),'CP'] ; 
    if cf_alpha~=0.05
        str_p = [str_p, '_', num2str(round(1000*cf_alpha))] ; 
    end
end

str_final = [str_model, str_inter,str_p] ; 


fn_data_tmp = [fn_res, 'data/'] ; 
fn_IO_fit = [fn_data_tmp, 'IO_fit', str_params, str_ampl_max, str_final, ...
    '.mat'] 

if select_subj
    str_params = ['_subj', num2str(subj_idx),str_LP, str_movmu] ; 
end

gray_col = [200,200,200]./255 ;
signif_col = [182, 0, 0]./255 ;% 210 
% ----- Topoplot colormap parameters
avg_cmap_levels = 0 ; 
numC_topo = 8 ;                 % if avg topoplots use filled levels 
topo_style = 'both' ;           % 'both', 'map', cfr test_scalp_plot
n_colors_def = 64 ;             % default nb of colors for the colormaps
if avg_cmap_levels
    cmap_topo = parula(numC_topo) ; %parula ; %jet ; %gnuplot2 ; %jet ; %plasma ; % viridis
else
    cmap_topo = parula(n_colors_def) ; %jet(n_colors_def) ; %parula ; 
    cmap_topo = jet(n_colors_def) ; % LVA
    cmap_topo = my_seismic(n_colors_def) ; 
    % perceptually uniform & colorblind friendly cmaps: viridis, inferno, plasma, magma
    % diverging: my_PiYG, my_bwr, my_seismic, my_coolwarm
    % best: seismic (white around 0) & coolwarm (no white, from blue to red)
end
add_hc_topo = 1 ; % colorbar for topoplots


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saved_data = load(fn_IO_fit) ;
indices_TP = saved_data.indices_TP ; 
time_vec = saved_data.time_vec;
fs = saved_data.fs;
xstep = saved_data.xstep;
chanlocs = saved_data.chanlocs;
y_I1 = saved_data.y_I1; % (n_chan, n_times, n_subj) 
y_I2 = saved_data.y_I2 ; 
time_vec_IO = saved_data.time_vec_IO ;

% compute GFP
GFP_I1 = squeeze(compute_GFP(y_I1)) ; % (n_times, n_subj)
GFP_I2 = squeeze(compute_GFP(y_I2)) ; 
mean_GFP_I1 = 2.*(squeeze(mean(GFP_I1,2)))' ;  mean_GFP_I2 = 2.*(squeeze(mean(GFP_I2,2)))' ; 
mean_GFP_I1 = 8.*(mean_GFP_I1-min(mean_GFP_I1(:)))./(max(mean_GFP_I1(:))-min(mean_GFP_I1(:))) ; 
mean_GFP_I2 = 8.*(mean_GFP_I2-min(mean_GFP_I2(:)))./(max(mean_GFP_I2(:))- min(mean_GFP_I2(:))) ; 
std_GFP_I1 = (squeeze(std(GFP_I1,[],2)))' ; std_GFP_I2 = (squeeze(std(GFP_I2,[],2)))' ; 

all_beta = squeeze(saved_data.all_beta) ; % (4, n_chan, n_times_IO, n_subj, 2)
all_R_squared = squeeze(saved_data.all_R_squared) ; % (2, n_chan, n_times_IO, n_subj, 2) 
all_MSE = squeeze(saved_data.all_MSE) ; % (2, n_chan, n_times_IO, n_subj, 2) 

all_names = {chanlocs.labels} ; 
mask_chan_topo = logical([chanlocs.topo_enabled]) ; 
chan_topo = find(mask_chan_topo) ;


% 1 color per channel that can be plotted on topo
if only_keep_topo_chans
    n_curves = length(chan_topo) ;
elseif selected_chan
    n_curves = n_chans_plot ;
else
    n_curves = length(chanlocs) ;
end
if n_curves==1
    cmap = [0, 128, 255]./255 ; 
    cmap = viridis(5+1) ;
    cmap = cmap(2,:) ; 
else
    cmap = parula(n_curves) ; %jet(n_curves) ;   % 
    cmap = viridis(n_curves+1) ; %jet(n_curves) ;   % 
    % perceptually uniform & colorblind friendly cmaps: 
    % viridis, inferno, plasma, magma
    cmap = cmap(2:end,:) ; 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exclude channels that cannot be plotted on topo for all future analyses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% need to be done BEFORE calling compute_CC... otherwise could have some 
% time points induced by chans not on topo
if only_keep_topo_chans
    [y_I1, y_I2, all_beta, all_R_squared, all_MSE, chanlocs, all_names] = ...
        select_chans(chan_topo, y_I1, y_I2, all_beta, all_R_squared, all_MSE, chanlocs, all_names) ; 
end
auto_order = 1 ; % auto-order channels on butterfly plot 
if show_signif_chans
    % define the chan order here (to keep it across connected components)
    chan_order = fliplr(define_chan_order_colors(chanlocs)) ;
    [cmap, y_I1, y_I2, all_beta, all_R_squared, all_MSE, chanlocs, ...
        all_names] = change_chan_order_topo(chan_order, cmap, y_I1, y_I2, all_beta, all_R_squared, ...
        all_MSE, chanlocs, all_names) ; 
    
    auto_order = 0 ; 
end
chanlocs_butterfly = chanlocs ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[~, order_selected] = define_chan_order_colors(chanlocs, chan_selected) ; 
% in sec
time_vec_used = time_vec_IO ; 
time_IO_ms = 1000.*time_vec_used ; 
time_vec_GFP = time_vec ; 
xlim_butterfly = [-0.25, 1] ; %[0,1] ; %[-0.25,0.75]  ; %[time_vec_used(1), time_vec_used(end)] ; %[-0.25,0.75] ; % 
xlab = 'time after stimulus (s)' ; 
if plot_ms
    % in seconds
    xlim_butterfly = 1000.*xlim_butterfly ;
    time_vec_used = time_IO_ms ; 
    time_vec_GFP = 1000.*time_vec_GFP ; 
    xlab = 'time after stimulus (ms)' ; 
end
show_sem = 1 ; 
% for plot_mean_std
x_double = [time_vec_used; time_vec_used] ;   
col_mean = [0, 128, 255; 102, 204, 0]./255 ; col_std = col_mean ; 

fn_res_IO = [fn_res, 'IO_fit/'] ; 
if ~exist(fn_res_IO, 'dir') ;  mkdir(fn_res_IO) ; end


fn_IO_param = [fn_res_IO, str_final, '/'] ; % create one folder for each IO param
if ~exist(fn_IO_param, 'dir') ;  mkdir(fn_IO_param) ; end

fn_Rs = [fn_IO_param, 'Rs', str_params] ; %str_ampl_max] ; 
if topo_colors ; fn_Rs = [fn_Rs, '_topo'] ; end

fn_Rs_chan = [fn_Rs, str_plot] ; 
fn_Rs_avg = [fn_Rs, '_avg_chan'] ; % avg on channels

if CB_PT ; alpha_level = 0.02 ; else ; alpha_level = 0.05 ; end

paired_tests = 0 ; mean_H0 = 0 ; always_t_test = 1 ;
two_sided = 1 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% R squared btw (r)conf & BPE (Bayesian prediction error)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VIF: = 1/(1-Rs), to check collinearities btw (r)conf & BPE.
% https://www.geeksforgeeks.org/detecting-multicollinearity-with-vif-python/
if isfield(saved_data, 'all_Rs_VIF')
   %all_Rs_VIF: (2, n_chan, n_subj, 2) ; % btw (r)conf & PE    
   % values are indep of channel, can take first channel
   all_Rs_VIF = squeeze(saved_data.all_Rs_VIF) ;   
   all_Rs_VIF = squeeze(all_Rs_VIF(:, 1, :,:)) ; % (2, n_subj, 2) % [rconf, conf; I1, I2]
   
   all_pvs = ones(2,2) ; 
   for iI=1:2
       curr_VIF = (1./(1-squeeze(all_Rs_VIF(:,:,iI))))' ; % (n_subj, 2)
       [p_values, h_mask, ~, ~, ~, t_stats] = ...
            compute_significance_matrix(curr_VIF, alpha_level, 10, ...
            paired_tests, 5, always_t_test, two_sided) ;
        all_pvs(:,iI) = p_values ; 
   end
   % Avg subjects
   all_Rs_VIF = squeeze(mean(all_Rs_VIF,2)) ; 
   % Variance Inflation Factor (VIF)
   all_VIF = 1./(1-all_Rs_VIF) % rows: rconf, conf; col: I1, I2
   all_Rs_VIF = 100.*all_Rs_VIF % % explained var btw conf & PE
   all_pvs = all_pvs % VIF signif diff from 5?
   
   % VIF should be smaller than 5 to avoid strong collinearities
   % https://en.wikipedia.org/wiki/Variance_inflation_factor#cite_note-Sheather_2009_p.-6
   % Sheather, Simon (2009). A modern approach to regression with R.
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% R squared (Proportion of explained variance by the regression)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfigs for I1 and I2 on same fig (for each regression)

fig_sz = [1 2 32 12] ;

taille = 20 ; 
n_Rs = size(all_R_squared,1) ; 
all_R_names = {'PE and R. Confidence', 'PE and confidence'} ; % for the titles
all_R_acro = {'_PE_rconf', '_PE_conf'} ; % for the filenames
ylab = 'explained variance (%)' ; 

if show_R2
    for idx_R=1:n_Rs
        % all_R_squared: (2, n_chan, n_times_IO, n_subj, 2) 
        Rs_tmp = 100.*squeeze(all_R_squared(idx_R,:,:,:,:)) ;
        if select_subj
            Rs_tmp = Rs_tmp(:,:,subj_idx,:) ; 
        end

        % Avg across subjects
        Rs_avg_subj = squeeze(mean(Rs_tmp,3)) ; 
        % (n_chan, n_times_IO, 2)
        if selected_chan
            Rs_I1 = squeeze(Rs_avg_subj(order_selected, :,1)) ; 
            Rs_I2 = squeeze(Rs_avg_subj(order_selected, :,2)) ; 
            chanlocs_butterfly = chanlocs(order_selected) ; 
        end

        % --*-- ----> I_1
        fig_tmp = figure('units','centimeters','outerposition',fig_sz,...
        'Name', all_R_names{idx_R}) ; 
        h_tmp(1) = subplot(1, 2, 1) ; hold on ; 
        plot_butterfly(time_vec_used, Rs_I1, 'xlim_val', ...
            xlim_butterfly, 'x_lab', xlab, 'topo_colors', ...
            topo_colors, 'chanlocs', chanlocs_butterfly, 'eps_fig', eps_fig, ...
            'fig_fig', fig_fig, 'pdf_fig', pdf_fig, 'create_fig', 0, 'save_fig', 0, ...
            'y_lab', ylab, 'cmap', cmap) ;
        title('$I_1$', 'FontSize', taille,'Interpreter','Latex') ;    
        % --*-- ----> I_2
        h_tmp(2) = subplot(1, 2, 2) ; hold on ; 
        plot_butterfly(time_vec_used, Rs_I2, 'xlim_val', ...
            xlim_butterfly, 'x_lab', xlab, 'topo_colors', ...
            topo_colors, 'chanlocs', chanlocs_butterfly, 'eps_fig', eps_fig, ...
            'fig_fig', fig_fig, 'pdf_fig', pdf_fig, 'create_fig', 0, 'save_fig', 0, ...
            'y_lab', ylab, 'cmap', cmap) ; 
        title('$I_2$', 'FontSize', taille,'Interpreter','Latex') ; 

        linkaxes(h_tmp,'y') ;    

        fn_Rs_tmp = [fn_Rs_chan, all_R_acro{idx_R}] ; 
        my_save_fig(fig_tmp, fn_Rs_tmp, eps_fig, fig_fig, pdf_fig) ; 

        % --*-- ----> Avg across all (topo) chans (+ sem across subjects)
        % Avg across topo chans
        Rs_avg_chan = squeeze(mean(Rs_tmp(chan_topo,:,:,:),1)) ; 
        % (n_times_IO, n_subj, 2)
        if select_subj
            mean_Rs = Rs_avg_chan' ; 
            std_Rs = zeros(size(mean_Rs)) ; 
        else
            mean_Rs = (squeeze(mean(Rs_avg_chan,2)))' ; 
            std_Rs = (squeeze(std(Rs_avg_chan,[],2)))' ;
        end

        if show_sem ; std_Rs = std_Rs./sqrt(n_subj) ; end
        caption_lgd = {'I_1', 'I_2'} ;
        fn_Rs_tmp = [fn_Rs_avg, all_R_acro{idx_R}] ; 
        plot_mean_std(x_double, mean_Rs, std_Rs, 'col_mean', ...
            col_mean, 'x_lab', xlab, 'y_lab', ylab,'fn_save',...
            fn_Rs_tmp, 'fig_fig', fig_fig, 'pdf_fig',pdf_fig, 'eps_fig', eps_fig, ...
            'add_prior', 0, 'shaded_std', 1, 'show_max', 0, 'col_std', col_std, ...
            'caption_lgd', caption_lgd, 'logx_curve', 0, 'xlim_val', xlim_butterfly) ;
    end
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Betas (Fixed effect coefficients).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if show_signif_chans
    fig_sz = [1, 2, 40, 20]; 
end
show_cluster_signif = 1 ; 

[n_betas, n_chan, n_times_IO, ~,~] = size(all_beta) ; 
% all_beta: (n_beta, n_chan, n_times_IO, n_subj, 2) 


t_crit = tinv(1-alpha_level/2, n_subj-1) ; 
if isfield(saved_data, 'all_t_beta') 
    all_t_beta = squeeze(saved_data.all_t_beta) ; % (n_betas, n_chan, nT_TF, 2)
    all_p_beta = squeeze(saved_data.all_p_beta) ; 
    if CB_PT
        disp(['--- CB corrected p-values'])
       % p_beta_CP, p_beta_CP_max: (n_betas, n_chan, nT_data, 2) ; 
       %all_p_beta = squeeze(saved_data.p_beta_CP) ; % saved_data.p_beta_CP_max ; 
       all_p_beta = squeeze(saved_data.p_beta_CP_max) ; 
    end
    if only_keep_topo_chans
        all_p_beta = all_p_beta(:,chan_topo,:,:) ; 
        all_t_beta = all_t_beta(:,chan_topo,:,:) ; 
    elseif show_signif_chans
        all_p_beta = all_p_beta(:,chan_order,:,:) ; 
        all_t_beta = all_t_beta(:,chan_order,:,:) ; 
    end
else
    % ======== Compute significance of beta coefficients
    % ======================================================================= %    
    all_t_beta = zeros(n_betas, n_chan, n_times_IO, 2) ; 
    all_p_beta = ones(n_betas, n_chan, n_times_IO, 2) ; 
    for i_beta = 1:n_betas
        for i_I=1:2
            curr_betas = squeeze(all_beta(i_beta, :, :, :,i_I)) ;
            curr_betas = permute(curr_betas, [3, 1, 2]) ; % (n_subj, n_chan, n_times_IO)
            curr_betas = reshape(curr_betas, n_subj, n_chan*n_times_IO) ; 
            [p_values, h_mask, ~, ~, ~, t_stats] = ...
                compute_significance_matrix(curr_betas, alpha_level, 10, ...
                paired_tests, mean_H0, always_t_test, two_sided) ;

            t_stats = reshape(t_stats, n_chan, n_times_IO) ; 
            all_t_beta(i_beta, :, :, i_I) = t_stats ; 
            p_values = reshape(p_values, n_chan, n_times_IO) ; 
            all_p_beta(i_beta,:,:,i_I) = p_values ; 
        end
    end
end

% Subfigs for I1 and I2 on same fig (for each regression)
fn_B = [fn_IO_param, 'Beta', str_params] ; %str_ampl_max] ; 
fn_t = [fn_IO_param, 'T', str_params] ; 
fn_Rs_IV = [fn_IO_param, 'Rs_IV', str_params] ; 
if topo_colors ; fn_B = [fn_B, '_topo'] ; end
fn_B_chan = [fn_B, str_plot] ; fn_t_chan = [fn_t, str_plot] ; fn_Rs_chan = [fn_Rs_IV, str_plot] ; 
fn_B_avg = [fn_B, '_avg_chan'] ; % avg on channels
fn_t_avg = [fn_t, '_avg_chan'] ;

if p_predictor
    all_b_names = {'PE (R. Conf. regress)','R. Confidence', 'Proba. estimates',...
        'PE', 'Confidence'} ; % for the titles
    all_b_acro = {'_PE_src', '_rconf', '_p', '_PE', '_conf'} ; % for the filenames
    % src = surprise + residual conf
else
    all_b_names = {'PE (R. Conf. regress)','R. Confidence', ...
        'PE', 'Confidence'} ; % for the titles
    all_b_acro = {'_PE_src', '_rconf', '_PE', '_conf'} ; % for the filenames
end
ylab = 'regression coefficient' ; 
ylab_t = 't-value' ; 


chan_used_I1 = 1:n_chan ; %chan_topo ; 
chan_used_I2 = 1:n_chan ; %chan_topo ; 
if ~(show_signif_chans && CB_PT)    
    if selected_chan
        chan_used_I1 = order_selected ;  
        chan_used_I2 = order_selected ;        
    end
    chanlocs_I1 = chanlocs(chan_used_I1) ; 
    chanlocs_I2 = chanlocs(chan_used_I2) ; 
    chan_idx_I1 = {chan_used_I1} ; n_CC_I1 = 1 ; 
    chan_idx_I2 = {chan_used_I2} ; n_CC_I2 = 1 ; 
end
% ======== FDR correction for MC
% ======================================================================= % 
if FDR_correct && ~CB_PT    
    all_h_beta = zeros(size(all_p_beta)) ; % (n_betas, n_chan, n_times_IO, 2)
    for idx_b=1:n_betas
        % =================== I1
        p_tmp = squeeze(all_p_beta(idx_b,chan_used_I1,:,1)) ; 
        h_tmp =  significance_FDR_corrected(p_tmp, alpha_level) ;
        all_h_beta(idx_b,chan_used_I1,:,1) = h_tmp ;
        
        % =================== I2
        p_tmp = squeeze(all_p_beta(idx_b,chan_used_I2,:,2)) ; 
        h_tmp =  significance_FDR_corrected(p_tmp, alpha_level) ; 
        all_h_beta(idx_b,chan_used_I2,:,2) = h_tmp ;
    end    
else
    if CB_PT
        all_h_beta = all_p_beta<=alpha_level/2 ; % one-sided test 
        %disp('------ using one-sided signif threshold')
    else
        all_h_beta = all_p_beta<=alpha_level ; 
    end
end


idx_b_topos = 1:4 ; %1:4 ; %3:4; 
b_ylim = [-0.4, 0.4;-0.7, 0.5] ; % ylim for beta coeff -- BPE and conf
t_ylim = [-8.5, 8.5; -8.5, 8.5] ; % ylim for tstats -- BPE and conf
t_ylim = [-9, 9; -9, 9] ; % IF8

if CB_PT
   % only keep significant cluster in [0, 1.5] time window
    ind_zero = msec_to_ind(time_vec_used, [0]) ; 
    all_h_beta(:,:,1:ind_zero,:) = 0 ;  
end

%%%%%%%%%%% 

for idx_b=1:n_betas
      
    if show_signif_chans
        % all_h_beta: (n_betas, n_chan, n_times_IO, 2)
        mask_tmp = squeeze(all_h_beta(idx_b,:,:,:)) ; % :
        
        % only keep significant cluster in [0, 1.5] time window
        ind_zero = msec_to_ind(time_vec_used, [0]) ;
        mask_tmp(:,1:ind_zero,:) = 0 ; 
        
        mask_tmp1 = squeeze(mask_tmp(:,:,1)) ; 
        
        mask_tmp2 = squeeze(mask_tmp(:,:,2)) ;
        % tstats to order clusters according to cluster effect size
        tstats_tmp = squeeze(all_t_beta(idx_b,:,:,1)) ; %all_t_beta: (n_betas, n_chan, nT_TF, 2)
        
        [chan_idx_I1, all_pts_I1] = compute_CC(mask_tmp1, n_times_IO, all_names, tstats_tmp) ;
        [chan_idx_I1, all_pts_I1] = keep_positive_clusters(chan_idx_I1, all_pts_I1,time_vec_used) ;% positive time points
        n_CC_I1 = length(chan_idx_I1) ; 
        nrows_I1 = ceil(n_CC_I1/4) ; ncols_I1 = ceil(n_CC_I1/nrows_I1) ;  

        tstats_tmp = squeeze(all_t_beta(idx_b,:,:,2)) ;
        [chan_idx_I2,all_pts_I2] = compute_CC(mask_tmp2, n_times_IO, all_names, tstats_tmp) ; 
        [chan_idx_I2, all_pts_I2] = keep_positive_clusters(chan_idx_I2, all_pts_I2,time_vec_used) ;
        n_CC_I2 = length(chan_idx_I2) ; 
        nrows_I2 = ceil(n_CC_I2/4) ; ncols_I2 = ceil(n_CC_I2/nrows_I2) ; 

        %%%% MIX all clusters for beta coeffs:
        mask_tmp1 = (sum(mask_tmp(:,:,1),2))' ; % row
        mask_tmp2 = (sum(mask_tmp(:,:,2),2))' ; 
        chan_used_I1 = find(mask_tmp1) ; %find(mask_tmp1 & mask_chan_topo) ; 
        chan_used_I2 = find(mask_tmp2) ;   
        chanlocs_I1 = chanlocs(chan_used_I1) ; 
        chanlocs_I2 = chanlocs(chan_used_I2) ;
        %signif_pts_I1 = find(sum(mask_tmp(:,:,1),1)) ; 
        %signif_pts_I2 = find(sum(mask_tmp(:,:,2),1)) ;
    else%if selected_chan || only_keep_topo_chans
        % all_h_beta: (n_betas, n_chan, n_times_IO, 2)
        mask_tmp = squeeze(all_h_beta(idx_b,:,:,:)) ; 

        all_pts_I1 = {find(sum(mask_tmp(chan_used_I1,:,1),1))} ; 
        all_pts_I2 = {find(sum(mask_tmp(chan_used_I2,:,2),1))} ; 
    end
    
    % all_beta: (n_beta, n_chan, n_times_IO, n_subj, 2)
    B_tmp = shiftdim(all_beta(idx_b,:,:,:,:),1) ; % (n_chan, n_times_IO, n_subj, 2)
    if select_subj
        B_tmp = B_tmp(:,:,subj_idx,:) ; 
    end
    
    B_I1 = squeeze(mean(B_tmp(chan_used_I1,:,:,1),3)) ; % (n_chan, n_times_IO)
    B_I2 = squeeze(mean(B_tmp(chan_used_I2,:,:,2),3)) ; 

    
    % --*-- ----> I_1
    if selected_chan || only_keep_topo_chans ; cmap_tmp = cmap ; else ; cmap_tmp = cmap(chan_used_I1,:) ; end
    fig_tmp = figure('units','centimeters','outerposition',fig_sz,...
        'Name', all_b_names{idx_b}) ; 
    h_tmp(1) = subplot(1, 2, 1) ; hold on ; 
    
    if ~isempty(chan_used_I1)
    plot_butterfly(time_vec_used, B_I1, 'xlim_val', ...
        xlim_butterfly, 'x_lab', xlab, 'topo_colors', ...
        topo_colors, 'chanlocs', chanlocs_I1, 'eps_fig', eps_fig, ...
        'fig_fig', fig_fig, 'pdf_fig', pdf_fig, 'create_fig', 0, 'save_fig', 0, ...
        'y_lab', ylab, 'cmap', cmap_tmp, 'order_top_colors', 1, 'auto_order', auto_order) ; ylim(b_ylim(mod(idx_b-1,2)+1,:)) ; 
    end
    %plot([0,0], ylim, 'k--') ; plot(xlim, [0,0], 'k--') ; 
    title('$I_1$', 'FontSize', taille,'Interpreter','Latex') ;  
    % --*-- ----> I_2
    if selected_chan || only_keep_topo_chans ; cmap_tmp = cmap ; else ; cmap_tmp = cmap(chan_used_I2,:) ; end
    h_tmp(2) = subplot(1, 2, 2) ; hold on ; 
    if ~isempty(chan_used_I2)
    plot_butterfly(time_vec_used, B_I2, 'xlim_val', ...
        xlim_butterfly, 'x_lab', xlab, 'topo_colors', ...
        topo_colors, 'chanlocs', chanlocs_I2, 'eps_fig', eps_fig, ...
        'fig_fig', fig_fig, 'pdf_fig', pdf_fig, 'create_fig', 0, 'save_fig', 0, ...
        'y_lab', ylab, 'cmap', cmap_tmp, 'order_top_colors', 1, 'auto_order', auto_order) ; ylim(b_ylim(mod(idx_b-1,2)+1,:)) ; 
    %plot([0,0], ylim, 'k--') ; plot(xlim, [0,0], 'k--') ;
    end
    title('$I_2$', 'FontSize', taille,'Interpreter','Latex') ; 
    linkaxes(h_tmp,'y') ;      
    fn_B_tmp = [fn_B_chan, all_b_acro{idx_b}] ; 
    my_save_fig(fig_tmp, fn_B_tmp, eps_fig, fig_fig, pdf_fig) ; 
    if show_avg
        % --*-- ----> Avg across all (topo) chans (+ sem across subjects)
        % Avg across topo chans
        B_avg_chan = squeeze(mean(B_tmp(chan_topo,:,:,:),1)) ; 
        % (n_times_IO, n_subj, 2)

        if select_subj
            mean_B = B_avg_chan' ; 
            std_B = zeros(size(mean_B)) ; 
        else
            mean_B = (squeeze(mean(B_avg_chan,2)))' ; 
            std_B = (squeeze(std(B_avg_chan,[],2)))' ;
        end

        if show_sem ; std_B = std_B./sqrt(n_subj) ; end
        caption_lgd = {'I_1', 'I_2'} ;
        fn_B_tmp = [fn_B_avg, all_b_acro{idx_b}] ; 
        plot_mean_std(x_double, mean_B, std_B, 'col_mean', ...
            col_mean, 'x_lab', xlab, 'y_lab', ylab,'fn_save',...
            fn_B_tmp, 'fig_fig', fig_fig, 'pdf_fig',pdf_fig, 'eps_fig', eps_fig, ...
            'add_prior', 0, 'shaded_std', 1, 'show_max', 0, 'col_std', col_std, ...
            'caption_lgd', caption_lgd, 'logx_curve', 0, 'xlim_val', xlim_butterfly) ;
    end
    
    % ====================== Idem with t-stats
    if ~select_subj
        % all_t_beta: (n_beta, n_chan, n_times_IO, 2) 
        fig_tmp = figure('units','centimeters','outerposition',fig_sz,...
            'Name', all_b_names{idx_b}) ; 
        for i_cc1=1:n_CC_I1
            % --*-- ----> I_1
            chan_used_tmp = chan_idx_I1{i_cc1} ;   
            nc_tmp = length(chan_used_tmp) ; 
                        
            t_I1 = reshape(all_t_beta(idx_b, chan_used_tmp, :,1), nc_tmp, n_times_IO) ; 
            chanlocs_tmp = chanlocs(chan_used_tmp) ; 
            
            if show_signif_chans
                h_t(i_cc1) = subplot(nrows_I1, ncols_I1, i_cc1) ; hold on ; 
            else
                h_t(1) = subplot(1, 2, 1) ; hold on ; 
            end
            if show_GFP
                % ===== SHOW GFP
                fill([time_vec_GFP, fliplr(time_vec_GFP)], ...
                    [mean_GFP_I1-std_GFP_I1, fliplr(mean_GFP_I1+std_GFP_I1)], gray_col, ...
                    'EdgeColor','none', 'FaceAlpha', 0.2) ; hold on ;      
                plot(time_vec_GFP,mean_GFP_I1,'color', gray_col, 'linewidth', 2) ;    
            end
            hold on ;
            if selected_chan || only_keep_topo_chans ; cmap_tmp = cmap ; else ; cmap_tmp = cmap(chan_used_tmp,:) ; end
            plot_butterfly(time_vec_used, t_I1, 'xlim_val', ...
                xlim_butterfly, 'x_lab', xlab, 'topo_colors', ...
                topo_colors, 'chanlocs', chanlocs_tmp, 'eps_fig', eps_fig, ...
                'fig_fig', fig_fig, 'pdf_fig', pdf_fig, 'create_fig', 0, 'save_fig', 0, ...
                'y_lab', ylab_t, 'add_HL', 1, 'HL_val', t_crit, 'neg_pos_HL', 1, ...
                'cmap', cmap_tmp, 'order_top_colors', 1, 'auto_order', auto_order) ; ylim(t_ylim(mod(idx_b-1,2)+1,:)) ; 
            %plot([0,0], ylim, 'k--') ; plot(xlim, [0,0], 'k--') ; 
            title('$I_1$', 'FontSize', taille,'Interpreter','Latex') ;
            if show_cluster_signif
                % signif_pts_I1 = find(sum(mask_tmp(:,:,1),1)) ; 
                signif_pts_I1 = all_pts_I1{i_cc1} ;    
                yL = get(gca, 'YLim') ; 
                npts = length(signif_pts_I1) ; 
                hold on ; 
                plot(time_vec_used(signif_pts_I1), ones(1,npts).*yL(1), 's', ...
                    'color', signif_col, 'MarkerSize', 1, 'Linewidth',6) ; 
            end
            if show_topos && any(idx_b==idx_b_topos)
                hold on ;     
                % Indicate times at which show topo
                yL = get(gca, 'ylim') ; 
                for t_I1=t_topo_I1{mod(idx_b-1,2)+1} 
                    [~,idx_topo] = min(abs(time_vec_GFP-t_I1));  
                    plot([t_I1, t_I1],[yL(1),mean_GFP_I1(idx_topo)],':','color', 'k', 'linewidth', 1) ; 
                    hold on ; 
                end
            end
        end
        
        
        if show_signif_chans && exist('h_t', 'var')
            linkaxes(h_t,'y') ;      
            fn_t_tmp = [fn_t_chan, all_b_acro{idx_b}, '_I1'] ; 
            my_save_fig(fig_tmp, fn_t_tmp, eps_fig, fig_fig, pdf_fig) ; 
        
            fig_tmp = figure('units','centimeters','outerposition',fig_sz,...
                'Name', all_b_names{idx_b}) ; 
        end
        
        for i_cc2=1:n_CC_I2
            % --*-- ----> I_2            
            chan_used_tmp = chan_idx_I2{i_cc2} ;      
            nc_tmp = length(chan_used_tmp) ; 
            t_I2 = reshape(all_t_beta(idx_b, chan_used_tmp, :,2), nc_tmp, n_times_IO) ;            
            chanlocs_tmp = chanlocs(chan_used_tmp) ;            
            
            if show_signif_chans
                h_t(i_cc2) = subplot(nrows_I2, ncols_I2, i_cc2) ; hold on ; 
            else
                h_t(2) = subplot(1, 2, 2) ; hold on ; 
            end
            if show_GFP
                % ===== SHOW GFP
                fill([time_vec_GFP, fliplr(time_vec_GFP)], ...
                    [mean_GFP_I2-std_GFP_I2, fliplr(mean_GFP_I2+std_GFP_I2)], gray_col, ...
                    'EdgeColor','none', 'FaceAlpha', 0.2) ; hold on ;      
                plot(time_vec_GFP,mean_GFP_I2,'color', gray_col, 'linewidth', 2) ;    
            end
            hold on ;  
            if selected_chan || only_keep_topo_chans ; cmap_tmp = cmap ; else ; cmap_tmp = cmap(chan_used_tmp,:) ; end
            plot_butterfly(time_vec_used, t_I2, 'xlim_val', ...
                xlim_butterfly, 'x_lab', xlab, 'topo_colors', ...
                topo_colors, 'chanlocs', chanlocs_tmp, 'eps_fig', eps_fig, ...
                'fig_fig', fig_fig, 'pdf_fig', pdf_fig, 'create_fig', 0, 'save_fig', 0, ...
                'y_lab', ylab_t,  'add_HL', 1, 'HL_val', t_crit, 'neg_pos_HL', 1, ...
                'cmap', cmap_tmp, 'order_top_colors', 1, 'auto_order', auto_order) ; ylim(t_ylim(mod(idx_b-1,2)+1,:)) ; 
            %plot([0,0], ylim, 'k--') ; plot(xlim, [0,0], 'k--') ; 
            title('$I_2$', 'FontSize', taille,'Interpreter','Latex') ; 
            if show_cluster_signif
                signif_pts_I2 = all_pts_I2{i_cc2} ; 
                yL = get(gca, 'YLim') ; 
                npts = length(signif_pts_I2) ; 
                hold on ; 
                plot(time_vec_used(signif_pts_I2), ones(1,npts).*yL(1), 's', ...
                    'color', signif_col, 'MarkerSize', 1, 'Linewidth',6) ; 
            end
            if show_topos && any(idx_b==idx_b_topos)
                hold on ;     
                % Indicate times at which show topo
                yL = get(gca, 'ylim') ; 
                for t_I2=t_topo_I2{mod(idx_b-1,2)+1} 
                    [~,idx_topo] = min(abs(time_vec_GFP-t_I2));  
                    plot([t_I2, t_I2],[yL(1),mean_GFP_I2(idx_topo)],':','color', 'k', 'linewidth', 1) ; 
                    hold on ; 
                end
            end
        end
        if exist('h_t', 'var')
        linkaxes(h_t,'y') ;   
        end
        if show_signif_chans
            fn_t_tmp = [fn_t_chan, all_b_acro{idx_b}, '_I2'] ; 
        else
            fn_t_tmp = [fn_t_chan, all_b_acro{idx_b}] ; 
        end
        my_save_fig(fig_tmp, fn_t_tmp, eps_fig, fig_fig, pdf_fig) ;  
        if CB_PT 
            disp('Press enter to close current figs and continue running the function...') ; 
            pause ; close all ; 
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Topoplots at selected times
        if show_topos && any(idx_b==idx_b_topos)
            for t_I1=t_topo_I1{mod(idx_b-1,2)+1}
                % Extract topo at time t_I1
                [~,idx_topo] = min(abs(time_IO_ms-t_I1));  

                % 1. standardize topo of t stats
                % all_t_beta: (n_beta, n_chan, n_times_IO, 2) 
                y_t = squeeze(all_t_beta(idx_b,:,idx_topo,1)) ; % (n_chan, 1) 
                y_t = zscore(y_t) ; %(y_t - mean(y_t))./( sqrt(sum(y_t.^2)) ) ; 

                fn_tmp = [fn_IO_param, 'Topo_I1_', num2str(t_I1), 'ms', str_params, all_b_acro{idx_b}] ;
                plot_single_topo(y_t, chanlocs, fn_tmp, ['I_1 - t = ', num2str(t_I1), ' ms'], ...
                    avg_cmap_levels, eps_fig, fig_fig, pdf_fig, topo_style, numC_topo, cmap_topo, add_hc_topo) ; 
            end
            % Idem I2
            for t_I2=t_topo_I2{mod(idx_b-1,2)+1}

                % Extract topo at time t_I1
                [~,idx_topo] = min(abs(time_IO_ms-t_I2));  

                % 1. standardize topo of each subject
                y_t = squeeze(all_t_beta(idx_b,:,idx_topo,2)) ; 
                y_t = zscore(y_t) ; 

                fn_tmp = [fn_IO_param, 'Topo_I2_', num2str(t_I2), 'ms', str_params,all_b_acro{idx_b}] ;
                plot_single_topo(y_t, chanlocs, fn_tmp, ['I_2 - t = ', num2str(t_I2), ' ms'], ...
                    avg_cmap_levels, eps_fig, fig_fig, pdf_fig, topo_style, numC_topo, cmap_topo, add_hc_topo) ; 
            end
        end
        
    end    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% R2 for each IV separately
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
show_indiv_R2 = 0 ; 
if isfield(saved_data, 'all_Rs_IV') && show_indiv_R2
    all_Rs_IV = 100.*squeeze(saved_data.all_Rs_IV) ; 
    % all_Rs_IV: (n_betas, n_chan, nT_data, n_subj, 2)
    ylab = 'explained variance (%)' ; 
    
    mean_GFP_I1 = 1.*(mean_GFP_I1-min(mean_GFP_I1(:)))./(max(mean_GFP_I1(:))-min(mean_GFP_I1(:))) ;
    mean_GFP_I2 = 1.*(mean_GFP_I2-min(mean_GFP_I2(:)))./(max(mean_GFP_I2(:))- min(mean_GFP_I2(:))) ;
    
    for idx_b=1:n_betas        
        Rs_tmp = shiftdim(all_Rs_IV(idx_b,:,:,:,:),1) ; % (n_chan, n_times_IO, n_subj, 2)
        if select_subj
            Rs_tmp = Rs_tmp(:,:,subj_idx,:) ; 
        end
        Rs_I1 = squeeze(mean(Rs_tmp(chan_used_I1,:,:,1),3)) ; % (n_chan, n_times_IO)
        Rs_I2 = squeeze(mean(Rs_tmp(chan_used_I2,:,:,2),3)) ;     
    
        % --*-- ----> I_1
        fig_tmp = figure('units','centimeters','outerposition',fig_sz,...
            'Name', all_b_names{idx_b}) ; 
        h_tmp(1) = subplot(1, 2, 1) ; hold on ; 
        if show_GFP
            % ===== SHOW GFP
            fill([time_vec_GFP, fliplr(time_vec_GFP)], ...
                [mean_GFP_I1-std_GFP_I1, fliplr(mean_GFP_I1+std_GFP_I1)], gray_col, ...
                'EdgeColor','none', 'FaceAlpha', 0.2) ; hold on ;      
            plot(time_vec_GFP,mean_GFP_I1,'color', gray_col, 'linewidth', 2) ;    
        end
        hold on ;
        
        plot_butterfly(time_vec_used, Rs_I1, 'xlim_val', ...
            xlim_butterfly, 'x_lab', xlab, 'topo_colors', ...
            topo_colors, 'chanlocs', chanlocs_I1, 'eps_fig', eps_fig, ...
            'fig_fig', fig_fig, 'pdf_fig', pdf_fig, 'create_fig', 0, 'save_fig', 0, ...
            'y_lab', ylab, 'cmap', cmap, 'order_top_colors', 1) ;
        title('$I_1$', 'FontSize', taille,'Interpreter','Latex') ;    
        % --*-- ----> I_2
        h_tmp(2) = subplot(1, 2, 2) ; hold on ; 
        
        if show_GFP
            % ===== SHOW GFP
            fill([time_vec_GFP, fliplr(time_vec_GFP)], ...
                [mean_GFP_I2-std_GFP_I2, fliplr(mean_GFP_I2+std_GFP_I2)], gray_col, ...
                'EdgeColor','none', 'FaceAlpha', 0.2) ; hold on ;      
            plot(time_vec_GFP,mean_GFP_I2,'color', gray_col, 'linewidth', 2) ;     
        end
        hold on ; 
        
        plot_butterfly(time_vec_used, Rs_I2, 'xlim_val', ...
            xlim_butterfly, 'x_lab', xlab, 'topo_colors', ...
            topo_colors, 'chanlocs', chanlocs_I2, 'eps_fig', eps_fig, ...
            'fig_fig', fig_fig, 'pdf_fig', pdf_fig, 'create_fig', 0, 'save_fig', 0, ...
            'y_lab', ylab, 'cmap', cmap, 'order_top_colors', 1) ; 
        title('$I_2$', 'FontSize', taille,'Interpreter','Latex') ; 
        linkaxes(h_tmp,'y') ;      
        fn_Rs_tmp = [fn_Rs_chan, all_b_acro{idx_b}] ; 
        my_save_fig(fig_tmp, fn_Rs_tmp, eps_fig, fig_fig, pdf_fig) ;         
    end
    
    
end



end

function [y_I1, y_I2, all_beta, all_R_squared, all_MSE, chanlocs, all_names] = ...
    select_chans(chan_select, y_I1, y_I2, all_beta, all_R_squared, all_MSE, chanlocs, all_names)
% select channel indices chan_topo in all input files. 

chanlocs = chanlocs(chan_select) ; 
all_names = all_names(chan_select) ; 
y_I1 = y_I1(chan_select,:,:) ; % (n_chan, n_times, n_subj) 
y_I2 = y_I2(chan_select,:,:) ;
all_beta = all_beta(:,chan_select,:,:,:) ; % (4, n_chan, n_times_IO, n_subj, 2)
all_R_squared = all_R_squared(:,chan_select,:,:,:) ; % (2, n_chan, n_times_IO, n_subj, 2)
all_MSE = all_MSE(:,chan_select,:,:,:) ; % (2, n_chan, n_times_IO, n_subj, 2)


end

function [cmap, y_I1, y_I2, all_beta, all_R_squared, all_MSE, chanlocs, ...
    all_names] = change_chan_order_topo(chan_order, cmap, y_I1, y_I2, all_beta, all_R_squared, ...
    all_MSE, chanlocs, all_names) 
% change the channel ordering to be suited for a topoplot. 

cmap = flipud(cmap) ; 
chanlocs = chanlocs(chan_order) ; 
all_names = all_names(chan_order) ; 
y_I1 = y_I1(chan_order,:,:) ; % (n_chan, n_times, n_subj) 
y_I2 = y_I2(chan_order,:,:) ;
all_beta = all_beta(:,chan_order,:,:,:) ; % (4, n_chan, n_times_IO, n_subj, 2)
all_R_squared = all_R_squared(:,chan_order,:,:,:) ; % (2, n_chan, n_times_IO, n_subj, 2)
all_MSE = all_MSE(:,chan_order,:,:,:) ; % (2, n_chan, n_times_IO, n_subj, 2)

end

