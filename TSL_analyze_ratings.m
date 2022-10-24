function [all_corrs, all_MSE, est_p1, IO_p1, IO_MAP_p1, est_conf, IO_conf, ...
    last_I_p1, last_I_conf, TPs_found, p1_poly, MAP_p1_poly, conf_poly] = ...
    TSL_analyze_ratings(varargin)
% Load and analyze the behavioral data of all the subjects, for one model 
% and one parameter set. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ======================================================================= %
% ===== * ===== General parameters
% ======================================================================= %
all_subj_idx    = setdiff(1:36,[1, 11, 15, 28, 33]) ;


str_stim        = '_TCS' ;
n_stim_test     = 100 ; 
fn_dir          = ['./data_ratings/'] ; 
show_plots      = 1 ; 
model_type      = 'IO' ;    % RW, IO
                            % RW = Rescorla-Wagner = delta rule
                            % IO = Ideal Observer
show_sanit_checks = 1 ;     % sanitary checks of number of stimuli received/TPs evaluated
% ======================================================================= %
% ===== * ===== RW + IO: Default parameters
% ======================================================================= %
learned_param   = 'transition' ; % 'frequency', 'transition', 'alternation'

% ======================================================================= %
% ===== * ===== RW: Default parameters
% ======================================================================= %
a_RW = 0.1 ;            % alpha = learning rate in RW model

% ======================================================================= %
% ===== * ===== Ideal Observer: Default parameters
% ======================================================================= %
% Memory parameters
leaky_int = 1 ;         %
leakyF_int = 0 ;        % integrate in flipped direction (more weights at sequence start)
window_int = 0 ;        %
windowF_int = 0 ;       % 

decay = 8 ;             % integration time constant
decayF = 0 ;            % Forward (flipped) integration time constant
weight = 1 ;            % weighting of backward observations compared to forward
window = 30 ;           % window length in backward direction (from last stimulus received)
windowF = 30 ;          % window length in forward direction (from first stimulus in the sequence)
% parameters to test effects of priors on more pain vs. less pain
prior_data = struct('priorN1', 1, 'priorN2', 1, 'factorN1', 1, 'factorN2', 1) ; 
%%%%%%%%%%%%%%%%%%%%%%%%%
MemParam = {} ;
% in.opt.MemParam = {'Decay', 16,'Limited',30};

if leaky_int
    MemParam(end+1:end+2) = {'Decay', decay} ; 
end
if window_int
    MemParam(end+1:end+2) = {'Limited', window} ; 
end
if windowF_int
    MemParam(end+1:end+2) = {'LimitedF', windowF} ;    
end
if leakyF_int
    MemParam(end+1:end+2) = {'DecayF', decayF} ;  
end
if leaky_int && leakyF_int && weight~=1
    MemParam(end+1:end+2) = {'Weight',weight} ;
end
AboutFirst      = 'WithoutFirst' ;



% ======================================================================= %
% ===== * ===== Parse varargin
% ======================================================================= %
nargs_left = length(varargin) ;
if nargs_left > 0
    if ~(round(nargs_left/2) == nargs_left/2)
        error('-- There should be an even number of varargin arguments')
    end
    for i = 1:2:nargs_left
        Name_arg = varargin{i};
        Val_arg = varargin{i+1};
        if ~ischar(Name_arg)
            error('Flag arguments must be strings')
        end
        Name_arg = lower(Name_arg);
        switch Name_arg
            case 'memparam'
                MemParam = Val_arg ;
            case 'aboutfirst'
                AboutFirst = Val_arg ;
            case 'learned_param'
                learned_param = Val_arg ; 
            case 'show_plots'
                show_plots = Val_arg ; 
            case 'all_subj_idx'
                all_subj_idx = Val_arg ; 
            case 'model_type'
                model_type = Val_arg ; 
            case 'a_rw'
                a_RW = Val_arg ; 
        end
    end
end


[window_int, windowF_int, window, windowF, leaky_int, leakyF_int, ...
    decay, decayF, weight, ~] = get_memParams(MemParam) ; 

str_model = get_model_param_str(AboutFirst, learned_param, decay, decayF, ...
    window, windowF, leaky_int, leakyF_int, window_int, windowF_int, weight, ...
    prior_data, model_type, a_RW)  

% ======================================================================= %
fn_res = ['./results_ratings/'] ; 
if strcmpi(model_type,'rw') ; fn_res = './results_ratings_d/' ; end
eps_fig = 0 ; fig_fig = 0 ; pdf_fig = 0 ; % otherwise: png

n_subj = length(all_subj_idx) ; 
str_subj = ['_',num2str(n_subj),'subj'] ; 
TPs_found = {} ;    % TPs found (str_TP)
all_TPs = [] ;      % generative TPs in final order

% Data stored in cells, with one entry for each subject --> vector of all 
% quantities (all TPs merged, otherwise not enough data).  
est_p1 = cell(n_subj,1) ;       % estimated TPs
est_conf = cell(n_subj,1) ;     % estimated confidence

RT_subj = cell(n_subj,1) ;      % reaction times
RT_conf_subj = cell(n_subj,1) ; % recation times for the confidence

IO_p1 = cell(n_subj,1) ;        % IO TPs (mean of posterior)
IO_MAP_p1 = cell(n_subj,1) ;    % IO TPs (MAP of posterior)
IO_conf = cell(n_subj,1) ;      % IO TPs
last_I_p1 = cell(n_subj,1) ;    % last intensity received (when p is estimated)
next_I_p1 = cell(n_subj,1) ;    % next intensity to be received (when p is estimated)
last_I_conf = cell(n_subj,1) ;  % last intensity received (when conf is estimated)

%%%%%%%%%%% Sanitary check data
N_stim = zeros(n_subj,2) ;      % nb of stim of I1 and I2
N_rated = zeros(n_subj,2) ;     % nb of rated I from I1 and from I2

%%%%%%%%%% Data to save
n_out = 4 ; 
n_conds_max = 10 ; 
ratings_IO_data = NaN(n_out, n_stim_test, n_conds_max, n_subj) ;% est_p1, est_conf, IO_p1, IO_conf
fn_ratingsIO = [fn_res, 'data/'] ; if ~exist(fn_ratingsIO, 'dir') ; mkdir(fn_ratingsIO) ; end
fn_ratingsIO = [fn_ratingsIO, 'Ratings', str_model, '.mat'] ; 
remove_missed_resp = 0 ; 


% ======================================================================= %
% Parameters %%%%%%%%%%%%%%%%%%%%%%%%
in.learned = learned_param ;        % estimate TPs
in.jump = 0 ;                       % without jump
%in.mode = 'HMM' ;% use the HMM algo (cfr Meyniel2016 p17) (not sampling, Meyniel2015)
in.opt.MemParam = MemParam ;        % memory limit
in.opt.AboutFirst = AboutFirst ;    % use p(y1|theta) = 0.5 instead of a fct 
                                    % of theta to have analytical solutions
n = 40 ;                            % resolution of the univariate probability grid
in.opt.pgrid = linspace(0,1,n) ;    % grid to return full distributions
in.opt.ReturnDist = 1 ;             % Return full posterior distributions
in.opt.ReturnDistUpdate = 1 ; 
eps_prior = 0.1 ; %0.01 ; 
in.opt.priorp1g2 = [prior_data.priorN1 prior_data.priorN2] + eps_prior ;
                                    % uniform Beta prior (when learning TP)
                                    % Add eps to 1 to the prior: avoid NaN
                                    % in MAP when learning TP and AF after
                                    % the first observation (cfr MAP in
                                    % mar_ComputeMAPandPrediction and
                                    % bern_ComputeMAPandPrediction). 
in.opt.priorp2g1 = [prior_data.priorN2, prior_data.priorN1] + eps_prior ;
                                    % uniform Beta prior (when learning TP)
if strcmpi(learned_param,'alternation')
    in.opt.priorp1 = [1 1] + eps_prior ;  
                                    % uniform Beta prior (when learning IF, AF)
else
    in.opt.priorp1 = [prior_data.priorN1 prior_data.priorN2] + eps_prior ;  
end     
in.opt.prior_data = prior_data ;    % weighting of observations according to identity
in.verbose = 1;                     % to check that no default values are used.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ======================================================================= %
% ===== * ===== Load data
% ======================================================================= %
for i_subj = 1:n_subj
    subj_idx = num2str(all_subj_idx(i_subj)) ; 
    subj_idx_padd = num2str_padd(all_subj_idx(i_subj), 2) ; 
    param_str = [str_stim, '_subj',subj_idx] ; 
    fn_dir_subj = [fn_dir, 'subj',subj_idx_padd,'/'] ; 
    fn_cond_order = [fn_dir_subj,'Cond',param_str,'.mat'] ; 
    
    if exist(fn_cond_order, 'file')
        cond_data = load(fn_cond_order) ; 
        all_TP_test = cond_data.all_TP_test ;  
        TP_order_full = cond_data.TP_order_full ; 
        TP_order = cond_data.TP_order ; 
        order_used = TP_order_full ;
        % ================================
        
        n_conds = length(order_used) ;  
        for idx_cond = 1:n_conds
            TP_tmp = all_TP_test(order_used(idx_cond),:) ; 
            str_TP = get_str_TP(TP_tmp) ; 
            
            fn_test = [fn_dir_subj, 'Test_cond',num2str(idx_cond), param_str, ...
                str_TP, '_',num2str(n_stim_test),'stim_1.mat'] ; 
            if exist(fn_test, 'file')
                % ===== * ===== TP already found ? 
                idx_TP = find(strcmpi(TPs_found, str_TP),1) ; 
                new_TP = 0 ; 
                if isempty(idx_TP)
                    % Add the TP to the list
                    TPs_found{end+1} = str_TP ;   
                    all_TPs = [all_TPs; TP_tmp] ;
                    idx_TP = length(TPs_found) ; 
                    new_TP = 1 ; 
                end
                
                % ===== * ===== Load and extract useful information 
                if strcmpi(model_type,'io')
                    [n_data, n_data_conf, est_conf_tmp, est_p1_tmp, IO_p1_tmp, ...
                        IO_MAP_p1_tmp, IO_conf_tmp, trials_resp_TP, trials_resp_conf, ...
                        s, RT_tmp, RT_conf_tmp, full_IO_p1, full_IO_conf] = ...
                        load_extract_ratings(fn_test, remove_missed_resp, in) ; 
                    % proba_same --> est_p1 ; conf_proba --> est_conf
                elseif strcmpi(model_type, 'rw')
                    [n_data, n_data_conf, est_conf_tmp, est_p1_tmp, IO_p1_tmp, ...
                        IO_conf_tmp, trials_resp_TP, trials_resp_conf, s, ...
                        RT_tmp, RT_conf_tmp, full_IO_p1, full_IO_conf] = ...
                        load_extract_delta(fn_test, remove_missed_resp, learned_param, a_RW) ; 
                end
                
                % ===== * ===== Store useful information   
                if new_TP
                    est_p1{i_subj,idx_TP} = [] ; est_conf{i_subj, idx_TP} = [] ; 
                    IO_p1{i_subj,idx_TP} = [] ; IO_MAP_p1{i_subj,idx_TP} = [] ; 
                    IO_conf{i_subj,idx_TP} = [] ; last_I_p1{i_subj,idx_TP} = [] ; 
                    last_I_conf{i_subj,idx_TP} = [] ; 
                    next_I_p1{i_subj, idx_TP} = [] ; 
                    
                    RT_subj{i_subj,idx_TP} = [] ; RT_conf_subj{i_subj,idx_TP} = [] ;
                end
                est_p1{i_subj,idx_TP} = [est_p1{i_subj,idx_TP}, ...
                    reshape(est_p1_tmp,1,n_data)] ;
                est_conf{i_subj,idx_TP} = [est_conf{i_subj, idx_TP}, ...
                    reshape(est_conf_tmp,1,n_data_conf)] ;
                
                RT_subj{i_subj,idx_TP} = [RT_subj{i_subj,idx_TP}, ...
                    reshape(RT_tmp,1,n_data)] ; 
                RT_conf_subj{i_subj,idx_TP} = [RT_conf_subj{i_subj,idx_TP}, ...
                    reshape(RT_conf_tmp,1,n_data_conf)] ;

                IO_p1{i_subj,idx_TP} = [IO_p1{i_subj,idx_TP}, 100.*IO_p1_tmp] ; 
                IO_conf{i_subj,idx_TP} = [IO_conf{i_subj,idx_TP}, IO_conf_tmp] ;   
                
                last_I_p1{i_subj,idx_TP} = [last_I_p1{i_subj,idx_TP}, ...
                    reshape(s(trials_resp_TP),1,n_data)] ;       
                last_I_conf{i_subj,idx_TP} = [last_I_conf{i_subj,idx_TP}, ...
                    reshape(s(trials_resp_conf),1,n_data_conf)] ;
                
                if strcmpi(model_type,'io')
                    IO_MAP_p1{i_subj,idx_TP} = [IO_MAP_p1{i_subj,idx_TP}, 100.*IO_MAP_p1_tmp] ;
                end
                
                idx_next = trials_resp_TP + 1 ; 
                n_out = sum(idx_next>n_stim_test) ; 
                idx_next(idx_next>n_stim_test) = [] ; 
                next_I_p1{i_subj, idx_TP} = [next_I_p1{i_subj, idx_TP}, ...
                    reshape([reshape(s(idx_next), 1, n_data-n_out), zeros(1,n_out)],1,n_data)] ; 
                
                N_stim(i_subj,:) = N_stim(i_subj,:) + [sum(s==1), sum(s==2)] ; % nb of stim of I1 and I2
                N_rated(i_subj,:) = N_rated(i_subj,:) + [...
                    sum(s(trials_resp_TP)==1), sum(s(trials_resp_TP)==2)] ;
                
                
                % ========================
                % (n_out, n_stim_test, n_conds_max, n_subj) 
                ratings_IO_data(1, trials_resp_TP, idx_cond, i_subj) = est_p1_tmp ;
                ratings_IO_data(2, trials_resp_conf, idx_cond, i_subj) = est_conf_tmp ;
                ratings_IO_data(3, :, idx_cond, i_subj) = full_IO_p1 ;
                ratings_IO_data(4, :, idx_cond, i_subj) = full_IO_conf ;
                
            end
        end
    end    
end
save(fn_ratingsIO, '-v7.3', 'ratings_IO_data') ; 
% ======================================================================= %
% ===== * ===== Compute linear fit and correlation *per subject*
% ======================================================================= %

% Params of the linear fits
p1_poly = NaN(n_subj,2) ; MAP_p1_poly = NaN(n_subj,2) ; conf_poly = NaN(n_subj,2) ; 

% Correlations
all_corrs = NaN(n_subj,3) ;
all_MSE = NaN(n_subj,3) ; 


for i_subj=1:n_subj
    curr_est_p1 = [est_p1{i_subj,:}] ;
    curr_est_conf = [est_conf{i_subj,:}] ;
    
    curr_IO_p1 = [IO_p1{i_subj,:}] ;      
    curr_MAP_p1 = [IO_MAP_p1{i_subj,:}] ;
    curr_IO_conf = [IO_conf{i_subj,:}] ;
            
    if length(curr_est_p1)>1 % min 2 values to fit line
        all_corrs(i_subj,1) = corr(curr_IO_p1',curr_est_p1', 'rows','complete') ;
        curr_p1_poly = polyfit(curr_IO_p1,curr_est_p1,1) ;
        p1_poly(i_subj,:) = curr_p1_poly ;
        
        hat_p1 = polyval(curr_p1_poly,curr_IO_p1) ;                
        all_MSE(i_subj,1) = compute_mse_rmse(hat_p1',curr_est_p1') ;        
        
        if strcmpi(model_type,'io')
            all_corrs(i_subj,2) = corr(curr_MAP_p1',curr_est_p1', 'rows','complete') ; 
            curr_MAP_poly = polyfit(curr_MAP_p1,curr_est_p1,1) ;
            MAP_p1_poly(i_subj,:) = curr_MAP_poly ;
            all_MSE(i_subj,2) = compute_mse_rmse(polyval(curr_MAP_poly,curr_MAP_p1),...
                curr_est_p1) ; 
        end
        
    end
    
    if length(curr_est_conf)>1
        curr_conf_poly = polyfit(curr_IO_conf,curr_est_conf,1) ;
        conf_poly(i_subj,:) = curr_conf_poly ;

        all_corrs(i_subj,3) = corr(curr_IO_conf',curr_est_conf','rows','complete') ; 
        all_MSE(i_subj,3) = compute_mse_rmse(polyval(curr_conf_poly,curr_IO_conf),...
            curr_est_conf) ;
    end
    
end


if show_plots
    % ======================================================================= %
    % ===== * ===== Display data: scatter plots (+linear fits)
    % ======================================================================= %
    n_TPs = length(TPs_found) ; str_TPs = ['_',num2str(n_TPs),'TPs'] ;
    str_params = [str_model, str_subj, str_TPs] ;
    str_params_no_IO = [str_subj, str_TPs] ;
    if ~exist(fn_res, 'dir') % 0 or 7 if it exists
        mkdir(fn_res) ; 
    end
    % Plot params
    fig_sz = [1 2 12 12] ; 
    taille = 18 ; taille_axis = 14 ;     
    warm_col = [255,0,0]./255 ; cool_col = [0,0,179]./255 ; 
    mk_plot = {'o'} ; sz_pts = 9 ;
    show_data = 1 ; 
    
    single_fit = 0 ; % show all lin regs or just the avg
    if single_fit ; show_data = 1 ; end
    
    % == * == Significance of the correlation across subjects
    alpha_level = 0.05 ; paired_tests = 0 ; mean_H0 = 0 ; always_t_test = 0 ;
    two_sided = 1 ;
    str_acc = '' ; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Correlation btw TP estimation accuracy & confidence accuracy? 
    corr_proba = all_corrs(:,1) ; 
    corr_conf = all_corrs(:,3) ; 
    data_struct = struct() ; 
    poly = polyfit(corr_proba,corr_conf,1) ; 
    data_struct(1).x = corr_proba ;
    data_struct(1).y = corr_conf ;
    data_struct(1).coeffs = poly ;
    [rho_acc,pv_acc] = corr(corr_proba, corr_conf) ; 
    % here effect size = corr because we are interested in effect size of
    % the relationship btw x & y not size of the difference btw the 2
    % samples
    str_acc = ['* Corr btw the accuracies of probability estimates and confidence reports: ', ...
        num2str(rho_acc), ' (p = ',num2str(pv_acc),')'] ; 
    disp(str_acc)
    
    curr_fn = [fn_res, 'Estim_acc', str_params] ;
    scatter_colormap_lin_fit(data_struct, 'eps_fig', eps_fig, 'fig_fig', fig_fig, ...
        'pdf_fig', pdf_fig, 'fn_save', curr_fn, 'save_fig',1, 'fig_name', ...
        'estimation acc', 'add_lin_fit', 1, 'use_colormap', 0, 'x_lab', ...
        'proba. estimation accuracy', 'y_lab', 'conf. accuracy', 'fig_sz', fig_sz, 'taille_tick', ...
        taille, 'taille_stick', taille_axis, 'sz_pts', sz_pts, ...
        'markers', mk_plot, 'single_fit', single_fit, 'add_lgd', 0, ...
        'aggregated_corr', 1, 'axis_equal', 1, ...
        'add_identity', 0, 'c_pts', [0,0,0], 'show_data', 1) ; %, 'xL', [0.05,0.9]) ; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % xdata & ydata: should be n_subj x n_TP cells
    fn_txt_res = [fn_res, 'All_corr', str_params, '.txt'] ;
    %%%%%%%%%% empty the file fn_txt_res 
    write_empty_lines(fn_txt_res, 0,'w') ;
    write_to_txt(fn_txt_res, str_acc,'a',2) ;
    
    % Check the fit depending on the type of transition?
    % --- Based on last I received
    % --- OR next stim to be received (when exist)
    
    IO_p1_g1 = cell(n_subj,n_TPs) ; 
    IO_p1_g2 = cell(n_subj,n_TPs) ; 
    est_p1_g1 = cell(n_subj,n_TPs) ; 
    est_p1_g2 = cell(n_subj,n_TPs) ; 
    
    IO_p1_n1 = cell(n_subj,n_TPs) ; % next = 1
    IO_p1_n2 = cell(n_subj,n_TPs) ; % next = 2
    est_p1_n1 = cell(n_subj,n_TPs) ; 
    est_p1_n2 = cell(n_subj,n_TPs) ; 
        
    % Generative proba
    gen_p1 = cell(n_subj,n_TPs) ;
    gen_p1_g1 = cell(n_subj,n_TPs) ; gen_p1_g2 = cell(n_subj,n_TPs) ; 
    all_p1 = compute_p1(all_TPs(:,1), all_TPs(:,2)) ; 
    
    IO_pnext = IO_p1 ; % proba to receive next stim
    est_pnext = est_p1 ; 
    
    for i_s = 1:n_subj
        for j_tp = 1:n_TPs
            IO_tmp = IO_p1{i_s, j_tp} ; 
            est_tmp = est_p1{i_s,j_tp} ; 
            s_tmp = last_I_p1{i_s,j_tp} ; 
            
            IO_p1_g1{i_s,j_tp} = IO_tmp(s_tmp==1) ; 
            IO_p1_g2{i_s,j_tp} = IO_tmp(s_tmp==2) ; 
            est_p1_g1{i_s,j_tp} = est_tmp(s_tmp==1) ; 
            est_p1_g2{i_s,j_tp} = est_tmp(s_tmp==2) ;
            
            s_next_tmp = next_I_p1{i_s, j_tp} ; 
            IO_p1_n1{i_s,j_tp} = IO_tmp(s_next_tmp==1) ; 
            IO_p1_n2{i_s,j_tp} = IO_tmp(s_next_tmp==2) ; 
            est_p1_n1{i_s,j_tp} = est_tmp(s_next_tmp==1) ; 
            est_p1_n2{i_s,j_tp} = est_tmp(s_next_tmp==2) ;
            
            TP_tmp = 100.*all_TPs(j_tp,:) ; p1 = 100.*all_p1(j_tp) ; 
            gen_p1{i_s, j_tp} = repmat(p1, 1, length(est_tmp)) ; 
            gen_p1_g1{i_s, j_tp} = repmat(100-TP_tmp(2), 1, length(est_p1_g1{i_s,j_tp})) ; 
            gen_p1_g2{i_s, j_tp} = repmat(TP_tmp(1), 1, length(est_p1_g2{i_s,j_tp})) ; 
            
            IO_pnext_tmp = IO_tmp ; 
            IO_pnext_tmp(s_next_tmp==2) = 100 - IO_tmp(s_next_tmp==2) ; 
            IO_pnext{i_s,j_tp} = IO_pnext_tmp ; 
            est_pnext_tmp = est_tmp ; 
            est_pnext_tmp(s_next_tmp==2) = 100 - est_tmp(s_next_tmp==2) ; 
            est_pnext{i_s,j_tp} = est_pnext_tmp ; 
        end
    end
    
    all_fn_names = {'RT_est_conf', 'RT_conf', 'Fit_post', ...
        'Fit_conf', 'Gen_p1_g1', 'Gen_p1_g2', 'Gen_p1', ...
        'Conf_RT_est_conf', 'Conf_RT_conf', ...
        'Fit_post_g1', 'Fit_post_g2', 'Fit_post_n1', 'Fit_post_n2', ...
        'p1_est_conf', 'est_p1_est_conf', 'pnext_est_conf', 'est_pnext_est_conf', 'p1_conf', 'pnext_conf', 'Fit_MAP'} ;
    all_gen_params = {str_params_no_IO, str_params, str_params, ...
        str_params, str_params_no_IO, str_params_no_IO, str_params_no_IO, str_params_no_IO, str_params, ...
        str_params, str_params, str_params, str_params, ...
        str_params, str_params_no_IO, str_params, str_params_no_IO, str_params, str_params, str_params} ; 
    all_xdata_names = {'est_conf', 'IO_conf', 'IO_p1', ...
        'IO_conf', 'gen_p1_g1', 'gen_p1_g2', 'gen_p1', ...
        'est_conf', 'IO_conf', ...
        'IO_p1_g1', 'IO_p1_g2', 'IO_p1_n1', 'IO_p1_n2', ...
        'IO_p1', 'est_p1', 'IO_pnext', 'est_pnext', 'IO_p1', 'IO_pnext', 'IO_MAP_p1'} ;
    all_ydata_names = {'RT_subj', 'RT_subj', 'est_p1', ...
        'est_conf', 'est_p1_g1', 'est_p1_g2', 'est_p1', ...
        'RT_conf_subj', 'RT_conf_subj', ...
        'est_p1_g1', 'est_p1_g2', 'est_p1_n1', 'est_p1_n2', ...
        'est_conf', 'est_conf', 'est_conf', 'est_conf', 'IO_conf', 'IO_conf', 'est_p1'} ;  
    all_figs_names = {'Corr proba. RT - est. conf', 'Corr proba. RT - IO conf', 'Estimations of the TPs', ...
        'Confidence estimates', 'Gen vs. est p1g1', 'Gen vs. est p1g2', 'Gen vs. est p1', ...
        'Corr conf. RT - est. conf', ...
        'Corr conf. RT - IO conf', 'Estimations of the TPs from I1', 'Estimations of the TPs from I2', ...
        'Estimations of the TPs to I1', 'Estimations of the TPs to I2', ...
        'Bayesian proba estimates - conf estimates', 'Proba estimates - conf estimates', ...
        'pnext - est conf.', 'est. pnext - est conf.', 'Bayesian proba estimates - conf', 'Bayesian proba estimates - conf', 'Estimations of the TPs (MAP)'} ;     
    all_xlab_names = {'confidence estimates', 'IO confidence', 'Bayesian proba. estimates', ...
        'Bayesian confidence', 'true p(I_1|I_1)', 'true p(I_1|I_2)', 'true p(I_1)', ...
        'Confidence estimates', ...
        'IO confidence', 'IO proba. estimates', 'IO proba. estimates', 'IO proba. estimates', 'IO proba. estimates', ...
        'Bayesian proba. estimates', 'subjective proba. estimates', 'Bayesian proba. estimates', 'subjective proba. estimates', ...
        'Bayesian proba. estimates', 'Bayesian proba. estimates', 'IO MAP proba. estimates'} ; 
    all_ylab_names = {{'RT for proba. estimates'}, {'RT for proba. estimates'}, {'subjective proba. estimates'}, ...
        {'confidence estimates'}, {'rated p(I_1|I_1)'}, {'rated p(I_1|I_2)'}, {'rated p(I_1)'}, ...
        {'RT for conf. estimates'}, ...
        {'RT for conf. estimates'}, {'Subjective proba. estimates'}, {'Subjective proba. estimates'}, ...
        {'Subjective proba. estimates'}, {'Subjective proba. estimates'}, ...
        {'confidence estimates'}, {'confidence estimates'}, {'confidence estimates'}, ...
        {'confidence estimates'}, {'Bayesian confidence'}, {'Bayesian confidence'}, {'Subjective proba. estimates'}} ;
    all_axis_equal = [0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1] ; 
    all_add_identity = [0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1] ; 
    plot_box = [0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] ; % when discrete x values
    x_box = {[], [], [], [], 100.*unique(1-all_TPs(:,2)), 100.*unique(all_TPs(:,1)), 100.*unique(all_p1), ...
         [], [], [], [], [], [], [], [], [], [], [], [], []} ; 
    
    order_fit = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 1] ; % polynomial order of the fitted curves

    name_corr = {'RT - est. conf', 'RT - IO conf', 'est. proba - IO proba', ...
        'est. conf - IO conf', 'gen vs. est p1g1', 'gen vs. est p1g2', 'gen vs. est p1', ...
        'conf. RT - est. conf', ...
        'conf. RT - IO conf', 'est. proba - IO proba g1', 'est. proba - IO proba g2', ...
        'est. proba - IO proba n1', 'est. proba - IO proba n2', ...
        'IO proba - est conf.', 'est proba - est conf', 'IO pnext - est conf.', 'est pnext - est conf', 'IO proba - IO conf.', 'IO pnext - IO conf.', 'est. proba - IO MAP proba'} ; 
    
    if strcmpi(model_type,'io')
        n_corrs = length(all_xdata_names) ;
    elseif strcmpi(model_type,'rw')
        n_corrs = length(all_xdata_names) - 1 ; % don't do last fit = with the MAP
    end
                    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    aggregated_corr = 0 ; % for title: single corr for all plotted data OR mean of the subject corr
    
    
    for idx_corr= [3, 4, 5, 6, 7, 15, 18] % 1:n_corrs
        curr_xdata = eval(all_xdata_names{idx_corr}) ; 
        curr_ydata = eval(all_ydata_names{idx_corr}) ;
        curr_fn_lab = all_fn_names{idx_corr} ; 
        curr_param = all_gen_params{idx_corr} ; 
        curr_xlab = all_xlab_names{idx_corr} ; 
        curr_ylab = all_ylab_names{idx_corr} ; 
        fig_name = all_figs_names{idx_corr} ; 
        axis_equal = all_axis_equal(idx_corr) ; 
        add_identity = all_add_identity(idx_corr) ; 
        show_box = plot_box(idx_corr) ; 
        if show_box 
            curr_show_data = 0 ; curr_single_fit = 1 ; 
        else
            curr_show_data = show_data ; curr_single_fit = single_fit ; 
        end 
        curr_x_box = x_box{idx_corr} ; 
        curr_order_fit = order_fit(idx_corr) ; 
        
        curr_fn = [fn_res, curr_fn_lab, curr_param] ;
        
        if idx_corr==1
            % Count once the number of times the ratings are asked to the
            % subjects
           n_data_per_subj = NaN(n_subj,1) ;  
           for idx_subj=1:n_subj
               n_data_per_subj(idx_subj) = length([curr_xdata{idx_subj,:}]) ; 
           end
           str_tmp = ['Mean (std) data points per subject: ', num2str(round(mean(n_data_per_subj),3)), '(',...
                num2str(round(std(n_data_per_subj),3)),')'] ; 
           disp(str_tmp) ; write_to_txt(fn_txt_res, str_tmp,'a',2) ;
        end
        
        data_struct = struct() ; 
        indiv_poly = NaN(n_subj,curr_order_fit+1) ;
        curr_corr = NaN(n_subj,1) ; 
        for idx_subj=1:n_subj
            curr_x = [curr_xdata{idx_subj,:}] ;
            curr_y = [curr_ydata{idx_subj,:}] ;
            
            curr_poly = polyfit(curr_x,curr_y,curr_order_fit) ;            
            indiv_poly(idx_subj,:) = curr_poly ; 
            if curr_order_fit==1
                curr_corr(idx_subj) =  corr(curr_x',curr_y','rows','complete') ; 
            else
                curr_y_hat = polyval(curr_poly,curr_x) ; 
                % 2 options, can compute:
                % 1. the coefficient of determination (R^2)
                mse_tmp = compute_mse_rmse(curr_y_hat,curr_y) ; 
                curr_corr(idx_subj) = 1 - mse_tmp./ mean( (curr_y-mean(curr_y)).^2 ) ; 
                
                % 2. correlation target - model output (better to see if negative)                
                %curr_corr(idx_subj) =  corr(curr_y_hat',curr_y','rows','complete') ; 
            end
            
            data_struct(idx_subj).x = curr_x ;
            data_struct(idx_subj).y = curr_y ;
            
            %if show_box ; curr_x_box = curr_x_box end
        end 
        mean_poly = nanmean(indiv_poly,1) ; 
        for idx_subj=1:n_subj
            if curr_single_fit
                data_struct(idx_subj).coeffs = mean_poly ;
            else
                data_struct(idx_subj).coeffs = indiv_poly(idx_subj,:) ;
            end
        end
        
        scatter_colormap_lin_fit(data_struct, 'eps_fig', eps_fig, 'fig_fig', fig_fig, ...
            'pdf_fig', pdf_fig, 'fn_save', curr_fn, 'save_fig',1, 'fig_name', ...
            fig_name, 'add_lin_fit', 1, 'use_colormap', 0, 'x_lab', ...
            curr_xlab, 'y_lab', curr_ylab, 'fig_sz', fig_sz, 'taille_tick', ...
            taille, 'taille_stick', taille_axis, 'sz_pts', sz_pts, ...
            'markers', mk_plot, 'single_fit', curr_single_fit, 'add_lgd', 0, ...
            'aggregated_corr', aggregated_corr, 'axis_equal', axis_equal, ...
            'add_identity', add_identity, 'show_data', curr_show_data, ...
            'show_box', show_box, 'x_box',curr_x_box) ; 
        
        % ======================================================================= %
        % ===== * ===== Correlations
        % ======================================================================= %
        curr_corr_name = name_corr{idx_corr} ; 
        if n_subj>2
            % significance of corr: t-tests vs. 0 

            [p_values, ~, ~, ~, ~, test_statistic, cohend] = ...
                compute_significance_matrix(curr_corr, alpha_level, 10, ...
                paired_tests, mean_H0, always_t_test, two_sided) ;

            str_tmp = ['Corr ',curr_corr_name,': ', num2str(round(nanmean(curr_corr),3)), ' (',...
                num2str(round(test_statistic,3)), ', p = ',...
                num2str(round(p_values,5)),', Cohen''s d = ', num2str(round(cohend,5)),')'] ; 
        else
            str_tmp = ['Corr ', curr_corr_name, ': ', num2str(round(nanmean(curr_corr),3))] ; 
        end
        disp(str_tmp) ; write_to_txt(fn_txt_res, str_tmp,'a',1) ;
    end
    
    % ======================================================================= %
    % ===== * ===== Sanitary checks
    % ======================================================================= %
    if show_sanit_checks && n_subj>=2
        curr_fn = [fn_res, 'N_stim', str_params_no_IO] ;
        plot_boxplot(N_stim, 'eps_fig', eps_fig, 'fig_fig', fig_fig, ...
            'pdf_fig', pdf_fig, 'fn_save', curr_fn, 'save_fig',1, 'fig_name', ...
            'Nb of stim received', 'x_lab', '', 'y_lab', 'count', 'labels', {'I_1', 'I_2'}) ; 

        curr_fn = [fn_res, 'N_rated', str_params_no_IO] ;
        plot_boxplot(N_rated, 'eps_fig', eps_fig, 'fig_fig', fig_fig, ...
            'pdf_fig', pdf_fig, 'fn_save', curr_fn, 'save_fig',1, 'fig_name', ...
            'Nb of stim rated', 'x_lab', '', 'y_lab', 'count', 'labels', {'I_2|I_1', 'I_1|I_2'}) ; 
    end
    if n_subj>2
        [p_values, h_mask, ~, ~, conf_intervals, test_statistic, cohend] = ...
            compute_significance_matrix(N_stim, alpha_level, 10, ...
            1, mean_H0, always_t_test, two_sided) ;

        str_tmp = ['Difference btw nb of I_1 and nb of I_2: ', ...
            num2str(round(diff(nanmean(N_stim)),3)), ' (',...
            num2str(round(test_statistic(2,1),3)), ', p = ',...
            num2str(round(p_values(2,1),5)),', Cohen''s d = ',num2str(round(cohend(2,1),5)),')'] ; 
        disp(str_tmp) ; write_to_txt(fn_txt_res, str_tmp,'a',1) ;
        
        [p_values, h_mask, ~, ~, conf_intervals, test_statistic, cohend] = ...
            compute_significance_matrix(N_rated, alpha_level, 10, ...
            1, mean_H0, always_t_test, two_sided) ;

        str_tmp = ['Difference btw nb of transitions rated from I_1 and from I_2: ', ...
            num2str(round(diff(nanmean(N_rated)),3)), ' (',...
            num2str(round(test_statistic(2,1),3)), ', p = ',...
            num2str(round(p_values(2,1),5)),', Cohen''s d =',num2str(round(cohend(2,1)),5),')'] ; 
        disp(str_tmp) ; write_to_txt(fn_txt_res, str_tmp,'a',1) ;
    end
end


end
