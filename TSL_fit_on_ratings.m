function [all_params] = TSL_fit_on_ratings( )
% Compute the fit of different models (with different parameters) for the
% behavioral ratings. 
% 

% ======================================================================= %
% ===*=== Parameters
compute_data        = 0 ;   % otherwise, reload the fits if data previously computed
save_each_step      = 0 ;   % can take up to 5-10 sec per param
reload_when_exist   = 1 ;   % reload data when already saved for the current parameter

show_plots_per_param= 0 ;   % produce the plots in TSL_analyze_ratings.m  
                            % for each tested parameters
all_subj_idx        = setdiff(1:36,[1, 11, 15, 28, 33]) ; 
% --------------------------------------------------------- % 
remove_AF           = 0 ;   % Remove AF from the BMC.  
use_PEP             = 0 ;   % protected exceedance probability 
% correct EPs for the possiblity that observed differences in model 
% evidences (over subjects) are due to chance (Rigoux et al., 2014). 

% ======================================================================= %
% ===== * ===== RW + IO: Default parameters
% ======================================================================= %
fit_confidence      = 0 ; % select the model which maximizes fit of confidence instead of proba
n_subj              = length(all_subj_idx) ; 
% ======================================================================= %
% ===== * ===== RW: Default parameters
% ======================================================================= %
all_params_RW       = get_all_params_RW() ; 
n_params_RW         = length(all_params_RW) ; 

base_fn = './results_ratings_d/' ; 
fn_res_RW = [base_fn, 'all_fits/'] ; 
if fit_confidence ; fn_best_RW = [base_fn, 'best_params_conf/'] ; 
else ; fn_best_RW = [base_fn, 'best_params/'] ; end
if ~exist(fn_res_RW, 'dir') ; mkdir(fn_res_RW) ; end
if ~exist(fn_best_RW, 'dir') ; mkdir(fn_best_RW) ; end

str_fit = '' ; 
if fit_confidence ; str_fit = '_opt_conf' ; end
str_fit_fig = str_fit ; str_AF = '' ; 
if remove_AF ; str_fit_fig = [str_fit_fig, '_noAF'] ; str_AF = '_noAF' ; end

names_models_RW = unique({all_params_RW(:).learned_param}) ; n_models_RW = length(names_models_RW) ; 
param_mod = ['_', num2str(n_models_RW), 'models'] ;
str_genp_RW = ['_', num2str(n_subj), 'subj_', num2str(n_params_RW), 'params', param_mod] ; 
str_params_fits_RW = [str_fit_fig, str_genp_RW] ; % for the plots
fn_all_data_RW = [fn_res_RW, 'Fit',str_genp_RW,'.mat'] ;

% ======================================================================= %
% ===== * ===== Ideal Observer: Default parameters
% ======================================================================= %
show_best_params    = 0 ; 
consider_flip       = 0 ; % flipped (forward) decays vs. backward decays 
combine_decay       = 0 ; % Combine the forward AND backward decays, with 2 distinct leaks
relative_weight     = 0 ; % when combine_decay, also consider a relative weigh btw fwd and bwd leaks
show_sem            = 1 ; % for posterior of leaks
% Parameters for the heatmap if combine_decay
map_opt = 1 ;
% option for map: 1) imagesc, 2) contourf (filled contours), 3) contour
interp_and_log = 1 ;    % log scales for the x and y axes
log_data = 1 ;          % log scale for the proba (data in the map)
%%%%%%%%%%%%%%%%%%%%%%%%%
str_map = ['_map', num2str(map_opt)] ; 
if interp_and_log ; str_map = [str_map, '_logAx'] ; end
if log_data ; str_map = [str_map, '_log'] ; end   

if combine_decay
    consider_flip = 0 ; % both decays change separately
end

all_params          = get_all_params_to_test(consider_flip, combine_decay,relative_weight) ; 
n_params            = length(all_params) ; 

fn_res = ['./results_ratings/all_fits/'] ; 
fn_dat = ['./results_ratings/data/'] ; 
if fit_confidence
    fn_best_params = ['./results_ratings/best_params_conf/'] ; 
else
    fn_best_params = ['./results_ratings/best_params/'] ; 
end
if ~exist(fn_res, 'dir') ; mkdir(fn_res) ; end % 0 or 7 if it exists
if ~exist(fn_dat, 'dir') ; mkdir(fn_dat) ; end

if ~exist(fn_best_params, 'dir') ; mkdir(fn_best_params) ; end
str_comb = '' ; 
if combine_decay ; str_comb = '_Comb' ; if relative_weight
    str_comb = [str_comb, '_weight'] ; end ; end
if consider_flip ; str_comb = '_flip' ; end

names_models = unique({all_params(:).learned_param}) ; n_models = length(names_models) ; 
param_mod = ['_', num2str(n_models), 'models'] ;
str_gen_param = ['_', num2str(n_subj), 'subj_', num2str(n_params), 'params', ...
    param_mod, str_comb] ; 
str_params_fits = [str_fit, str_gen_param] ;
fn_all_data = [fn_res, 'Fit',str_params_fits,'.mat'] ; % best params dpd on fit_confidence
str_params_fits = [str_fit_fig, str_gen_param] ; % for the plots

% Params when comparing IO and RW together
str_genp_IO_RW = ['_', num2str(n_subj), 'subj_', num2str(n_params_RW+n_params), ...
    'params_', num2str(n_models+n_models_RW), 'models'] ; 
str_fits_IO_RW = [str_fit_fig, str_genp_IO_RW] ; % for the plots


% General plot params
eps_fig = 0 ; fig_fig = 0 ; pdf_fig = 0 ;
taille = 18 ; taille_axis = 14 ; 

if compute_data
    % =================================================================== % 
    % ===*=== Compute data: IO
    % =================================================================== %
    MSE_all_subj        = NaN(n_subj,3, n_params) ; 
    % need to have all MSE before computing the BIC since need to select
    % best decay w. 
    corr_all_subj       = NaN(n_subj,3, n_params) ; 

    % Data structures to keep IO data leading to the best MSE for each
    % subject
    best_IO_p1          = cell(n_subj,1) ;      % IO p1 (mean of posterior)
    best_IO_MAP_p1      = cell(n_subj,1) ;      % IO p1 (MAP of posterior)
    best_IO_conf        = cell(n_subj,1) ;      % IO conf
    
    best_p1_poly        = NaN(n_subj,2) ; 
    best_MAP_p1_poly    = NaN(n_subj,2) ; 
    best_conf_poly      = NaN(n_subj,2) ; 

    best_MSE            = inf*ones(n_subj,1) ;  % MSE: should be minimized
    idx_best_params     = NaN(n_subj,1) ; 
    % =================================================================== %
    for idx_param=1:n_params
        
        AboutFirst = all_params(idx_param).AboutFirst ; 
        MemParam = all_params(idx_param).MemParam ; 
        learned_param = all_params(idx_param).learned_param ;         
        str_model = all_params(idx_param).str_model ; 
       
        fn_data_param = [fn_res, 'Fit_',num2str(n_subj), 'subj', str_model,'.mat'] ;
        
        if reload_when_exist && exist(fn_data_param, 'file')
            disp(['--- Reloading fit for param set ', num2str(idx_param), '/', ...
                num2str(n_params),' ', str_model, ' ---'])
            data_param = load(fn_data_param) ;
            all_corrs = data_param.all_corrs;
            all_MSE = data_param.all_MSE ;
            est_p1 = data_param.est_p1 ;
            IO_p1 = data_param.IO_p1 ;
            IO_MAP_p1 = data_param.IO_MAP_p1 ;
            est_conf = data_param.est_conf ;
            IO_conf = data_param.IO_conf ;
            last_I_p1 = data_param.last_I_p1 ;
            last_I_conf = data_param.last_I_conf ;
            TPs_found = data_param.TPs_found ;
            p1_poly = data_param.p1_poly ; 
            MAP_p1_poly = data_param.MAP_p1_poly ; 
            conf_poly = data_param.conf_poly ; 
        else

            disp(['--- Computing fit for param set ', num2str(idx_param), '/', ...
                num2str(n_params),' ',str_model,' ---'])
            t_st = tic ; 
            [all_corrs, all_MSE, est_p1, IO_p1, IO_MAP_p1, est_conf, IO_conf, ...
                last_I_p1, last_I_conf, TPs_found, p1_poly, MAP_p1_poly, conf_poly] = ...
                TSL_analyze_ratings('show_plots',show_plots_per_param, ...
                'MemParam', MemParam, 'AboutFirst', AboutFirst, 'learned_param', ...
                learned_param, 'all_subj_idx', all_subj_idx, 'model_type', 'IO') ;            
            
            % Save data for the current param
            save(fn_data_param, '-v7.3', 'all_corrs', 'all_MSE', 'est_p1', ...
                'IO_p1', 'IO_MAP_p1', 'est_conf', 'IO_conf', 'last_I_p1', ...
                'last_I_conf', 'TPs_found', 'p1_poly', 'MAP_p1_poly', 'conf_poly') ;
            t_elps = toc(t_st) ; 
            disp(['Time elapsed for ', num2str(n_subj),' subj: ', num2str(t_elps)])
        end
        
        MSE_all_subj(:,:,idx_param) = all_MSE ;        
        corr_all_subj(:,:,idx_param) = all_corrs ; 
        
        % ------- Update the best MSE per subject?
        if fit_confidence
            curr_MSE = all_MSE(:,3) ; 
        else
            curr_MSE = all_MSE(:,1) ; 
        end
        idx_new_best = best_MSE>curr_MSE ; 

        best_MSE(idx_new_best) = curr_MSE(idx_new_best) ; 
        idx_best_params(idx_new_best) = idx_param ; 
        n_TPs = size(IO_p1, 2) ; 
        best_IO_p1(idx_new_best,1:n_TPs) = IO_p1(idx_new_best,:) ; 
        best_IO_MAP_p1(idx_new_best,1:n_TPs) = IO_MAP_p1(idx_new_best,:) ; 
        best_IO_conf(idx_new_best,1:n_TPs) = IO_conf(idx_new_best,:) ; 
        
        best_p1_poly(idx_new_best,:) = p1_poly(idx_new_best,:) ;
        best_MAP_p1_poly(idx_new_best,:) = MAP_p1_poly(idx_new_best,:) ;
        best_conf_poly(idx_new_best,:) = conf_poly(idx_new_best,:) ;
        
        if save_each_step
            save(fn_all_data, '-v7.3', 'MSE_all_subj', 'corr_all_subj', 'all_params', ...
                'all_subj_idx', 'best_IO_p1', 'best_IO_MAP_p1', 'best_IO_conf', ...
                'idx_best_params', 'est_p1', 'est_conf', 'TPs_found', 'last_I_p1', ...
                'last_I_conf', 'best_p1_poly', 'best_MAP_p1_poly', 'best_conf_poly') ;
        end
    end 
    save(fn_all_data, '-v7.3', 'MSE_all_subj', 'corr_all_subj', 'all_params', ...
        'all_subj_idx', 'best_IO_p1', 'best_IO_MAP_p1', 'best_IO_conf', ...
        'idx_best_params', 'est_p1', 'est_conf', 'TPs_found', 'last_I_p1', ...
        'last_I_conf', 'best_p1_poly', 'best_MAP_p1_poly', 'best_conf_poly') ;

    % =================================================================== % 
    % ===*=== Compute data: RW
    % =================================================================== %
    MSE_RW      = NaN(n_subj, 2, n_params_RW) ; 
    corr_RW     = NaN(n_subj, 2, n_params_RW) ; 
    % =================================================================== %
    
    for idx_param=1:n_params_RW
        a_RW = all_params_RW(idx_param).a_RW ; 
        learned_param = all_params_RW(idx_param).learned_param ;         
        str_model = all_params_RW(idx_param).str_model ; 
       
        fn_data_param = [fn_res_RW, 'Fit_',num2str(n_subj), 'subj', str_model,'.mat'] ;
        
        if reload_when_exist && exist(fn_data_param, 'file')
            disp(['--- Reloading fit for param set ', num2str(idx_param), '/', ...
                num2str(n_params_RW),' ', str_model, ' ---'])
            data_param = load(fn_data_param) ;
            all_corrs = data_param.all_corrs;
            all_MSE = data_param.all_MSE ;
        else

            disp(['--- Computing fit for param set ', num2str(idx_param), '/', ...
                num2str(n_params_RW),' ',str_model,' ---'])
            t_st = tic ; 
            [all_corrs, all_MSE] = TSL_analyze_ratings('show_plots',show_plots_per_param, ...
                'learned_param', learned_param, 'all_subj_idx', all_subj_idx, ...
                'model_type', 'RW', 'a_RW', a_RW) ;            
            
            % Save data for the current param
            save(fn_data_param, '-v7.3', 'all_corrs', 'all_MSE') ;
            t_elps = toc(t_st) ; 
            disp(['Time elapsed for ', num2str(n_subj),' subj: ', num2str(t_elps)])
        end
        
        MSE_RW(:,:,idx_param) = all_MSE(:,[1,3]) ;        
        corr_RW(:,:,idx_param) = all_corrs(:,[1,3]) ; 

    end 
    save(fn_all_data_RW, '-v7.3', 'MSE_RW', 'corr_RW', 'all_params_RW', 'all_subj_idx') ;
    
    
else
    
    % =================================================================== %
    % ===*=== Reload data: IO
    % =================================================================== %
    if ~exist(fn_all_data, 'file')
        error(['   ==== Data not found ====   '])
    else
        data_reloaded   = load(fn_all_data) ; 
        MSE_all_subj    = data_reloaded.MSE_all_subj ; %./10000 ; 
        corr_all_subj   = data_reloaded.corr_all_subj ; 
        all_params      = data_reloaded.all_params ; 
        all_subj_idx    = data_reloaded.all_subj_idx ; 
        best_IO_p1      = data_reloaded.best_IO_p1 ; 
        best_IO_MAP_p1  = data_reloaded.best_IO_MAP_p1 ; 
        best_IO_conf    = data_reloaded.best_IO_conf ; 
        idx_best_params = data_reloaded.idx_best_params ;
        est_p1          = data_reloaded.est_p1 ; 
        est_conf        = data_reloaded.est_conf ; 
        TPs_found       = data_reloaded.TPs_found ; 
        last_I_p1       = data_reloaded.last_I_p1 ; 
        last_I_conf     = data_reloaded.last_I_conf ; 
        best_p1_poly    = data_reloaded.best_p1_poly ; 
        best_MAP_p1_poly= data_reloaded.best_MAP_p1_poly ; 
        best_conf_poly  = data_reloaded.best_conf_poly ; 
        
        mean_p1_poly = nanmean(best_p1_poly,1) ; % [1,1] ;%
        mean_MAP_p1_poly = nanmean(best_MAP_p1_poly,1) ; % [1,1] ; %
        mean_conf_poly = nanmean(best_conf_poly,1) ; % [1,1] ; %
    end
    % =================================================================== %
    % ===*=== Reload data: RW
    % =================================================================== %
    if ~exist(fn_all_data_RW, 'file')
        error(['   ==== RW data not found ====   '])
    else
        data_reloaded   = load(fn_all_data_RW) ; 
        MSE_RW    = data_reloaded.MSE_RW ; %./10000 ; % MSE btw the proba (remove *100 factor) 
        corr_RW   = data_reloaded.corr_RW ; 
        all_params_RW   = data_reloaded.all_params_RW ; 
        all_subj_idx    = data_reloaded.all_subj_idx ; 
    end
    
    
    
    % =============================================================== %
    % ===*=== Display results
    % =============================================================== %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% (1) best parameter *per subject* (IO only)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % ------- Distribution of the best leaks   
    % --------------------------------------------------------------- %
    % histograms of the time windows selected for each learned_param

    fn_hist_leak = [fn_best_params, 'Distrib_decays', str_params_fits] ;
    fn_txt_res = [fn_best_params, 'Opt_params', str_params_fits, '.txt'] ;
    % empty the file fn_txt_res 
    write_empty_lines(fn_txt_res, 0,'w') ;
    fig_sz = [1 2 18 12] ;

    learned_param = {all_params(:).learned_param} ; 
    names_models = unique(learned_param) ; 
    if remove_AF 
        names_models = names_models(~strcmpi(names_models, 'alternation')) ; 
    end    
    n_models = length(names_models) ;  

    all_leaky_int = logical([all_params(:).leaky_int]) ;
    all_leakyF_int = logical([all_params(:).leakyF_int]) ;

    all_decays = [all_params(:).decay] ;
    unique_decays = unique(all_decays(all_leaky_int)) ; 
    unique_decays_true = unique_decays ; % with inf values
    n_decays = length(unique_decays) ;

    all_decayFs = [all_params(:).decayF] ; 
    unique_decayFs = unique(all_decayFs(all_leakyF_int)) ; 
    n_decayFs = length(unique_decayFs) ;  

    sorted_decays = sort(unique_decays) ; sorted_decays(sorted_decays==Inf) = [] ;

    %%%%%%%%%%% Remove Inf values for the hist and other plots
    decay_val_for_inf = sorted_decays(end) + ...
        2*(sorted_decays(end)-sorted_decays(end-1)) ;
    idx_decay_inf = unique_decays==Inf ;
    idx_decay_0 = unique_decays==0 ; 
    unique_decays(idx_decay_inf) = decay_val_for_inf ; 
    all_decays(all_decays==Inf) = decay_val_for_inf ; 
    idx_decayF_inf = unique_decayFs==Inf ; 
    idx_decayF_0 = unique_decayFs==0 ; 
    unique_decayFs(idx_decayF_inf) = decay_val_for_inf ; 
    all_decayFs(all_decayFs==Inf) = decay_val_for_inf ; 

    
    best_learned_params = learned_param(idx_best_params) ;
    best_decays = all_decays(idx_best_params) ;
    if show_best_params
        % Find the optimal decays and learned param for each subject 
        best_learned_params = learned_param(idx_best_params) ; 

        best_decays = all_decays(idx_best_params) ; 
        best_decayFs = all_decayFs(idx_best_params) ; 
        best_leakyF_int = all_leakyF_int(idx_best_params) ; 
        best_leaky_int = all_leaky_int(idx_best_params) ;

        best_decays_neg_flipped = best_decays ; 
        best_decays_neg_flipped(best_leakyF_int & ~best_leaky_int) = ...
            - best_decayFs(best_leakyF_int & ~best_leaky_int) ; 

        data_struct = struct() ;
        for idx_model = 1:n_models
            idx_curr_model = cellfun(@(x) strcmpi(names_models{idx_model},x), ...
                best_learned_params) ; 
            if consider_flip
                % Adapt the values of the best decays when they are flipped                
                data_struct(idx_model).x = best_decays_neg_flipped(idx_curr_model) ;            
            elseif combine_decay
                data_struct(idx_model).x = best_decays(idx_curr_model) ; 
                data_struct(n_models + idx_model).x = best_decayFs(idx_curr_model) ; 
            else
                % only backward decays
                data_struct(idx_model).x = best_decays(idx_curr_model) ;                
            end
        end

        if consider_flip
            caption_lgd = names_models ;
        elseif combine_decay
            caption_lgd = cell(1, 2*n_models) ;
            caption_lgd(1:n_models) = names_models ; 
            caption_lgd(n_models+1:end) = strcat(names_models, '-F') ;
        else
            caption_lgd = names_models ; 
        end

        my_histograms(data_struct, 'eps_fig', eps_fig, 'fig_fig', fig_fig, 'pdf_fig', pdf_fig, ...
            'fn_save', fn_hist_leak, 'save_fig',1, 'fig_name', 'Distribution of the decays',...
            'axis_equal', 0, 'x_lab', 'Decay', 'y_lab', {'Count'}, ...
            'fig_sz', fig_sz, 'taille_tick', taille, 'taille_stick', taille_axis, ...
            'caption_lgd', caption_lgd) ;

        % ------- Fitting btw model and behavioral data  
        % --------------------------------------------------------------- %
        % (using the best model for each subject)

        % scatter_colormap_lin_fit(data_struct, varargin)
        % 1 entry of struct per subject (to have one color per subject)
        % with use_colormap = 0

        % Plot params
        fig_sz = [1 2 12 12] ; 
        taille = 18 ; taille_axis = 14 ; 
        mk_plot = {'o'} ; sz_pts = 9 ;

        % == * == p1 estimates (compared to IO mean p1)
        fn_p1 = [fn_best_params, 'Fit_post', str_params_fits] ;
        data_struct = struct() ; 
        for idx_subj=1:n_subj
            data_struct(idx_subj).x = [best_IO_p1{idx_subj, :}] ;
            data_struct(idx_subj).y = [est_p1{idx_subj,:}] ;
            data_struct(idx_subj).coeffs = mean_p1_poly ;
        end         
        scatter_colormap_lin_fit(data_struct, 'eps_fig', eps_fig, 'fig_fig', fig_fig, 'pdf_fig', pdf_fig, ...
            'fn_save', fn_p1, 'save_fig',1, 'fig_name', 'Estimations of p1',...
            'add_lin_fit', 1, 'use_colormap', 0, 'axis_equal', 1, 'add_identity', 1, ...
            'x_lab', 'IO proba. estimates', 'y_lab', {'subjective proba. estimates'}, ...
            'fig_sz', fig_sz, 'taille_tick', taille, 'taille_stick', taille_axis, ...
            'sz_pts', sz_pts, 'markers', mk_plot, 'single_fit', 1, ...
            'add_lgd', 0) ; 

        % == * == p1 estimates (compared to IO MAP p1)
        fn_MAPs = [fn_best_params, 'Fit_MAP', str_params_fits] ;         
        data_struct = struct() ; 
        for idx_subj=1:n_subj
            data_struct(idx_subj).x = [best_IO_MAP_p1{idx_subj, :}] ;
            data_struct(idx_subj).y = [est_p1{idx_subj,:}] ;
            data_struct(idx_subj).coeffs = mean_MAP_p1_poly ;
        end
        scatter_colormap_lin_fit(data_struct, 'eps_fig', eps_fig, 'fig_fig', fig_fig, 'pdf_fig', pdf_fig, ...
            'fn_save', fn_MAPs, 'save_fig',1, 'fig_name', 'Estimations of p1 (MAP)',...
            'add_lin_fit', 1, 'use_colormap', 0, 'axis_equal', 1, 'add_identity', 1, ...
            'x_lab', 'IO MAP proba. estimates', 'y_lab', {'subjective proba. estimates'}, ...
            'fig_sz', fig_sz, 'taille_tick', taille, 'taille_stick', taille_axis, ...
            'sz_pts', sz_pts, 'markers', mk_plot, 'single_fit', 1, ...
            'add_lgd', 0) ; 

        % == * == Confidence estimates
        fn_conf = [fn_best_params, 'Fit_conf', str_params_fits] ; 
        data_struct = struct() ; 
        for idx_subj=1:n_subj
            data_struct(idx_subj).x = [best_IO_conf{idx_subj, :}] ;
            data_struct(idx_subj).y = [est_conf{idx_subj,:}] ;
            data_struct(idx_subj).coeffs = mean_conf_poly ;
        end
        scatter_colormap_lin_fit(data_struct, 'eps_fig', eps_fig, 'fig_fig', fig_fig, 'pdf_fig', pdf_fig, ...
            'fn_save', fn_conf, 'save_fig',1, 'fig_name', 'Confidence estimations',...
            'add_lin_fit', 1, 'use_colormap', 0, 'sz_pts', sz_pts, ...
            'x_lab', 'IO confidence', 'y_lab', {'confidence estimates'}, ...
            'fig_sz', fig_sz, 'taille_tick', taille, 'taille_stick', taille_axis,...
            'axis_equal', 0, 'add_identity', 0, 'markers', mk_plot, ...
            'single_fit', 1, 'add_lgd', 0) ; 
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% (2) best parameter from Bayesian model selection + avg
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % == * == .A. Bayesian model comparison (BMC)
    % =================================================================== %
    % ======== IO model
    % =================================================================== %
    
    % Compute quantities required for the BIC
    error_var   = NaN(n_models, n_subj) ; 
    n_data      = NaN(1, n_subj) ; 
    if combine_decay
        if relative_weight
            q_params    = 5 ; % 2 for lin reg + 2 for decay & decayF + Weight
        else
            q_params    = 4 ; % 2 for lin reg + 2 for decay & decayF
        end
    else
        q_params    = 3 ; % 2 for lin reg + 1 for decay
    end

    for idx_subj=1:n_subj
        n_data(idx_subj) = length([est_p1{idx_subj,:}]) ; 
    end       

    for idx_model=1:n_models
        idx_curr_param = cellfun(@(x) strcmpi(names_models{idx_model},x), ...
            learned_param) ; 
        if fit_confidence
            curr_MSE_subj = (squeeze(MSE_all_subj(:,3,idx_curr_param)))' ;
        else
            curr_MSE_subj = (squeeze(MSE_all_subj(:,1,idx_curr_param)))' ;
        end
        % min MSE over all decays (including the flipped ones, if any) 
        error_var(idx_model,:) = min(curr_MSE_subj,[],1) ;             
    end        
    all_BIC_model = repmat(n_data,n_models,1).*log(error_var) + q_params.*...
        log(repmat(n_data,n_models,1)) ; 
    % ensure all BICs are >=0 & as close to 0 as possible
    
    log_model_evidence = -all_BIC_model./2 ; % = marginal likelihood
    options = struct() ; 
    options.DisplayWin = 0 ; options.verbose = 0 ; 
    tS = tic ;
    % cfr https://mbb-team.github.io/VBA-toolbox/wiki/BMS-for-group-studies/
    [posterior,out] = VBA_groupBMC(log_model_evidence, options) ;
    % arg1: Kxn array of log-model evidences (K models; n subjects)
    tEnd = toc(tS) ; disp(['=== Time needed for the model selection: ', num2str(tEnd)])

    model_freqs = out.Ef ; % n_models x 1
    EP = out.ep ; % "Exceedance probability" % 1 x n_models 
    if use_PEP
        EP = (1-out.bor)*EP + out.bor/n_models ; 
    end 
    
    % count the number of participants selecting each model 
    [~, i_opt_mod] = min(all_BIC_model,[],1) ; % minimize the BIC over the models 
    nb_opt_mod = zeros(1,n_models) ; 
    for imod = 1:n_models
        nb_opt_mod(imod) = sum(i_opt_mod==imod) ; 
    end
    
    % names of optimal model per subject (to reload & fit to EEG)
    names_indiv_mod = names_models(i_opt_mod) ;    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save EP and model_freqs
    nrows_tab = 3 ; 
    write_latex_table_latencies({'Model frequencies','Exceedance proba', '#participants'}, ...
        [make_row_vector(model_freqs); make_row_vector(EP); nb_opt_mod], fn_txt_res, ...
        names_models, 'Models', 'Bayesian model selection', 0, ...
        zeros(nrows_tab,n_models), zeros(nrows_tab,n_models),...
        'a','BMS_params') ; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % ---> Index of optimal model
    [~, idx_best_model] = max(EP) ; 
    best_model = names_models{idx_best_model} ; %%%%%%      
    idx_model_BMC = cellfun(@(x) strcmpi(best_model,x), learned_param) ; 

    % ---> Show model proba (MP) + EP
    fn_EP = [fn_best_params, 'MP', str_params_fits] ; 
    x_descr = get_short_names(names_models) ;
    x_descr_phi = x_descr ; 
    x_descr_phi{idx_best_model} = [x_descr_phi{idx_best_model}, ...
        ' - \phi = ', num2str(round(EP(idx_best_model),3))] ; 
    c_bars = parula(100) ; 
    plot_bars(make_col_vector(model_freqs), 'x_lab', '', 'y_lab', 'model probability','fn_save',...
        fn_EP, 'fig_fig', fig_fig, 'pdf_fig',pdf_fig, 'eps_fig', eps_fig, ...
        'x_descr', x_descr_phi, 'c_bars',c_bars, 'use_colormap', 1, 'view_3D', 0, ...
        'add_prior', 1) ; 
    
    % =================================================================== %
    % ======== RW model
    % =================================================================== %
    if remove_AF 
        names_models_RW = names_models_RW(~strcmpi(names_models_RW, 'alternation')) ; 
    end    
    n_models_RW = length(names_models_RW) ;  
    
    fn_txt_RW = [fn_best_RW, 'Opt_params', str_params_fits_RW, '.txt'] ;
    % empty the file fn_txt_res 
    write_empty_lines(fn_txt_RW, 0,'w') ;

    % Compute quantities required for the BIC
    error_var_RW    = NaN(n_models_RW, n_subj) ; 
    q_params_RW     = 3 ; % 2 for lin reg + 1 for learning rate
    
    learned_param_RW = {all_params_RW(:).learned_param} ; 
    all_a_RW = [all_params_RW(:).a_RW] ;    
    unique_a_RW = unique(all_a_RW) ;
    n_a_RW = length(unique_a_RW) ;

    for idx_model=1:n_models_RW
        idx_curr_param = cellfun(@(x) strcmpi(names_models_RW{idx_model},x), ...
            learned_param_RW) ; 
        if fit_confidence
            curr_MSE_subj = (squeeze(MSE_RW(:,2,idx_curr_param)))' ;
        else
            curr_MSE_subj = (squeeze(MSE_RW(:,1,idx_curr_param)))' ;
        end
        % min MSE over all learning rates a_RW 
        error_var_RW(idx_model,:) = min(curr_MSE_subj,[],1) ;             
    end        

    all_BIC_RW = repmat(n_data,n_models_RW,1).*log(error_var_RW) + q_params_RW.*...
        log(repmat(n_data,n_models_RW,1)) ; 
    log_ME_RW = -all_BIC_RW./2 ; % = marginal likelihood
    options = struct() ; options.DisplayWin = 0 ; options.verbose = 0 ; 
    % cfr https://mbb-team.github.io/VBA-toolbox/wiki/BMS-for-group-studies/
    [~,out] = VBA_groupBMC(log_ME_RW, options) ;
    % arg1: Kxn array of log-model evidences (K models; n subjects)

    model_freqs_RW = out.Ef ; % n_models_RW x 1
    EP_RW = out.ep ; % "Exceedance probability" % 1 x n_models_RW 
    if use_PEP
        EP_RW = (1-out.bor)*EP_RW + out.bor/n_models_RW ; 
    end   
    
    % count the number of participants selecting each model 
    [~, i_opt_mod] = min(all_BIC_RW,[],1) ; % minimize the BIC over the models 
    nb_opt_mod = zeros(1,n_models) ; 
    for imod = 1:n_models
        nb_opt_mod(imod) = sum(i_opt_mod==imod) ; 
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save EP and model_freqs
    nrows_tab = 3 ; 
    write_latex_table_latencies({'Model frequencies','Exceedance proba', '#participants'}, ...
        [make_row_vector(model_freqs_RW); make_row_vector(EP_RW); nb_opt_mod], fn_txt_RW, ...
        names_models_RW, 'Models', 'Bayesian model selection', 0, ...
        zeros(nrows_tab,n_models), zeros(nrows_tab,n_models),...
        'a','BMS_params') ; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % ---> Index of optimal model
    [~, idx_best_RW] = max(EP_RW) ; 
    best_model_RW = names_models_RW{idx_best_RW} ; %%%%%%       

    % ---> Show model proba (MP) + EP
    fn_EP_RW = [fn_best_RW, 'MP', str_params_fits_RW] ; 
    x_descr_RW = strcat('d', get_short_names(names_models_RW)) ;% dIF, dAF, dTP
    
    x_descr_RW_phi = x_descr_RW ; 
    x_descr_RW_phi{idx_best_RW} = [x_descr_RW{idx_best_RW}, ...
        ' - \phi = ', num2str(round(EP_RW(idx_best_RW),3))] ; 
    c_bars = parula(100) ; 
    plot_bars(make_col_vector(model_freqs_RW), 'x_lab', '', 'y_lab', 'model probability','fn_save',...
        fn_EP_RW, 'fig_fig', fig_fig, 'pdf_fig',pdf_fig, 'eps_fig', eps_fig, ...
        'x_descr', x_descr_RW_phi, 'c_bars',c_bars, 'use_colormap', 1, 'view_3D', 0, ...
        'add_prior', 1) ;     
    
    % =================================================================== %
    % ======== IO + RW models
    % =================================================================== %
    
    fn_txt_both = [fn_best_params, 'Opt_params_RW', str_fits_IO_RW, '.txt'] ;
    % empty the file fn_txt_res 
    write_empty_lines(fn_txt_both, 0,'w') ;
    names_models_both = [x_descr,x_descr_RW] ;     
    log_ME_both = [log_model_evidence; log_ME_RW] ;
    
    [~,out] = VBA_groupBMC(log_ME_both, options) ;
    % arg1: Kxn array of log-model evidences (K models; n subjects)

    model_freqs_both = out.Ef ; % n_models+n_models_RW x 1
    EP_both = out.ep ; % "Exceedance probability" % 1 x n_models_RW 
    if use_PEP
        EP_both = (1-out.bor)*EP_both + out.bor/(n_models+n_models_RW) ; 
    end
    
    % count the number of participants selecting each model 
    [~, i_opt_mod] = max(log_ME_both,[],1) ; % max model evidence over models 
    nb_opt_mod = zeros(1,n_models+n_models_RW) ; 
    for imod = 1:(n_models+n_models_RW)
        nb_opt_mod(imod) = sum(i_opt_mod==imod) ; 
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save EP and model_freqs
    nrows_tab = 3 ; 
    write_latex_table_latencies({'Model frequencies','Exceedance proba', '#participants'}, ...
        [make_row_vector(model_freqs_both); make_row_vector(EP_both); nb_opt_mod], fn_txt_both, ...
        names_models_both, 'Models', 'Bayesian model selection', 0, ...
        zeros(nrows_tab,n_models+n_models_RW), zeros(nrows_tab,n_models+n_models_RW),...
        'a','BMS_params') ; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % display individual model probabilities (model evidence normalized by
    % subject)    
    indiv_models = (exp(log_ME_both))' ; % (n_subj, n_models+n_models_RW) 
    
    indiv_models = indiv_models./(repmat(sum(indiv_models, 2), ...
        1, n_models+n_models_RW)) ; % (n_subj, n_models+n_models_RW)
    fn_indiv = [fn_best_params, 'MPmap', str_fits_IO_RW] ; 
    if remove_AF
        fs_MP = [5 10 10 12] ;
    else
        fs_MP = [5 10 15 12] ; 
    end
    plot_heatmap('individual model proba', 1:(n_models+n_models_RW), ...
        1:n_subj, indiv_models, flipud(gray), ['proba'], 'subject', ...
        '', 'filename', fn_indiv, 'eps_fig', eps_fig, 'fig_fig', fig_fig, 'pdf_fig', pdf_fig, ...
        'xticks', 1:(n_models+n_models_RW), ... % , 'yticks', 1:n_subj
        'xtick_labs', names_models_both, 'clim_val',[0,1], 'fig_size', fs_MP) ; 
    % flipud(gray), autumn, summer 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % ---> Index of optimal model
    [~, idx_best_both] = max(EP_both) ; 
    best_model_RW_IO = names_models_both{idx_best_both} ; %%%%%%      

    % ---> Show model proba (MP) + EP
    fn_EP_both = [fn_best_params, 'MP', str_fits_IO_RW] ; 
    x_descr_both = names_models_both ; 
    
    x_descr_both{idx_best_both} = [x_descr_both{idx_best_both}, ...
        ' - \phi = ', num2str(round(EP_both(idx_best_both),3))] ; 
    c_bars = parula(100) ; 
    plot_bars(make_col_vector(model_freqs_both), 'x_lab', '', 'y_lab', 'model probability','fn_save',...
        fn_EP_both, 'fig_fig', fig_fig, 'pdf_fig',pdf_fig, 'eps_fig', eps_fig, ...
        'x_descr', x_descr_both, 'c_bars',c_bars, 'use_colormap', 1, 'view_3D', 0, ...
        'add_prior', 1, 'fig_sz',[1 2 22 12]) ;   

    % == * == .B. Bayesian model averaging (BMA)
    % =================================================================== %  
    
    % =================================================================== %
    % ======== RW models: posterior of learning rate a_RW
    % =================================================================== %
    post_a_RW   = zeros(n_a_RW, n_subj) ; % posterior of a_RW
    
    if fit_confidence
        all_BICs = repmat(n_data,n_params_RW,1).*log( (squeeze(MSE_RW(:,2,:)))' ) + ...
            q_params_RW.*log(repmat(n_data,n_params_RW,1)) ; 
    else
        all_BICs = repmat(n_data,n_params_RW,1).*log( (squeeze(MSE_RW(:,1,:)))' ) + ...
            q_params_RW.*log(repmat(n_data,n_params_RW,1)) ; 
    end
    % Ensure BICs are >=0 and as close to 0 as possible (for exponential)
    all_BICs = all_BICs - min(all_BICs(:)) ;
    all_model_evidence = exp(-all_BICs./2) ;% model evidences
    
    for idx_a = 1:n_a_RW
        idx_curr_a = unique_a_RW(idx_a)==all_a_RW ; 
        if remove_AF
            idx_alt = cellfun(@(x) strcmpi('alternation',x), learned_param_RW) ;
            idx_curr_a = idx_curr_a & ~idx_alt ; 
        end

        post_a_RW(idx_a,:) = sum(...
            all_model_evidence(idx_curr_a,:),1) ; % idx_curr_a
    end
    
    norm_for_post = repmat(nansum(post_a_RW,1), n_a_RW,1) ; 
    post_a_RW = post_a_RW./norm_for_post ;
    fn_post_a = [fn_best_RW, 'Post_a', str_params_fits_RW] ;
    
    % ---> Index of optimal learning rate
    mean_post_a = nanmean(post_a_RW,2) ;
    [best_post_a, idx_best_a] = max(mean_post_a) ; 
    best_a_RW = unique_a_RW(idx_best_a) ; 
    str_LR = ['Optimal LR for the RW models: ', num2str(best_a_RW)] ;
    %disp(str_LR) ;
    write_empty_lines(fn_txt_both, 2,'a') ; write_to_txt(fn_txt_both, str_LR,'a',2) ;
    write_empty_lines(fn_txt_RW, 2,'a') ; write_to_txt(fn_txt_RW, str_LR,'a',2) ;

    % ---> Show posterior proba of LR (mean pm std)
    x_to_plot = make_row_vector(unique_a_RW) ; 
    y_to_plot = make_row_vector(mean_post_a) ; 
    y_std_to_plot = (nanstd(post_a_RW,[],2))' ; 
    col_mean = [0,0,0] ; 
    col_std = 0.6*ones(1,3) ; 
    caption_lgd = {} ; 
    if show_sem
        y_std_to_plot = y_std_to_plot./sqrt(n_subj) ;
    end
    figs_post = [1 2 14 12] ; 
    plot_mean_std(x_to_plot, y_to_plot, y_std_to_plot, 'col_mean', ...
        col_mean, 'x_lab', 'learning rate', 'y_lab', 'probability','fn_save',...
        fn_post_a, 'fig_fig', fig_fig, 'pdf_fig',pdf_fig, 'eps_fig', eps_fig, ...
        'add_prior', 1, 'shaded_std', 1, 'show_max', 1, 'col_std', col_std, ...
        'caption_lgd', caption_lgd, 'logx_curve', 1, 'fig_sz', figs_post) ;
    
    % =================================================================== %
    % ======== IO models: for each decay (flipped or not) 
    % =================================================================== %

    if combine_decay
        % posterior of decays
        post_decays   = zeros(n_decayFs, n_decays, n_subj) ; 
    else
        post_decays   = zeros(n_decays, n_subj) ; % posterior of decays
        post_flipped_decays = zeros(n_decayFs, n_subj) ;
    end

    if fit_confidence
        all_BICs = repmat(n_data,n_params,1).*log( (squeeze(MSE_all_subj(:,3,:)))' ) + ...
            q_params.*log(repmat(n_data,n_params,1)) ; 
    else
        all_BICs = repmat(n_data,n_params,1).*log( (squeeze(MSE_all_subj(:,1,:)))' ) + ...
            q_params.*log(repmat(n_data,n_params,1)) ; 
    end
    all_model_evidence = exp(-all_BICs./2) ; 

    idx_alt = cellfun(@(x) strcmpi('alternation',x), learned_param) ;
    if combine_decay          
        for idx_decay = 1:n_decays
            idx_curr_decay = unique_decays(idx_decay)==all_decays ; 
            if remove_AF
                idx_curr_decay = idx_curr_decay & ~idx_alt ; 
            end
            for idx_decayF = 1:n_decayFs
                idx_curr_decayF = unique_decayFs(idx_decayF)==all_decayFs ; 

                post_decays(idx_decayF, idx_decay, :) = ...
                    sum(all_model_evidence(idx_curr_decay & idx_curr_decayF,:),1) ;
            end
        end
        norm_for_post = repmat( (nansum(nansum(post_decays, 1),2)), n_decays, ...
            n_decayFs, 1) ;
    else
        % backward decays:
        for idx_decay = 1:n_decays
            idx_curr_decay = unique_decays(idx_decay)==all_decays ; 
            if remove_AF
                idx_curr_decay = idx_curr_decay & ~idx_alt ; 
            end

            post_decays(idx_decay,:) = sum(...
                all_model_evidence(idx_curr_decay & all_leaky_int,:),1) ; 
        end
        % forward decays: 
        for idx_decayF = 1:n_decayFs
            idx_curr_decayF = unique_decayFs(idx_decayF)==all_decayFs ; 
            if remove_AF
                idx_curr_decayF = idx_curr_decayF & ~idx_alt ; 
            end

            post_flipped_decays(idx_decayF,:) = sum(...
                all_model_evidence(idx_curr_decayF & all_leakyF_int,:),1) ; 
        end
        % n_decays, n_subj
        % attention: for the normalization, don't count perfect
        % integration twice (if encoded both with leaky_int &
        % leakyF_int)
        norm_for_post = repmat( nansum(post_decays,1)+...
            nansum(post_flipped_decays(setdiff(1:n_decayFs, idx_decayF_inf),:),1), ...
            n_decays,1) ; 
        if consider_flip
            post_flipped_decays = post_flipped_decays./norm_for_post ;
        end
    end
    post_decays = post_decays./norm_for_post ; 

    fn_post_w = [fn_best_params, 'Post_w', str_params_fits] ;
    fn_post_w_insets = [fn_best_params, 'Post_w_in', str_params_fits] ;
    if combine_decay
        % post_decays: n_decayFs, n_decays, n_subj
        mean_post_decays = squeeze(nanmean(post_decays,3)) ; 
        [best_post, lin_idx_best] = max(mean_post_decays(:)) ; 
        [idx_best_decayF, idx_best_decay] = ind2sub([n_decayFs, n_decays], lin_idx_best) ; 

        best_decay = unique_decays(idx_best_decay) ;            
        best_decayF = unique_decayFs(idx_best_decayF) ;  

        idx_decay_BMA = best_decay==all_decays & best_decayF==all_decayFs ;
        idx_int_BMA = all_leaky_int & all_leakyF_int ; 

        % ---> Show posterior proba of decays (heatmap)
        x_to_plot = make_row_vector(unique_decays) ; 
        y_to_plot = make_row_vector(unique_decayFs) ; 
        data_to_plot = mean_post_decays ; 
        eps_data = 0.00001 ; % avoid -Inf values with log
        if log_data
           data_to_plot = log(data_to_plot+eps_data) ;  
        end

        % ===*=== Interpolate data if necessary
        if interp_and_log       
            mX = min(x_to_plot) ; mY = min(y_to_plot) ; 
            [xmesh, ymesh] = meshgrid(x_to_plot, y_to_plot) ;
            % /!\ Assuming integers values for x_data & y_data
            x_to_plot = mX:max(x_to_plot) ; 
            y_to_plot = mY:max(y_to_plot) ;
            [xmesh_plot, ymesh_plot] = meshgrid(x_to_plot, y_to_plot) ; 

            data_to_plot = interp2(xmesh, ymesh, data_to_plot, xmesh_plot, ymesh_plot) ;
            log_xscale = 1 ; log_yscale = 1 ; 
            tick_show = [1, 2, 5, 10, 20, 40, 75, 150] ; 
        else
            log_xscale = 0 ; log_yscale = 0 ; 
            tick_show = [1,50, 100, 150] ; 
        end

        max_post = max(data_to_plot(:)) ; min_post = min(data_to_plot(:)) ; 
        % accomodate negative min values (using range)
        [f, power10] = get_first_signif_digit(max_post-min_post) ;
        if f<4
            step_levels = 5*10^(power10-1) ; % 5
        else
            step_levels = 2*10^power10 ; 
        end

        start_level = floor(min_post/step_levels)*step_levels ; 
        levels_proba = start_level:step_levels:(max_post+0.9*step_levels) ; 
        n_levels = length(levels_proba) ; 
        % levels for the proba to plot if map_opt in {2,3}

        if log_data
            hc_ticklabels = round(exp(levels_proba)-eps_data,2, 'significant') ;
            hc_ticklabels(hc_ticklabels<=10^(-8)) = 0 ; 
        else        
            hc_ticklabels = levels_proba ; %(1:2:end) ; % could be strings if in radians...
        end
        hc_ticks = levels_proba ; 
        clim_val = [levels_proba(1), levels_proba(end)] ;

        switch map_opt
            case 1
                cmap_used = parula ; 
            case 2
                cmap_used = magma(n_levels-1) ; 
            case 3
                cmap_used = magma(n_levels) ; 
                % gap_color: data range spanned by one color
                gap_color = (levels_proba(end)-levels_proba(1))./n_levels ; 
                %hc_ticks = levels_proba(1)+0.5*gap_color:2*gap_color:(levels_proba(end)-0.5*gap_color) ;         
                %hc_ticklabels = hc_ticklabels(1:2:end) ;

                hc_ticks = levels_proba(1)+0.5*gap_color:gap_color:(levels_proba(end)-0.5*gap_color) ; %0.5*gap_color:gap_color:(levels_proba(end)-0.5*gap_color) ;         
                hc_ticklabels = hc_ticklabels(1:end) ;
        end

        % Define ticks and ticklabels
        n_x = length(x_to_plot) ; n_y = length(y_to_plot) ;
        xticks_map = x_to_plot ; yticks_map = y_to_plot ; 
        xtick_labs = cell(n_x, 1) ; xtick_labs(:) = {''} ; 
        ytick_labs = cell(n_x, 1) ; ytick_labs(:) = {''} ; 
        for idx_x=1:n_x
            if any(tick_show==x_to_plot(idx_x))
                xtick_labs{idx_x} = num2str(x_to_plot(idx_x)) ; 
            end
        end
        for idx_y=1:n_y
            if any(tick_show==y_to_plot(idx_y))
                ytick_labs{idx_y} = num2str(y_to_plot(idx_y)) ; 
            end
        end
        if interp_and_log
            % Keep the labels but remove the 0's for the log axes
            if mX==0
                x_to_plot = x_to_plot + 1 ; xticks_map = xticks_map + 1 ; 
            end
            if mY==0
                y_to_plot = y_to_plot + 1 ; yticks_map = yticks_map + 1 ; 
            end
        end
        plot_heatmap('Posterior proba of fwd & bwd leaks', x_to_plot, ...
            y_to_plot, data_to_plot, cmap_used, ['Probability'], 'Decay-F', ...
            'Decay-B', 'log_xscale', log_xscale, 'log_yscale', log_yscale, ...
            'filename', [fn_post_w, str_map], 'levels', levels_proba, ...
            'map_opt', map_opt, 'hc_ticklabels',hc_ticklabels, 'hc_ticks', hc_ticks, ...
            'eps_fig', eps_fig, 'fig_fig', fig_fig, 'pdf_fig', pdf_fig, ...
            'xticks', xticks_map, 'yticks', yticks_map, 'xtick_labs', xtick_labs, ...
            'ytick_labs', ytick_labs, 'clim_val',clim_val) ; 

        % --- Data for the insets --- (when leaks are = 0)
        % (1) Right inset: when leak = 0. 
        x_vec_R = make_row_vector(unique_decayFs) ; 
        idx_dF = x_vec_R>0 ; 
        x_vec_R = x_vec_R(idx_dF) ;  % avoiding considering the case (0,0)
        mean_R = (mean_post_decays(idx_dF,1))' ;
        std_post_decays = squeeze(nanstd(post_decays,[], 3)) ;
        std_R = (std_post_decays(idx_dF, 1))' ; 
        if show_sem
            std_R = std_R./sqrt(n_subj) ; 
        end
        % (1) Top inset: when leakF = 0. 
        x_vec_T = make_row_vector(unique_decays) ; 
        idx_d = x_vec_T>0 ; x_vec_T = x_vec_T(idx_d) ; 
        mean_T = mean_post_decays(1, idx_d) ; 
        std_T = std_post_decays(1, idx_d) ; 
        if show_sem
            std_T = std_T./sqrt(n_subj) ; 
        end
        if interp_and_log
            % Remove the 0's for the log axes
            if mX==0 ; x_vec_T = x_vec_T + 1 ; end
            if mY==0 ; x_vec_R = x_vec_R + 1 ; end
        end

        plot_heatmap_insets('Posterior proba of fwd & bwd leaks', x_to_plot, ...
            y_to_plot, data_to_plot, cmap_used, ['Probability'], 'Decay-F', ...
            'Decay-B', x_vec_R, mean_R, std_R, x_vec_T, mean_T, std_T, ...
            'log_xscale', log_xscale, 'log_yscale', log_yscale, ...
            'filename', [fn_post_w_insets, str_map], 'levels', levels_proba, ...
            'map_opt', map_opt, 'hc_ticklabels',hc_ticklabels, 'hc_ticks', hc_ticks, ...
            'xticks', xticks_map, 'yticks', yticks_map, 'xtick_labs', xtick_labs, ...
            'ytick_labs', ytick_labs, 'clim_val', clim_val, 'add_prior', ...
            1, 'show_max', 1, 'prior_val', 1/numel(mean_post_decays), ...
            'eps_fig', eps_fig, 'fig_fig', fig_fig, 'pdf_fig', pdf_fig) ;
        %%%%%%%%%%%%%%%%%%%%%%%%
        % (i) add inset above (along x-axis): showing post of decayB
        %       when decayF = Inf (mean pm sem). 
        % (ii) add vertical inset on the right showing post of decayF
        %       when decayB = Inf (mean pm sem). 
        %%%%%%%%%%%%%%%%%%%%
    else
        % ---> Index of optimal decay
        mean_post_decays = nanmean(post_decays,2) ; 
        [best_post, idx_best_decay] = max(mean_post_decays) ; 
        best_decay = unique_decays(idx_best_decay) ;        

        % ---> Index of optimal flipped decay (decayF)
        mean_post_flipped_decays = nanmean(post_flipped_decays,2) ; 
        [best_flipped_post, idx_best_flipped_decay] = max(mean_post_flipped_decays) ; 
        best_decayF = unique_decays(idx_best_flipped_decay) ;

        if isempty(best_flipped_post) || best_post>=best_flipped_post
            idx_decay_BMA = best_decay==all_decays ;
            idx_int_BMA = all_leaky_int ;
        else
            idx_decay_BMA = best_decayF==all_decayFs ;
            idx_int_BMA = all_leakyF_int ;
        end

        % ---> Show posterior proba of decays (mean pm std)
        if consider_flip
            x_to_plot = [make_row_vector(unique_decays); ...
                make_row_vector(unique_decayFs)] ; 
            y_to_plot = [make_row_vector(mean_post_decays); ...
                make_row_vector(mean_post_flipped_decays)] ; 
            y_std_to_plot = [(nanstd(post_decays,[],2))'; ...
                (nanstd(post_flipped_decays,[],2))'] ;            
            col_mean = [0, 128, 255; 102, 204, 0]./255 ; 
            col_std = col_mean ; 
            caption_lgd = {'Bwd', 'Fwd'} ;
        else
            x_to_plot = make_row_vector(unique_decays) ; 
            y_to_plot = make_row_vector(mean_post_decays) ; 
            y_std_to_plot = (nanstd(post_decays,[],2))' ; 
            col_mean = [0,0,0] ; 
            col_std = 0.6*ones(1,3) ; 
            caption_lgd = {} ; 
        end

        if show_sem
            y_std_to_plot = y_std_to_plot./sqrt(n_subj) ;
        end
        plot_mean_std(x_to_plot, y_to_plot, y_std_to_plot, 'col_mean', ...
            col_mean, 'x_lab', 'decay', 'y_lab', 'probability','fn_save',...
            fn_post_w, 'fig_fig', fig_fig, 'pdf_fig',pdf_fig, 'eps_fig', eps_fig, ...
            'add_prior', 1, 'shaded_std', 1, 'show_max', 1, 'col_std', col_std, ...
            'caption_lgd', caption_lgd, 'fig_sz', figs_post) ;
    end        
    str_decay = ['Optimal decay for the IO models: ', num2str(best_decay)] ;
    str_decayF = ['Optimal decayF for the IO models: ', num2str(best_decayF)] ;
    write_to_txt(fn_txt_both, str_decay,'a',2) ; write_to_txt(fn_txt_both, str_decayF,'a',2) ;
    
    
    % == * == save optimal model for each participant (to reload & fit to EEG)    
    % optimal decay per subject
    [~, idx_indiv_dec] = max(post_decays,[],1) ;% (n_decays, n_subj)
    indiv_decays = unique_decays_true(idx_indiv_dec) ;
    names_indiv_mod = names_indiv_mod ;
    
    %best_vs_indiv_decays = [best_decays', indiv_decays'] % slight diff may appear because sum over models or not
    
    fn_indiv_models = [fn_dat, 'opt_models_',num2str(n_subj),'subj.mat'] ; 
    save(fn_indiv_models, '-v7.3','names_indiv_mod', 'best_decay', 'indiv_decays') ;
    % can always use the best avg decay, or select the optimal one per
    % participant
    
    
    % == * == .B-b. Bayesian model averaging (BMA)
    % for each weight, if used as parameter
    % === === === === === === === === === === === === === === === === %
    if combine_decay && relative_weight
        all_weights = [all_params(:).weight] ;
        unique_weights = unique(all_weights) ; 
        n_weights = length(unique_weights) ;

        % Posterior for 'Weight':
        post_weight = zeros(n_weights, n_subj) ;
        for idx_weight = 1:n_weights
            idx_curr_weight = unique_weights(idx_weight)==all_weights ; 

            post_weight(idx_weight,:) = sum(...
                all_model_evidence(idx_curr_weight,:),1) ; 
        end
        norm_for_post_w = repmat( (nansum(post_weight, 1)), n_weights, 1) ;
        post_weight = post_weight./norm_for_post_w ; 

        mean_post_w = nanmean(post_weight,2) ; 
        [best_post_w, idx_best_w] = max(mean_post_w) ; 
        best_weight = unique_weights(idx_best_w) ;   

        idx_decay_BMA = idx_decay_BMA & best_weight==all_weights ;    

        % ---> Show posterior proba of weights (mean pm std)
        x_to_plot = make_row_vector(unique_weights) ; 
        y_to_plot = make_row_vector(mean_post_w) ; 
        y_std_to_plot = (nanstd(post_weight,[],2))' ; 
        col_mean = [0,0,0] ; 
        col_std = 0.6*ones(1,3) ; 
        caption_lgd = {} ; 
        if show_sem
            y_std_to_plot = y_std_to_plot./sqrt(n_subj) ;
        end
        % weight: A (>< decay w)
        fn_post_A = [fn_best_params, 'Post_A', str_params_fits] ;
        plot_mean_std(x_to_plot, y_to_plot, y_std_to_plot, 'col_mean', ...
            col_mean, 'x_lab', 'Bwd Weighting', 'y_lab', 'Probability','fn_save',...
            fn_post_A, 'fig_fig', fig_fig, 'pdf_fig',pdf_fig, 'eps_fig', eps_fig, ...
            'add_prior', 1, 'shaded_std', 1, 'show_max', 1, 'col_std', col_std, ...
            'caption_lgd', caption_lgd, 'logx_curve', 0) ;
    end
end


end

function all_params = get_all_params_to_test(consider_flip, combine_decay, ...
    relative_weight)


all_AboutFirst  = {'WithoutFirst'} ; 
% since our sequences have length 100, decay should not be larger than
% 100...
%all_decays      =  [1:20] ; 
% Maheu2019: 55 values tested (cfr fig 5E & p18)

%all_decays      =  [1:20, 22:2:200] ;
all_weights = [1] ; 

if combine_decay
    % include 0 to also have the cases of only F or B integration!
    if relative_weight
        all_decays  = [0:20, Inf] ; % 22 values
        all_weights = [0.5, 1, 2, 3] ; % weights for the bwd leaks (compared to fwd leaks)
    else
        all_decays      = [0:20, 22:2:40, 45:5:100, Inf] ; % 44 values
        all_decays      = [0:20, 22:2:40, 45:5:90, Inf] ; % 42 values
        % 42*42*3 = 5292
    end
elseif consider_flip
    all_decays      = [1:20, 22:2:40, 45:5:150, Inf] ; % 53 values
else
    all_decays      = [1:20, 22:2:40, 45:5:200, Inf] ; % 63 values
    all_decays      = [1:20, 22:2:40, 45:5:250, Inf] ; % 73 values
    all_decays      = [1:20, 22:2:40, 45:5:300, Inf] ; % 83 values
    all_decays      = [1:20, 22:2:40, 45:5:350, Inf] ; % 93 values
    all_decays      = [1:20, 22:2:40, 45:5:400, Inf] ; % 103 values
    % 63*3 = 189; 
    % 73*3 = 219 ; 83*3 = 249; 93*3 = 279; 103*3 = 309
end

all_windows     = [] ;
all_learned_param = {'transition', 'frequency', 'alternation'} ; 
% 'frequency', 'transition', 'alternation' (cfr [NEWW] in IdealObserver)
if consider_flip
    all_flips = [0,1] ; 
else
    all_flips = 0 ; 
end
n_flips = length(all_flips) ; 
n_first = length(all_AboutFirst) ; 
n_decays = length(all_decays) ; 
n_weights = length(all_weights) ; 
n_windows = length(all_windows) ; 
n_learned_params = length(all_learned_param) ; 

if combine_decay
    all_windows     = [] ; %all_windowsF = [] ; 
    all_flips       = 0 ;
    
    all_decaysF     = all_decays ;
    n_params = n_first*n_learned_params*n_weights*(n_decays^2) ; 
else
    all_decaysF = [0] ; 
    n_params = n_first*n_learned_params*(n_flips*(n_windows + n_decays)) ; 
end

n_flips = length(all_flips) ; 

all_params = repmat(struct('AboutFirst','WithoutFirst', 'learned_param', 'transition',...
    'decay',[Inf],'decayF', [Inf], 'window', [Inf], 'windowF', [Inf], ...
    'leaky_int', 0, 'leakyF_int', 0, 'window_int',0, 'windowF_int', 0, ...
    'MemParam', [], 'weight', 1, 'str_model', '_TP_WithoutFirst'), n_params,1) ;  
% Preallocate the array which will contain all the sets of params


idx_param = 1 ; 

for idx_first = 1:n_first    
    AboutFirst = all_AboutFirst{idx_first} ; 
    for idx_learned_param = 1:n_learned_params
        learned_param = all_learned_param{idx_learned_param} ;      
        for weight = all_weights
            for idx_flip=1:n_flips
                curr_flip = all_flips(idx_flip) ; 

                % ===> Leaky integration
                for decay = all_decays
                    for decayF=all_decaysF   
                        all_params(idx_param).AboutFirst = AboutFirst ;
                        all_params(idx_param).learned_param = learned_param ;
                        all_params(idx_param).weight = weight ;

                        MemParam = {} ;  

                        if combine_decay
                            dB = decay ; dF = decayF ; 
                            leaky_bool = 1 ; leakyF_bool = 1 ; 
                        else
                            % can have all_flips = [0,1] ; when ~combine_decay
                            if curr_flip
                                % use decay (NOT decayF) for the flipped decay!!
                                dB = 0 ; dF = decay ; 
                                leaky_bool = 0 ; leakyF_bool = 1 ;                            
                            else
                                dB = decay ; dF = 0 ; 
                                leaky_bool = 1 ; leakyF_bool = 0 ;
                            end                        
                        end

                        all_params(idx_param).decay = dB ;
                        all_params(idx_param).decayF = dF ; 
                        all_params(idx_param).leaky_int = leaky_bool ;
                        all_params(idx_param).leakyF_int = leakyF_bool ;                      
                        MemParam(end+1:end+2) = {'Decay', dB} ;
                        MemParam(end+1:end+2) = {'DecayF', dF} ;
                        if weight~=1 ; MemParam(end+1:end+2) = {'Weight',weight} ; end
                        all_params(idx_param).MemParam = MemParam ;

                        str_model = get_model_param_str(AboutFirst, learned_param, ...
                            dB, dF, Inf, Inf, leaky_bool, leakyF_bool, 0, 0, weight) ; 
                        all_params(idx_param).str_model = str_model ; 

                        idx_param = idx_param + 1 ; 
                    end
                end

                % ===> Windowed integration
                for window = all_windows
                    all_params(idx_param).window = window ;
                    all_params(idx_param).AboutFirst = AboutFirst ;
                    all_params(idx_param).learned_param = learned_param ; 
                    all_params(idx_param).window_int = 1 ;
                    %all_params(idx_param).weight = weight ; % combined
                    %windows not used

                    MemParam = {} ;
                    if curr_flip
                        % (AboutFirst, learned_param, decay, decayF, window, 
                        %  windowF, leaky_int, leakyF_int, window_int, windowF_int)
                        str_model = get_model_param_str(AboutFirst, learned_param, ...
                            Inf, Inf, Inf, window, 0, 0, 0, 1) ; 

                        MemParam(end+1:end+2) = {'LimitedF', window} ; 
                    else
                        str_model = get_model_param_str(AboutFirst, learned_param, ...
                            Inf, Inf, window, Inf, 0, 0, 1, 0) ; 
                        MemParam(end+1:end+2) = {'Limited', window} ; 
                    end  
                    %if weight~=1 ; MemParam(end+1:end+2) = {'Weight',weight} ; end
                    all_params(idx_param).str_model = str_model ; 
                    all_params(idx_param).MemParam = MemParam ;

                    idx_param = idx_param + 1 ; 
                end

            end
        end
    end
end


end

function all_params_RW = get_all_params_RW()

% learning rates, in [0, 1]
all_a = 0.05:0.05:1 ; % 20 values; 3*21: 63
all_a = [0.005:0.005:0.04, 0.05:0.05:(1-0.05)] ; % 27 values
all_a = [0.005:0.01:(1-0.05)] ; % 95 values
all_a = [0.005:0.005:0.04, 0.05:0.01:(1-0.05)] ; % 99 values
n_a = length(all_a) ; 
all_learned_param = {'transition', 'frequency', 'alternation'} ; 
n_learned_params = length(all_learned_param) ;

n_params = n_a*n_learned_params ; 


all_params_RW = repmat(struct('learned_param', 'transition',...
    'a_RW', 0.1, 'str_model', '_dTP_a10'), n_params,1) ;  
% Preallocate the array which will contain all the sets of params


idx_param = 1 ; 


for idx_learned_param = 1:n_learned_params
    learned_param = all_learned_param{idx_learned_param} ;      
    for a_RW = all_a
        all_params_RW(idx_param).learned_param = learned_param ;
        all_params_RW(idx_param).a_RW = a_RW ;
        
        str_model = get_model_param_str('', learned_param,0,0, 0, 0, 0, 0, ...
            0, 0, 0,[],'RW',a_RW) ; 
        all_params_RW(idx_param).str_model = str_model ; 

        idx_param = idx_param + 1 ;
    end
end


end



function all_params = get_all_extreme_params_to_test()
% Extreme params when combine_decay (leaks of 0 or Inf). 


all_AboutFirst  = {'WithoutFirst'} ; 
% since our sequences have length 100, decay should not be larger than
% 100...
all_decays      = [0:20, 22:2:40, 45:5:150, Inf] ; % 54 values
all_extremes    = [0, Inf] ; 

all_windows     = [] ;
all_learned_param = {'transition', 'frequency', 'alternation'} ; 
% 'frequency', 'transition', 'alternation' (cfr [NEWW] in IdealObserver)

all_flips = 0 ; n_flips = length(all_flips) ; 

n_first = length(all_AboutFirst) ; 
n_decays = length(all_decays) ; n_extremes = length(all_extremes) ;
n_windows = length(all_windows) ; 
n_learned_params = length(all_learned_param) ; 

n_params = n_first*n_learned_params*(2*n_decays*n_extremes) ; 

all_params = repmat(struct('AboutFirst','WithoutFirst', 'learned_param', 'transition',...
    'decay',[Inf],'decayF', [Inf], 'window', [Inf], 'windowF', [Inf], ...
    'leaky_int', 0, 'leakyF_int', 0, 'window_int',0, 'windowF_int', 0, ...
    'MemParam', [], 'str_model', '_TP_WithoutFirst'), n_params,1) ;  
% Preallocate the array which will contain all the sets of params


idx_param = 1 ; 

for idx_first = 1:n_first    
    AboutFirst = all_AboutFirst{idx_first} ; 
    for idx_learned_param = 1:n_learned_params
        learned_param = all_learned_param{idx_learned_param} ;         
            
        % ===> Leaky-B integration with extremes decayF
        for decay = all_decays
            for decayF = all_extremes   
                all_params(idx_param).AboutFirst = AboutFirst ;
                all_params(idx_param).learned_param = learned_param ;

                all_params(idx_param).decay = decay ;
                all_params(idx_param).decayF = decayF ; 
                all_params(idx_param).leaky_int = 1 ;
                all_params(idx_param).leakyF_int = 1 ;  
                MemParam = {} ; 
                MemParam(end+1:end+2) = {'Decay', decay} ;
                MemParam(end+1:end+2) = {'DecayF', decayF} ;
                all_params(idx_param).MemParam = MemParam ;

                str_model = get_model_param_str(AboutFirst, learned_param, ...
                    decay, decayF, Inf, Inf, 1, 1, 0, 0) ; 

                all_params(idx_param).str_model = str_model ; 

                idx_param = idx_param + 1 ; 
            end
        end
        
        % ===> Leaky-F integration with extremes decay
        for decayF = all_decays
            for decay = all_extremes   
                all_params(idx_param).AboutFirst = AboutFirst ;
                all_params(idx_param).learned_param = learned_param ;

                all_params(idx_param).decay = decay ;
                all_params(idx_param).decayF = decayF ; 
                all_params(idx_param).leaky_int = 1 ;
                all_params(idx_param).leakyF_int = 1 ;  
                MemParam = {} ; 
                MemParam(end+1:end+2) = {'Decay', decay} ;
                MemParam(end+1:end+2) = {'DecayF', decayF} ;
                all_params(idx_param).MemParam = MemParam ;

                str_model = get_model_param_str(AboutFirst, learned_param, ...
                    decay, decayF, Inf, Inf, 1, 1, 0, 0) ; 

                all_params(idx_param).str_model = str_model ; 

                idx_param = idx_param + 1 ; 
            end
        end
            
    end
end


end




