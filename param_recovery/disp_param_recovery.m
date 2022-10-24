% Display the results of the parameter recovery analysis
% - correlation between true and fitted parameters
% - percentages of models recovered


% GOAL: 
% [param recov] have one correlation plot for each model (6 models)
% [model recov] histogram of % of each recovered model, per model (eg: true
% model = TP --> hist of the 3 models proportions recovered). 


model_type      = 'IO' ;    % RW (Rescorla-Wagner = delta rule), IO (Bayesian)
% RW: doesn't make sense as we can't sample proba estimates... recovery
% will be perfect

eps_fig = 0 ; fig_fig = 0 ; pdf_fig = 0 ;  
if strcmpi(model_type, 'RW') ; log_axes = 0 ; else ; log_axes = 1 ; end 

% ======================================================================= %
% ===*=== filenames etc
% ======================================================================= %
fn_save = ['./data_simu/'] ;
fn_plot = ['./figures/'] ; 
if ~exist(fn_plot, 'dir') ; mkdir(fn_plot) ; end
if strcmpi(model_type, 'RW')
    all_fn_simu = {'RW_dAF_10rep_21w', 'RW_dIF_10rep_20w', 'RW_dTP_10rep_21w'} ; 
    % ndata to compute BIC for model recovery
    all_fn_fit = strcat(all_fn_simu, '_297fits') ;     
else
    all_fn_simu = {'IO_AF_30rep_14w', 'IO_IF_30rep_14w', 'IO_TP_30rep_17w'} ;
    all_fn_fit = strcat(all_fn_simu, '_309fits') ;     
end
models_sim = {'alternation', 'frequency', 'transition'} ; % simulated models

n_models_sim = length(all_fn_simu) ; 

% ======================================================================= %
% ===*=== (1) Parameter recovery
% ======================================================================= %

for imod = 1:n_models_sim   
    
    fit_data = load([fn_save, all_fn_fit{imod}, '.mat']) ; 
    all_MSE = fit_data.all_MSE ; % NaN(n_sim,n_w_mod,n_params) ;
    all_params = fit_data.all_params ; 
    if strcmpi(model_type, 'RW')
        all_w = [all_params(:).a_RW] ; % hyper-parameters 
    else
        all_w = [all_params(:).decay] ;
    end
    
    true_params = fit_data.true_params ; % (n_sim, n_w_mod) 
    learned_param = {all_params(:).learned_param} ; 
    [n_sim, n_w_mod] = size(true_params) ; 
    fitted_params = NaN(n_sim, n_w_mod) ;     
    
    % for param recovery, only consider true model
    idx_curr_mod = cellfun(@(x) strcmpi(models_sim{imod},x), learned_param) ;    
    MSE_tmp = all_MSE(:,:,idx_curr_mod) ;
    curr_w = all_w(idx_curr_mod) ; 
    
    for isim=1:n_sim
        % identify best fitting parameter
        [~,iopt] = min(squeeze(MSE_tmp(isim,:,:)),[],2) ;
        fitted_params(isim,:) = curr_w(iopt) ; 
    end
    
    
    % =================================================================== %
    % ===*=== prepare & plot data
    % =================================================================== %
    if strcmpi(model_type, 'IO')
        fitted_params(fitted_params==Inf) = 410 ; 
        true_params(true_params==Inf) = 410 ; 
    end
    fig_sz = [1 2 12 12] ; taille = 18 ; taille_axis = 14 ; 
    mk_plot = {'o'} ; sz_pts = 9 ;

    fn_fig = [fn_plot, all_fn_fit{imod}] ;
    data_struct = struct() ;
    if log_axes
        full_poly = polyfit(log(true_params(:)),log(fitted_params(:)),1) ; 
    else
        full_poly = polyfit(true_params(:),fitted_params(:),1) ; 
    end
    for isim=1:n_sim
        data_struct(isim).x = true_params(isim,:) ;
        data_struct(isim).y = fitted_params(isim,:) ;
        data_struct(isim).coeffs = full_poly ;
    end   
    
    scatter_colormap_lin_fit(data_struct, 'eps_fig', eps_fig, 'fig_fig', fig_fig, 'pdf_fig', pdf_fig, ...
        'fn_save', fn_fig, 'save_fig',1, 'fig_name', all_fn_fit{imod},...
        'add_lin_fit', 1, 'use_colormap', 0, 'axis_equal', 1, 'add_identity', 1, ...
        'x_lab', 'true parameter', 'y_lab', {'fitted parameter'}, ...
        'fig_sz', fig_sz, 'taille_tick', taille, 'taille_stick', taille_axis, ...
        'sz_pts', sz_pts, 'markers', mk_plot, 'single_fit', 1, ...
        'add_lgd', 0, 'show_data', 1, 'log_axes', log_axes) ;
end


%%
% ======================================================================= %
% ===*=== (2) Model recovery
% ======================================================================= %
n_models_fit = n_models_sim ; 
q_params = 3 ; % 2 for lin reg + 1 for decay
x_descr = get_short_names(models_sim) ;

all_rec_models = zeros(n_models_sim, n_models_fit) ; % rows: true models, columns: fitted models  

for imod = 1:n_models_sim    
    fit_data = load([fn_save, all_fn_fit{imod}, '.mat']) ; 
    all_MSE = fit_data.all_MSE ; % (n_sim,n_w_mod,n_params) ;
    all_params = fit_data.all_params ; 
    n_params = length(all_params) ; 
    n_data = fit_data.n_data ; % (n_sim, n_w_mod) ; 
    [n_sim, n_w_mod] = size(n_data) ; 
    learned_param = {all_params(:).learned_param} ; 
    
    model_BIC = zeros(n_sim, n_w_mod, n_models_fit) ; 
    
    for imod_fit = 1:n_models_fit
        idx_curr_mod = cellfun(@(x) strcmpi(models_sim{imod_fit},x), learned_param) ; 
        MSE_tmp = all_MSE(:,:,idx_curr_mod) ; 
        for isim=1:n_sim
            % select best fitting parameter to define BIC
            MSE_min = (min(squeeze(MSE_tmp(isim,:,:)),[],2))' ;% (n_w_mod, 1)
            model_BIC(isim, :, imod_fit) = n_data(isim,:).*log(MSE_min) + q_params.*...
                log(n_data(isim,:)) ; 
        end
    end
    
    
    for isim=1:n_sim
        BIC_tmp = squeeze(model_BIC(isim,:,:)) ; 
        [~, mod_idx] = min(BIC_tmp,[],2) ; 
        for imod_fit = 1:n_models_fit       
            all_rec_models(imod, imod_fit) = all_rec_models(imod, imod_fit) + ...
                sum(mod_idx==imod_fit) ; 
        end
    end
    all_rec_models(imod,:) = 100.*all_rec_models(imod,:)./(n_sim*n_w_mod) ; 
    
    disp(['[true ',x_descr{imod},'] Percent recovered: ', num2str(all_rec_models(imod,imod))])
end

% =================================================================== %
% ===*=== plot model recovery
% =================================================================== %
str_tmp = all_fn_fit{1} ; str_tmp(4:6) = []; 
fn_fig = [fn_plot, str_tmp] ;

c_bars = viridis(n_models_fit+1) ; c_bars = c_bars(1:(end-1),:) ; % 2:end
% perceptually uniform & colorblind friendly cmaps: 
% viridis, inferno, plasma, magma

plot_bars(all_rec_models, 'x_lab', 'true models', 'y_lab', 'percentage recovered','fn_save',...
    fn_fig, 'fig_fig', fig_fig, 'pdf_fig',pdf_fig, 'eps_fig', eps_fig, ...
    'x_descr', x_descr, 'use_colormap', 0, 'view_3D', 0, ...
    'add_prior', 1, 'prior_val', 100/n_models_fit, 'gr_labels', x_descr, 'c_bars', c_bars) ; 

















