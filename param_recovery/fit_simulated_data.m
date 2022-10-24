function [] = fit_simulated_data()
% Compute the quality of fit on simulated data computed in
% simulate_behavior.m.

% ======================================================================= %
% ===*=== General parameters
% ======================================================================= %
model_type      = 'IO' ;    % RW (Rescorla-Wagner = delta rule), IO (Bayesian)
% simulate data for the different models of IO or RW
% RW: doesn't make sense as we can't sample proba estimates... recovery
% will be perfect

fn_save = ['./data_simu/'] ;
if strcmpi(model_type, 'RW')
    all_fn_simu = {'RW_dAF_10rep_21w', 'RW_dIF_10rep_20w', 'RW_dTP_10rep_21w'} ; 
    
    % parameters to fit
    all_params  = get_all_params_RW() ; 
    n_params    = length(all_params) ; 
else
    all_fn_simu = {'IO_AF_30rep_14w', 'IO_IF_30rep_14w', 'IO_TP_30rep_17w'} ;
    
    % parameters to fit
    all_params  = get_all_params_to_test(0, 0,0) ; 
    n_params    = length(all_params) ; 
end

n_models_sim = length(all_fn_simu) ; 



% ======================================================================= %
% ===*=== Fitting
% ======================================================================= %
in = get_default_IO_in() ;
n_conds_max = 10 ;
n_stim_test = 100 ; 

% loop over (true) simulated models
for imod = 1:n_models_sim
    ts = tic ; 
    simu_data = load([fn_save, all_fn_simu{imod}, '.mat']) ;
    model_p1 = simu_data.model_p1 ; 
    % model_p1{isim, ip}.p1 % NaN(n_stim_test, n_conds_max)
    % model_p1{isim, ip}.s
    % model_p1{isim, ip}.ndata
    
    [n_sim, n_w_mod] = size(model_p1) ;
    all_MSE = NaN(n_sim,n_w_mod,n_params) ;
    true_params = NaN(n_sim, n_w_mod) ;
    n_data = NaN(n_sim, n_w_mod) ; 
    
    for iparam=1:n_params
        disp(['[',all_fn_simu{imod},']    Starting param ',num2str(iparam),'/', num2str(n_params)])
        curr_param = all_params(iparam) ;
        learned_param = curr_param.learned_param ;
        
        if strcmpi(model_type, 'rw')
            a_RW = curr_param.a_RW ;            
        else
            in.learned = learned_param ;
            in.opt.MemParam = curr_param.MemParam ;
            in.opt.AboutFirst = curr_param.AboutFirst ;
        end
        
        for isim=1:n_sim
            for ip = 1:n_w_mod
                all_s = model_p1{isim, ip}.s ; 
                p1_est = model_p1{isim, ip}.p1 ; 
                true_params(isim, ip) = model_p1{isim, ip}.param ; 
                n_data(isim, ip) = model_p1{isim, ip}.ndata ; 
                
                p1_mod_tmp = NaN(n_stim_test, n_conds_max) ; 
                for icond = 1:n_conds_max
                    s = (all_s(:,icond))' ; 
                    trials_resp = find(~isnan(p1_est(:,icond))) ; 
                    
                    % with current sequence, run model
                    if strcmpi(model_type, 'rw')
                        %
                        [p1_mod, ~, ~] = RescorlaWagner(s, learned_param, a_RW) ; 
                    else
                        % get IO outcomes (p1 and p1_std) --> to be able to sample from beta
                        [~, ~, ~, ~, p1_mod, ~, ~, ~] = get_IO_outcomes(in, s) ;
                        % only keep model predictions for some trials
                    end
                    p1_mod_tmp(trials_resp,icond) = p1_mod(trials_resp) ;
        
                end % n_cond
                % for curr ip, isim, iparam, merging all conds:
                % regress estimates on model p1
                y_tmp = make_col_vector(p1_est(~isnan(p1_est))) ; 
                x_tmp = make_col_vector(p1_mod_tmp(~isnan(p1_est))) ;
                [~, mse, ~, ~] = regress_Y_on_X(y_tmp, x_tmp) ; 
                
                all_MSE(isim, ip, iparam) = mse ; 
                
            end % n_w_mod
        end % n_sim
    end
    
    % Save for each model (one learned param per file, w/ all true and fitted params)
    fn_mod_tmp = [fn_save, all_fn_simu{imod}, '_', num2str(n_params),'fits.mat'] ; 
    save(fn_mod_tmp, '-v7.3', 'all_MSE', 'all_params', 'n_data', 'true_params') ;
    
    disp(['=== Time to fit all params for simulated data ', all_fn_simu{imod}, ': ', ...
        num2str(toc(ts)), ' sec'])
end


end

