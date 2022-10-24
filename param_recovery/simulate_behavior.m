function simulate_behavior()
% Simulate behavior using a range of parameters consistent with the ones
% observed in our data set. Parameter recovery will be plotted in
% disp_param_recovery.m. 
% Simulations are performed for the 6 models. 
%
% Running this function may take a few minutes. 

% ======================================================================= %
% ===*=== General parameters
% ======================================================================= %
% number of simulations for each parameter
n_sim           = 30 ;      % for each parameter (Heald2021: 10)
opt_params      = 1 ;       % 1: all indiv optimal param; 
                            % 2. equally spaced params in range defined by subj
                            % issue with 2.: true param may not be in
                            % searched grid for the fitting...
n_params_sim    = 10 ;      % # of parameters to consider for each model (equally spaced)

model_type      = 'IO' ;    % RW (Rescorla-Wagner = delta rule), IO (Bayesian)
% simulate data for the different models of IO or RW
% RW: doesn't make sense as we can't sample proba estimates... recovery is always "perfect"

fn_save = ['./data_simu/'] ;
if ~exist(fn_save, 'dir') ; mkdir(fn_save) ; end

% ======================================================================= %
% ===*=== Reload fit of behavioral data
% ======================================================================= %
% Need data saved in TSL_fit_on_ratings.m to cover realistic range of parameter values. 
if strcmpi(model_type, 'RW')
    % RW models; can replace by data saved in TSL_fit_on_ratings 
    %fn_behav = '../results_ratings_d/all_fits/Fit_31subj_297params_3models.mat' ;
    fn_behav = './Fit_31subj_297params_3models.mat' ;
    data_reloaded   = load(fn_behav) ; 
    all_MSE         = data_reloaded.MSE_RW ; 
    all_MSE         = squeeze(all_MSE(:,1,:)) ; % (n_subj, n_params_RW)
    all_params      = data_reloaded.all_params_RW ;  % n_params = length(all_params) ; 
    all_w           = [all_params(:).a_RW] ; % hyper-parameters 
    %w_val = unique(all_w) ; nW = length(w_val) ;    
    learned_param = {all_params(:).learned_param} ; names_models = unique(learned_param) ; 
    n_models = length(names_models) ; 
    x_descr = strcat('d', get_short_names(names_models)) ;% dAF, dIF, dTP

else
    % IO models
    %fn_behav = '../results_ratings/all_fits/Fit_31subj_309params_3models.mat' ;
    fn_behav = './Fit_31subj_309params_3models.mat' ;
    data_reloaded   = load(fn_behav) ; 
    all_MSE         = data_reloaded.MSE_all_subj ; 
    all_MSE         = squeeze(all_MSE(:,1,:)) ; % (n_subj,n_params)
    all_params      = data_reloaded.all_params ; 
    all_w           = [all_params(:).decay] ;
    % w_val = unique(all_w) ; nW = length(w_val) ;     
    learned_param = {all_params(:).learned_param} ; names_models = unique(learned_param) ; 
    n_models = length(names_models) ; x_descr = get_short_names(names_models) ; % AF, IF, TP
end
in = get_default_IO_in() ;

% =================================================================== %
% ===*=== Identify range of parameters to simulate
% =================================================================== %
true_params = cell(n_models,1) ; 
for imod=1:n_models
    %disp(x_descr{imod})
    iparams = cellfun(@(x) strcmpi(names_models{imod},x), learned_param) ; 
    
    MSE_mod = all_MSE(:,iparams) ;
    w_mod = all_w(iparams) ; 
    [~,iopt] = min(MSE_mod,[],2) ; 
    opt_w = unique(w_mod(iopt)) ; 
    
    % ===*=== true parameters to use
    switch opt_params
        case 1
            % 1. optimal parameters of the participants --  idem Heald2021 p13
            true_params{imod} = opt_w ; 
        case 2
            % 2. equally spaced w among range defined by indiv optimal params
            minp = min(opt_w) ; 
            % identify if inf is present
            any_inf = any(isinf(opt_w)) ; 
            maxp = max(opt_w(~isinf(opt_w))) ; 
            if any_inf
                params_tmp = [linspace(minp, maxp, n_params_sim-1), inf] ; 
            else
                params_tmp = linspace(minp, maxp,n_params_sim) ; 
            end
            true_params{imod} = params_tmp ; 
    end  
    %params = true_params{imod}
end


% ======================================================================= %
% ===*=== Simulate behavioral data
% ======================================================================= %
% p1 & s: (nstim_test, nconds_max) (p1 with nan when no response actually asked)
rng(14) ;
n_stim_test = 100 ; 
n_conds_max = 10 ; 
IRI = 12:18 ; 
diff_TPs = get_TPs_config(0.3,0) ;  
all_TPs = [diff_TPs; diff_TPs] ; % (n_conds_max, 2)

for imod = 1:n_models
    disp(['Starting ',num2str(n_sim),' simulations for model ', model_type, ' ', x_descr{imod}, ...
        '...'])
    ts = tic ;
    
    curr_learned_param = names_models{imod} ; 
    w_mod = true_params{imod} ; 
    n_w_mod = length(w_mod) ; % number of hyper param for current model
    model_p1 = cell(n_sim, n_w_mod) ;
    
    in.learned = curr_learned_param ;  
    
    for ip = 1:n_w_mod
        % loop over hyper parameters
        w_tmp = w_mod(ip) ;
        in.opt.MemParam = {'Decay', [w_tmp]} ; 
        
        for isim = 1:n_sim
            model_p1{isim, ip}.param = w_tmp ; % true parameter (leak or learning rate)
            alls = NaN(n_stim_test, n_conds_max) ;
            allp1 = NaN(n_stim_test, n_conds_max) ;
            ndata = 0 ; % number of p1 for current simu
            for icond=1:n_conds_max
                s = GenRandSeq(n_stim_test, all_TPs(icond,:)) ; % entries in {1, 2}, (n_stim_test, 1)
                % Timing of the questions.
                trials_resp = sample_lists(n_stim_test, 3, 8, IRI) ;
                ndata = ndata + length(trials_resp) ; 
                
                if strcmpi(model_type, 'rw')
                    %
                    [p1_est, ~, ~] = RescorlaWagner(s, curr_learned_param, w_tmp) ; 
                    p1_est = p1_est(trials_resp) ; 
                else
                    % get IO outcomes (p1 and p1_std) --> to be able to sample from beta
                    [~, ~, ~, ~, p1_mean, ~, ~, p1_sd] = get_IO_outcomes(in, s) ;
                    % only keep model predictions for some trials
                    p1_mean = p1_mean(trials_resp) ; p1_sd = p1_sd(trials_resp) ;
                    [p1_a, p1_b] = BetaMomentsToParams(p1_mean, p1_sd);
                    p1_est = betarnd(p1_a, p1_b) ; 
                    % alternative: p1_est = p1_mean
                end
                allp1(trials_resp, icond) = p1_est ; 
                alls(:,icond) = s ; 
            end
            model_p1{isim, ip}.p1 = allp1 ; % simulated p1
            model_p1{isim, ip}.s = alls ; % full seq to apply different model afterwards
            model_p1{isim, ip}.ndata = ndata ; 
        end
    end
    
    
    % Save for each model (one learned param per file, w/ all tested params)
    fn_mod_tmp = [fn_save, model_type, '_', x_descr{imod},'_', num2str(n_sim),...
        'rep_',num2str(n_w_mod), 'w.mat'] ; 
    save(fn_mod_tmp, '-v7.3', 'model_p1') ;
    
    disp(['=== Time to simulate data for model ', model_type, ' ', x_descr{imod}, ...
        ': ', num2str(toc(ts)), ' s'])
end


% fitting: amounts to comparing the MSE for different models (cfr TSL_fit_on_ratings)
end

