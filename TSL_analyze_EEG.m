function [all_TPs_str,indices_TP] = TSL_analyze_EEG(varargin)
% Analyze EEG data of the experiments and fit learning models. 
% Data are saved and can be relaoded and plotted in other functions (such
% as TSL_plot_avg_EEG.m, TSL_plot_IO_fit.m, etc). 

%close all ; 

% ======================================================================= %
% ===== * ===== General parameters
% ======================================================================= %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
all_subj_idx    = setdiff(1:36,[1, 11, 15, 28, 33]) ; 
bcorr           = 1 ;       % baseline correct the signals
LP_filter       = 1 ;       % Low-pass filter the epochs
LP_freq         = 30 ;      % frequency cutoff for the LP filtering
butter_order    = 4 ;       % order of the Butterworth filter
reject_ep       = 1 ;       % 
reject_ampl     = 1 ;       % reject epochs based on amplitude
max_ampl_rej    = 80 ;      % amplitude criterion to reject epochs
                            % applied after LP AND after bc... % 150 when
                            % no LP, 80 otherwise
reject_usual    = 0 ;       % reject usual epochs based on signal with LP30 >80
save_avg_eeg    = 0 ;       % save the avg EEG per condition, per param, 
                            % and merging all conditions
save_AUC_VW     = 0 ;       % AUC of N2 and P2 + amplitudes of N2/P2 peaks
save_avg_TF     = 1 ;       % save TF avg across epochs and pre-alpha power (return after that)
if save_avg_TF; LP_filter = 0 ; max_ampl_rej = 150 ; end
save_EEG_model  = 0 ;       % save all the EEG and regressors info necessary to fit linear models
moving_avg      = 0 ;       % moving avg on each EEG *trial*
                            % ONLY for BMC (IO_fit_opt = 2)
movmu_ms        = 80 ;      % msec for the moving avg
reref           = 0 ;       % re-reference the EEG
elecs_ref = {'T7','T8'} ; %'PZ' ; 
reject_first    = 0 ;       % reject the first epochs
reject_perI     = 1 ;       % n_first_rej first epochs are rejected for each I
n_first_rej     = 2 ;       % # of initial epochs to reject
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IO_fit_opt      = 1 ;       % 0: no IO fitting
                            % 1: fitting of one parameter (selected below)
                            % 2: fitting all considered params to do BMC. 
without_AF      = 1 ;       % Don't consider AF in the params (when IO_fit_opt = 2)
single_elec     = 0 ;       % IO fitting at a single electrode
elec_IO         = 'FCZ' ;   % FCZ
random_inter    = 1 ;       % random intercepts for the different conditions (TPs).
conf_prev       = 1 ;       % use confidence in t-1 as predictor in t
save_ep_kept    = 0 ;       % Save indices of ep_kept for each subject
CB_PT           = 0 ;       % cluster-based permutation testing to correct for MC
n_perm          = 2000 ;    % nb of permutations for CB-PT
cf_alpha        = 0.05 ;    % cluster-forming p-value  
rconf_split     = 1 ;       % for the avg_split
rconf_BMC       = 1 ;       % rconf or conf for BMC (IO_fit_opt=2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TF parameters to analyze oscillations
TF_regress      = 0 ;       % do the IO regressions on the TF data
                            % only for IO_fit_opt = 1. 
use_stft        = 1 ;       % otherwise: cwt
win_width_sec   = 0.2 ;     % seconds (if use_stft)
                            % 0.3: 6 oscillations at 20Hz
stft_with_conv = 1 ;        % if use_stft % 5 oscillations @70Hz = 0.07sec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_begin = 6 ; f_step = 2 ; f_stop = 140 ; 

% TO KEEP TO 0 %%%%%%%%%%%%%%
plot_avg_eeg    = 0 ;       % plot the avg EEG responses (avg for I1 and I2)
                            % can do it in TSL_plot_avg_EEG.m
p_predictor     = 0 ;       % Add the IO estimate of TP as predictor
rconf_no_logN   = 0 ;       % Add log(N) as a regressor to compute residual conf
% would lead to save too heavy files
save_data       = 0 ;       % save the merged EEG data of all subjects
reload_data     = 0 ;       % if file previously saved, reload the merged
                            % EEG data of all subjects.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if use_stft
    use_cwt = 0 ; 
else
    use_cwt = 1 ; 
end

% ======================================================================= %
% ===== * ===== Ideal Observer: Default parameters
% ======================================================================= %
% === when IO_fit_opt = 1
use_indiv_mod   = 0 ;       % models fitted individually on ratings
use_indiv_leak  = 0 ;       % leaks fitted individually on ratings, only if use_indiv_mod

MemParam = {'Decay', [8]} ; 
AboutFirst      = 'WithoutFirst' ;
learned_param   = 'transition' ; %'transition' ; % 'frequency', 'transition', 'alternation'

[window_int, windowF_int, window, windowF, leaky_int, leakyF_int, ...
    decay, decayF, weight, ~] = get_memParams(MemParam) ; 

str_model = get_model_param_str(AboutFirst, learned_param, decay, decayF, ...
    window, windowF, leaky_int, leakyF_int, window_int, windowF_int, weight) ;
n_subj = length(all_subj_idx) ; 
if use_indiv_mod
    data_indiv = load(['./results_ratings/data/opt_models_',num2str(n_subj),'subj.mat']) ; 
    names_indiv_mod = data_indiv.names_indiv_mod ; indiv_decays = data_indiv.indiv_decays ;     
    str_model = ['_opt', str_model(4:end)] ; 
    if use_indiv_leak
        idx_leak = strfind(str_model, 'leak') ;
        str_model = [str_model(1:(idx_leak-1)), 'opt'] ;        
    end        
end

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
in.opt.priorp1g2 = [1 1] +eps_prior;% uniform Beta prior (when learning TP)
                                    % Add eps to 1 to the prior: avoid NaN
                                    % in MAP when learning TP and AF after
                                    % the first observation (cfr MAP in
                                    % mar_ComputeMAPandPrediction and
                                    % bern_ComputeMAPandPrediction). 
in.opt.priorp2g1 = [1, 1]+eps_prior;% uniform Beta prior (when learning TP)
in.opt.priorp1 = [1,1]+eps_prior ;  % uniform Beta prior (when learning IF, AF)
in.verbose = 1;                     % to check that no default values are used.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======================================================================= %

str_stim        = '_TCS' ;
n_stim_test     = 100 ; 
n_conds_max     = 10 ;  
fn_dir_EEG      = './data_EEG/' ; 
fn_dir_ratings  = './data_ratings/' ; 


fn_res = ['./results_EEG/'] ; str_lica = '' ; 
if reject_first
    if reject_perI
        fn_res = [fn_res(1:end-1), '_wo', num2str(n_first_rej), 'I/'] ; 
    else
        fn_res = [fn_res(1:end-1), '_wo', num2str(n_first_rej), '/'] ; 
    end
end

if reref 
    str_ref = '' ; 
    for ie=1:length(elecs_ref)
        str_ref = [str_ref,'_', elecs_ref{ie}] ;
    end
    fn_res = [fn_res, str_ref(2:end), '/'] ; 
end

all_TPs = [0.5, 0.5; 0.3, 0.7; 0.7, 0.3; 0.3, 0.3; 0.7, 0.7] ; 
n_TPs = size(all_TPs,1) ; all_TPs_str = cell(1, n_TPs) ; 
for i_tp = 1:n_TPs
    all_TPs_str{i_tp} = get_str_TP(all_TPs(i_tp,:)) ;
end


% General plot params --------------------------------------------------- %
eps_fig = 0 ; fig_fig = 0 ; pdf_fig = 0 ; % otherwise: png
taille = 18 ; taille_axis = 14 ; 
% Specific plot params -------------------------------------------------- %
plot_ms = 1 ;               % msec time unit

% Define channels to show on plot
only_plot_topo_chan = 1 ;   % don't show channels out of topoplot in butterfly plot
topo_colors = 1 ;           % show color codes for the avg signals on topo
selected_chan = 1 ;         % only plot chan_selected on the butterfly plot
chan_selected = {'FCZ', 'C3', 'C4', 'CPz', 'Cz'} ;
% ----------------------------------------------------------------------- %


% === Parse varargin
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
            case 'all_subj_idx'
                all_subj_idx = Val_arg ;
            case 'only_plot_topo_chan'
                only_plot_topo_chan = Val_arg ;
            case 'selected_chan'
                selected_chan = Val_arg ;
            case 'chan_selected'
                chan_selected = Val_arg ;
        end
    end
end
% ======================================================================= %
% ======================================================================= %

str_TF = '' ; 
if TF_regress
   moving_avg = 0 ;  LP_filter = 0 ; 
   str_TF = '_TF' ; 
   %single_elec = 1 ; 
   max_ampl_rej = 150 ; %150 ; % without LP30, must increase max ampl!
   if use_stft
       str_TF = [str_TF, '_w', num2str(round(100*win_width_sec))] ; % window width
   end
end
if IO_fit_opt==2
   moving_avg = 1 ;  
end

n_chans_plot = length(chan_selected) ;

if ~exist(fn_res, 'dir') % 0 or 7 if it exists
    mkdir(fn_res) ; 
end
 
% ======================================================================= %
% ===== * ===== Initialize data structures
% ======================================================================= %
all_s = NaN(n_stim_test, n_conds_max, n_subj) ; % sequences, in {1, 2}
all_trials_resp = false(n_stim_test, n_conds_max, n_subj) ; 
indices_TP = NaN(n_conds_max,n_subj); % indices for all_TPs and all_TPs_str 
% giving condition for each bloc
% all_EEG_data & all_GSR: def based on nb of time samples

switch IO_fit_opt
    case 1
        n_IO_out = 7 ; 
        % prediction error - confidence - update - unpredictability - p1 -
        % surprise - N_used (= N observations for current
        % inference)
        all_IO_data = NaN(n_IO_out, n_stim_test, n_conds_max, n_subj) ;
        n_params = 1 ; 
    case 2
        n_IO_out = 4 ; 
        % PE - confidence - p1 (p1: to compute residual confidence) -
        % N_used
        
        all_params = get_all_params_EEG(0,0,0, without_AF) ; 
        n_params = length(all_params) ;

        all_IO_data = NaN(n_IO_out, n_stim_test, n_conds_max, n_subj, n_params) ;
end




if n_subj>1
    str_subj = ['_',num2str(n_subj),'subj'] ; 
else
    str_subj = ['_subj',num2str(num2str_padd(all_subj_idx(1), 2))] ;
end
str_LP = '' ; 
if LP_filter
    str_LP = ['_LP', num2str(LP_freq)] ; 
end
str_plot = '' ; 
if selected_chan
    str_plot = ['_',num2str(n_chans_plot),'chans'] ; 
elseif only_plot_topo_chan
    str_plot = '_topo_chans' ; 
end

fn_data_tmp = [fn_res, 'data/'] ; 
if ~exist(fn_data_tmp, 'dir') % 0 or 7 if it exists
    mkdir(fn_data_tmp) ; 
end

str_movmu = '' ; 
if moving_avg
    str_movmu = ['_mmu', num2str(movmu_ms)] ;
end
fn_ep_usual = [fn_data_tmp, 'Ep_kept_31subj_LP30',str_movmu,'_maxA80.mat'] ; 

str_ampl_max = '' ; 
if reject_ampl
    if reject_usual && exist(fn_ep_usual, 'file')
        str_ampl_max = ['_maxAU'] ; 
    else
        str_ampl_max = ['_maxA', num2str(max_ampl_rej)] ; 
    end
end

str_params = [str_subj, str_LP, str_movmu] ;

fn_saved_data = [fn_data_tmp, 'Data', str_params, str_ampl_max, '.mat'] ; 
fn_ep_kept = [fn_data_tmp, 'Ep_kept', str_params, str_ampl_max, '.mat'] ; 


if reject_usual && exist(fn_ep_usual, 'file')
    disp(['--- Rejecting usual epochs (based on LP30 data with maxA 80'])
    tmp = load(fn_ep_usual) ; 
    all_ep_kept = tmp.all_ep_kept ; % used in TSL_plot_BMC
    compute_ep_rej = 0 ; save_ep_kept = 0 ;
else
    compute_ep_rej = 1  ; save_ep_kept = 1 ; 
    all_ep_kept = true(n_stim_test, n_conds_max, n_subj) ; % logical to keep or not each epoch!
end

fn_params_data = [fn_data_tmp, 'Params', str_params, str_ampl_max, '.mat'] ; % without the EEG data
fn_avg_eeg = [fn_data_tmp, 'Avg_eeg', str_params, str_ampl_max, '.mat'] ; 

str_inter = '' ; 
if random_inter
   str_inter = '_rinter' ;  
end

str_p = str_TF ; 
if p_predictor && IO_fit_opt~=2 ; str_p = [str_p,'_pp'] ; end
if rconf_no_logN ; str_p = [str_p, '_noN'] ; end % noN = no Nb of observations % && IO_fit_opt~=2
if conf_prev ; str_p = [str_p, '_pc'] ; end % pc = previous confidence used
if single_elec ; str_p = [str_p, '_', elec_IO] ; end

str_conf_split = '' ; 
if rconf_split ; str_conf_split = ['_rc',str_p] ; end
fn_avg_split = [fn_data_tmp, 'Avg_split',str_params, str_ampl_max, str_model, ...
    str_inter, str_conf_split, '.mat'] ; 

if CB_PT && IO_fit_opt~=2 
    str_p = [str_p, '_',num2str(n_perm),'CP'] ; 
    if cf_alpha~=0.05
        str_p = [str_p, '_', num2str(round(1000*cf_alpha))] ; 
    end
end

if IO_fit_opt==2
    str_rc_bmc = '' ; 
    if rconf_BMC ; str_rc_bmc = '_rc' ; end
    fn_IO_fit = [fn_data_tmp, 'IO_fit', str_params, str_ampl_max, str_rc_bmc, ...
        '_', num2str(n_params), 'params', str_inter, str_p, '.mat'] 
    % [~, chan_desired] = define_chan_order_colors(chanlocs, chan_selected) ;
else
    fn_IO_fit = [fn_data_tmp, 'IO_fit', str_params, str_ampl_max, str_model, ...
        str_inter, str_p, '.mat'] 
end


if reload_data && exist(fn_saved_data, 'file')
% ======================================================================= %
% ===== * ===== Load pre-processed data
% ======================================================================= %
    disp(['--- Reloading EEG data of all subjects', ' ---'])
    t_s = tic ; 
    saved_eeg_data = load(fn_saved_data) ;
    disp(['=== Time to reload data: ', num2str(toc(t_s))])
    all_EEG_data = saved_eeg_data.all_EEG_data;
    all_s = saved_eeg_data.all_s;
    all_ep_kept = saved_eeg_data.all_ep_kept;
    indices_TP = saved_eeg_data.indices_TP;
    time_vec = saved_eeg_data.time_vec;
    fs = saved_eeg_data.fs;
    xstep = saved_eeg_data.xstep;
    chanlocs = saved_eeg_data.chanlocs ; 
    
    [n_chan, n_times, n_stim_test, n_conds_max, n_subj] = size(all_EEG_data) ; 
else
% ======================================================================= %
% ===== * ===== Load raw data
% ======================================================================= %
    t_IO = 0 ;
    for i_subj = 1:n_subj
        subj_number = all_subj_idx(i_subj) ; 
        subj_idx = num2str(subj_number) ; 
        subj_idx_padd = num2str_padd(subj_number, 2) ; 
        param_str = [str_stim, '_subj',subj_idx] ; 
        str_subj_dir = ['subj',subj_idx_padd,'/'] ;
        str_subj_dir_EEG = '' ; 

        % --------------------------------------------------------------- %
        % --*-- EEG data
        % --------------------------------------------------------------- %
        fn_EEG_data = [fn_dir_EEG, str_subj_dir_EEG, 'Subject', subj_idx_padd, str_lica,'.mat'] ; 
        disp(['===== Loading data ', fn_EEG_data])
        lwdata_subj = load(fn_EEG_data) ;
        lwdata_subj = lwdata_subj.lwdata ;
        disp(['===== Data loaded for subject ', subj_idx, ' (',num2str(i_subj),...
            '/',num2str(n_subj),')'])

        % Extract the time vec (once for all data)
        if i_subj==1
            header = lwdata_subj(1,1).header ;
            xstep = header.xstep ; 
            xstart = header.xstart ;
            fs = 1/xstep ;	
            idx_tmp = 1:header.datasize(1,6) ; 
            time_vec = xstart + xstep.*(idx_tmp-1) ;
            
            nbins_mmu = 2*floor(fs*movmu_ms/1000/2)+1 ;

            % ==*== Initialize matrix to store EEG data
            chanlocs = header.chanlocs ; 
            n_chan = header.datasize(1,2) ; 
            n_times = length(time_vec) ;

            all_names = {chanlocs.labels} ;
            idx_GSR = find(contains(all_names, 'GSR')) ;  
            idx_elecs = setdiff(1:n_chan, idx_GSR) ; 
            chanlocs = chanlocs(idx_elecs) ; 
            if ~isempty(idx_GSR)
                n_chan = n_chan - 1 ;     
            end

            all_EEG_data = NaN(n_chan, n_times, n_stim_test, n_conds_max, n_subj) ; 
            %all_GSR = NaN(n_times, n_stim_test, n_conds_max, n_subj) ; 
        end

        % --------------------------------------------------------------- %
        % --*-- Test conditions and order
        % --------------------------------------------------------------- %
        fn_cond_subj = [fn_dir_ratings, str_subj_dir] ;     
        fn_cond_data = [fn_cond_subj,'Cond',param_str,'.mat'] ;
        cond_data = load(fn_cond_data) ; 
        all_TP_test = cond_data.all_TP_test ;  
        TP_order_full = cond_data.TP_order_full ; 
        order_used = TP_order_full ; % TP_order ;   
        n_conds = length(order_used) ;  

        % --------------------------------------------------------------- %
        % --*-- Load data of all conditions
        % --------------------------------------------------------------- %
        for i_cond = 1:n_conds
            cond_idx_padd = num2str_padd(i_cond, 2) ; 
            TP_tmp = all_TP_test(order_used(i_cond),:) ; 
            str_TP = get_str_TP(TP_tmp) ;
            fn_test = [fn_cond_subj, 'Test_cond',num2str(i_cond), param_str, ...
                str_TP, '_',num2str(n_stim_test),'stim_1.mat'] ; 

            % ===== * ===== TP already found ? 
            idx_TP = find(strcmpi(all_TPs_str, str_TP),1) ; 
            indices_TP(i_cond,i_subj) = idx_TP ; 

            % ===== * ===== Load and extract useful information

            % ------> from the EEG
            % EEG names: 'Testii_TPjj_kk'
            % --> COULD want to analyze EEG as a fct of session index OR TP...
            % cfr find_idx_modality
            % str_to_find: str_TP OR Testii
            idx_EEG_cond = find_idx_conds(lwdata_subj, ['Test', cond_idx_padd]) ; 

            if length(idx_EEG_cond)>1
               error(['==[Subj ', subj_idx, ']==  ', num2str(length(idx_EEG_cond)),...
                   ' test conditions of index ', cond_idx_padd,']']) ; 
            end

            EEG_data = squeeze(lwdata_subj(1,idx_EEG_cond).data) ;
            %curr_all_epoch = squeeze(curr_data(:,chan_num, 1,1,1,:)) ;
            % lwdata: [epoch, channels, index, z,y,x] 
            %     --> (epoch, channels, x)
            EEG_tmp = permute(EEG_data(:,idx_elecs,:), [2,3,1]) ; % Size: (chan, x, epoch)

            % ===== * ===== LP filtering
            if LP_filter  
                if i_cond==1 ; t_end_filter = 0 ; end            
                t_filter = tic ; 

                [z,p,k] = butter(butter_order,LP_freq./(fs/2),'low');
                [sos,G] = zp2sos(z,p,k);
                %[b,a]=butter(butter_order,LP_freq./(fs/2),'low');

                n_ep = size(EEG_tmp,3) ; 
                for i_chan=1:n_chan
                    EEG_chan_tmp = squeeze(EEG_tmp(i_chan,:,:)) ;
                    for i_ep = 1:n_ep
                        EEG_tmp(i_chan,:,i_ep) = filtfilt(sos, G, double(EEG_chan_tmp(:,i_ep)));
                        %EEG_tmp(i_chan,:,i_ep) = filtfilt(b, a,EEG_chan_tmp(:,i_ep));                        
                    end
                end

                t_end_filter = t_end_filter + toc(t_filter) ;             
                if i_cond==n_conds
                    disp(['    [Signals filtered (LP',num2str(LP_freq),')]: ',...
                        num2str(t_end_filter),' seconds for ',num2str(n_conds),' conds ']) ;
                end
            end
            % ===== * ===== Baseline correction
            if bcorr
                idx_base = time_vec<=0 ; 
                EEG_base = mean(EEG_tmp(:,idx_base,:),2) ;
                EEG_tmp = EEG_tmp - repmat(EEG_base,[1,n_times,1]) ; 
            end

            % ------> Sequence & index of behavioral questions
            [s, trials_resp] = load_extract_sequence_info(fn_test) ;
            
            % ===== * ===== Identify epochs to reject
            if compute_ep_rej
                curr_ep_rej = get_rejected_epochs(subj_number, i_cond) ; 
                if reject_ampl
                    if i_cond==1 ; str_n_rej = '' ; end 
                    idx_ampl = squeeze(sum(sum(abs(EEG_tmp)>=max_ampl_rej,2),1))>0 ; 
                    n_rej = sum(idx_ampl) ; str_n_rej = [str_n_rej, ' ', num2str(n_rej)] ; 
                    if i_cond==n_conds
                        disp(['    [#ep w ampl > ',num2str(max_ampl_rej),']: ', ...
                            str_n_rej]) ; 
                    end
                    curr_ep_rej = [curr_ep_rej, find(idx_ampl)'] ; 
                end  
                if reject_first
                    if reject_perI
                        idx_I1_tmp = make_row_vector(find(s==1, n_first_rej)) ; 
                        idx_I2_tmp = make_row_vector(find(s==2, n_first_rej)) ; 
                        curr_ep_rej = [curr_ep_rej, idx_I1_tmp, idx_I2_tmp] ; 
                    else
                        curr_ep_rej = [curr_ep_rej, 1:n_first_rej] ; 
                    end
                end
                all_ep_kept(curr_ep_rej,i_cond, i_subj) = false ; 
            end


            % all_EEG_data: (n_chan, n_times, n_stim_test, n_conds_max, n_subj) 
            if moving_avg
                all_EEG_data(:,:,:,i_cond,i_subj) = movmean(EEG_tmp,nbins_mmu,2) ; 
            else
                all_EEG_data(:,:,:,i_cond,i_subj) = EEG_tmp ; 
            end
            if reref
               [~, idx_elec] = define_chan_order_colors(chanlocs, elecs_ref) ; 
               all_EEG_data(:,:,:,i_cond,i_subj) = EEG_tmp - repmat(mean(EEG_tmp(idx_elec,:,:),1), n_chan, 1, 1) ;
            end
            
            %all_GSR(:, :, i_cond, i_subj) = (squeeze(EEG_data(:,idx_GSR,:)))' ; 

            

            all_codes = {lwdata_subj(1,idx_EEG_cond).header.events.code} ; 
            s_from_EEG = cellfun(@(x) isempty(x), strfind(all_codes, 'I1')) ; 
            s_from_EEG = s_from_EEG + 1 ;
            if sum(abs(s-s_from_EEG))>0
               error('--- sequence s is not correctly encoded in the triggers ---') 
            end
            all_s(:, i_cond, i_subj) = s ;  %NaN(n_stim_test, n_conds_max, n_subj) ;
            all_trials_resp(trials_resp, i_cond, i_subj) = true ; 


            % ===== * ===== IO: extract outcomes of interest
            switch IO_fit_opt
                case 1
                    ts = tic ; 
                    if use_indiv_mod
                        in.learned = names_indiv_mod{i_subj} ; 
                        if use_indiv_leak
                            MemParam = {'Decay', [indiv_decays(i_subj)]} ; 
                            in.opt.MemParam = MemParam ;
                        end
                    end
                    
                    [surprise, distUpdate, conf_p1, unpred, p1_mean, ...
                        pred_error, N_used] = get_IO_outcomes(in, s) ; 
                    if conf_prev
                        % Shift all confidence & p1 to the right & N
                        conf_p1 = conf_p1([1,1:end-1]) ;
                        p1_mean = p1_mean([1,1:end-1]) ;
                        N_used = N_used([1,1:end-1]) ; %
                        %stupid_N_for_conf: was only shifted once, for 1st
                        %cond!!!
                    end
                    % /!\ switch PE and surprise
                    all_IO_data(1, :, i_cond, i_subj) = pred_error ; % prediction error 
                    all_IO_data(2, :, i_cond, i_subj) = conf_p1 ;
                    all_IO_data(3, :, i_cond, i_subj) = distUpdate ;
                    all_IO_data(4, :, i_cond, i_subj) = unpred ;
                    all_IO_data(5, :, i_cond, i_subj) = p1_mean ;
                    all_IO_data(6, :, i_cond, i_subj) = surprise ;
                    all_IO_data(7, :, i_cond, i_subj) = N_used ; 
                    % PE: cfr Bayesian predictive coding in Aitchison2017            
                    % surprise - confidence - update - unpredictability - p1 -
                    % prediction error
                    t_IO = t_IO + toc(ts) ; 
                    
                case 2
                    ts = tic ; 
                    for i_param=1:n_params
                        curr_param = all_params(i_param) ; 

                        in.learned = curr_param.learned_param ;  
                        in.opt.MemParam = curr_param.MemParam ; 
                        in.opt.AboutFirst = curr_param.AboutFirst ; 
                                            
                        [surprise, ~, conf_p1, ~, p1_mean, pred_error, N_used] = ...
                            get_IO_outcomes(in, s) ; 
                        if conf_prev
                            % Shift all confidence & p1 to the right
                            conf_p1 = conf_p1([1,1:end-1]) ;
                            p1_mean = p1_mean([1,1:end-1]) ;
                            N_used = N_used([1,1:end-1]) ; 
                        end
                        
                        all_IO_data(1, :, i_cond, i_subj, i_param) = pred_error ; %surprise ;
                        all_IO_data(2, :, i_cond, i_subj, i_param) = conf_p1 ;
                        all_IO_data(3, :, i_cond, i_subj, i_param) = p1_mean ;
                        all_IO_data(4, :, i_cond, i_subj, i_param) = N_used ;
                        % surprise - confidence - p1 - N_used
                    end
                    t_IO = t_IO + toc(ts) ; 
            end

        end    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    if IO_fit_opt~=0
        disp(['=== Time for IO calls (',num2str(n_params),' params): ', num2str(t_IO), ' sec'])
    end
    
    %save(fn_params_data, '-v7.3', 'indices_TP', 'time_vec', 'fs', 'xstep', ...
    %    'chanlocs', 'all_TPs', 'all_TPs_str') ; 
    if save_data
        % ===== * ===== Save useful information  
        t_s = tic ;
        all_EEG_data = single(all_EEG_data) ; 
        save(fn_saved_data, '-v7.3', 'all_EEG_data', 'all_s', ...
            'all_ep_kept', 'indices_TP', 'time_vec', 'fs', 'xstep', 'chanlocs') ;
        t_save = toc(t_s) ; 
        disp(['=== Time to save data: ', num2str(t_save)]) ; 
    end
    if save_ep_kept
        save(fn_ep_kept, '-v7.3','all_ep_kept', 'all_s') ;
    end
% ======================================================================= %
end

N_obs = n_stim_test*n_conds_max ; 


% ======================================================================= %
% ===== * ===== Extract average EEG of each I
% ======================================================================= %

% all_EEG_data: (n_chan, n_times, n_stim_test, n_conds_max, n_subj)
% all_s:                         (n_stim_test, n_conds_max, n_subj)
y_I1 = NaN(n_chan, n_times, n_subj) ;
y_I2 = NaN(n_chan, n_times, n_subj) ;


if save_avg_eeg
    % ----- idem for each TP separately (to save)
    y_I1_TP = NaN(n_chan, n_times, n_subj, n_TPs) ; 
    y_I2_TP = NaN(n_chan, n_times, n_subj, n_TPs) ; 

    % ----- idem for each cond separately
    y_I1_cond = NaN(n_chan, n_times, n_subj, n_conds_max) ; 
    y_I2_cond = NaN(n_chan, n_times, n_subj, n_conds_max) ; 
        
end
if save_AUC_VW
    AUC_VW = NaN(n_chan, n_stim_test, n_conds_max,n_subj, 3) ; % AUC of N2, AUC of P2, AUC of N2-P2 
    ampl_VW = NaN(n_chan, n_stim_test, n_conds_max,n_subj, 3) ; % amplitudes of N2, P2 and N2-P2 waves
    mampl_VW = NaN(n_chan, n_stim_test, n_conds_max,n_subj, 3) ; % idem for mean amplitudes
    % base on avg EEG
    inter_AUC_I1 = [170, 236; 240, 410; 170, 410] ; % time intervals for N2, P2 and N2-P2 complex -- cfr TSL_plot_avg_split
    inter_AUC_I2 = [316, 450; 432, 668; 316, 668] ; 
    % based on TSL_plot_IO_fit.m --> singificant intervals
    inter_AUC_I1 = [208, 256; 284, 396; 208, 396] ; % time intervals for N2, P2 and N2-P2 complex -- cfr TSL_plot_avg_split
    inter_AUC_I2 = [300, 404; 456, 588; 300, 588] ; 
    ind_AUC_I1 = msec_to_ind(time_vec.*1000, inter_AUC_I1) ; 
    ind_AUC_I2 = msec_to_ind(time_vec.*1000, inter_AUC_I2) ; 
end
for i_subj=1:n_subj
    for i_chan=1:n_chan
        y_tmp = squeeze(all_EEG_data(i_chan, :, :, :, i_subj)) ; 
        % (n_times, n_stim_test, n_conds_max)
        s_tmp = squeeze(all_s(:,:,i_subj)) ; 
        % (n_stim_test, n_conds_max)
        
        if save_avg_eeg
            % ==*== Selecting each condition
            for i_cond = 1:n_conds_max
                y_tmp_cond = reshape(y_tmp(:,:,i_cond),[n_times,n_stim_test]) ; 
                s_tmp_cond = reshape(s_tmp(:,i_cond),[1,n_stim_test]) ; 

                if reject_ep
                    ep_kept_tmp = reshape(squeeze(all_ep_kept(:,i_cond,i_subj)), ...
                        [1, n_stim_test]) ; 
                    y_tmp_cond = y_tmp_cond(:,ep_kept_tmp) ; 
                    s_tmp_cond = s_tmp_cond(ep_kept_tmp) ; 
                end
                % Mean over all the I1/I2 epochs
                y_I1_cond(i_chan, :,i_subj, i_cond) = mean(y_tmp_cond(:,s_tmp_cond==1),2) ; 
                y_I2_cond(i_chan, :,i_subj, i_cond) = mean(y_tmp_cond(:,s_tmp_cond==2),2) ;
            end


            % ==*== Selecting each TP
            for i_tp = 1:n_TPs
                i_conds = find(indices_TP(:,i_subj)==i_tp) ;
                y_tmp_tp = reshape(y_tmp(:,:,i_conds),[n_times,n_stim_test*2]) ; 
                s_tmp_tp = reshape(s_tmp(:,i_conds),[1,n_stim_test*2]) ; 

                if reject_ep
                    ep_kept_tmp = reshape(squeeze(all_ep_kept(:,i_conds,i_subj)), ...
                        [1, n_stim_test*2]) ; 
                    y_tmp_tp = y_tmp_tp(:,ep_kept_tmp) ; 
                    s_tmp_tp = s_tmp_tp(ep_kept_tmp) ; 
                end
                % Mean over all the I1/I2 epochs
                y_I1_TP(i_chan, :,i_subj, i_tp) = mean(y_tmp_tp(:,s_tmp_tp==1),2) ; 
                y_I2_TP(i_chan, :,i_subj, i_tp) = mean(y_tmp_tp(:,s_tmp_tp==2),2) ;
            end
        end
        
        % ==*== Merging all conditions
        y_tmp = reshape(y_tmp,[n_times,N_obs]) ; 
        s_tmp = reshape(s_tmp,[1,N_obs]) ; 
        
        if reject_ep
            ep_kept_tmp = reshape(squeeze(all_ep_kept(:,:,i_subj)), ...
                [1, N_obs]) ; 
            % all_ep_kept = NaN(n_stim_test, n_conds_max, n_subj) ; 
            y_tmp = y_tmp(:,ep_kept_tmp) ;             
            
            lin_idx_I1 = find(ep_kept_tmp & s_tmp==1) ;
            lin_idx_I2 = find(ep_kept_tmp & s_tmp==2) ;
            s_tmp = s_tmp(ep_kept_tmp) ; 
        end
        
        % Mean over all the I1/I2 epochs
        y_I1(i_chan, :,i_subj) = mean(y_tmp(:,s_tmp==1),2) ; 
        y_I2(i_chan, :,i_subj) = mean(y_tmp(:,s_tmp==2),2) ;
        
        if save_AUC_VW
            % AUC_VW: (n_chan, n_stim_test, n_conds_max,n_subj, 3)
            
            % ======================================== I1
            % I1 & I2: different windows to consider
            y_tmp_I1 = abs(y_tmp(:,s_tmp==1)) ; 
            nep_tmp = size(y_tmp_I1,2) ; 
            [ind_stim, ind_cond] = ind2sub([n_stim_test, n_conds_max], lin_idx_I1) ; 
            
            % AUC of N2 -- fills AUC_VW(i_chan, :, :, i_subj, 1)
            AUC_tmp = trapz(y_tmp_I1(ind_AUC_I1(1,1):ind_AUC_I1(1,2),:)) ; 
            AUC_VW(sub2ind([n_chan, n_stim_test, n_conds_max,n_subj, 3], ...
                repmat(i_chan, 1, nep_tmp), ind_stim, ind_cond, repmat(i_subj, 1, nep_tmp), repmat(1, 1, nep_tmp) )) = AUC_tmp ;                
            % AUC of P2 -- AUC_VW(i_chan, :, :, i_subj, 2)
            AUC_tmp = trapz(y_tmp_I1(ind_AUC_I1(2,1):ind_AUC_I1(2,2),:)) ; 
            AUC_VW(sub2ind([n_chan, n_stim_test, n_conds_max,n_subj, 3], ...
                repmat(i_chan, 1, nep_tmp), ind_stim, ind_cond, repmat(i_subj, 1, nep_tmp), repmat(2, 1, nep_tmp) )) = AUC_tmp ; 
            % AUC of N2-P2 -- AUC_VW(i_chan, :, :, i_subj, 3)
            AUC_tmp = trapz(y_tmp_I1(ind_AUC_I1(3,1):ind_AUC_I1(3,2),:)) ; 
            AUC_VW(sub2ind([n_chan, n_stim_test, n_conds_max,n_subj, 3], ...
                repmat(i_chan, 1, nep_tmp), ind_stim, ind_cond, repmat(i_subj, 1, nep_tmp), repmat(3, 1, nep_tmp) )) = AUC_tmp ; 
            
            %%% idem for mean amplitudes 
            mampl_N2 = mean(y_tmp_I1(ind_AUC_I1(1,1):ind_AUC_I1(1,2),:),1) ; % N2
            mampl_VW(sub2ind([n_chan, n_stim_test, n_conds_max,n_subj, 3], ...
                repmat(i_chan, 1, nep_tmp), ind_stim, ind_cond, repmat(i_subj, 1, nep_tmp), repmat(1, 1, nep_tmp) )) = mampl_N2 ;   
            mampl_P2 = mean(y_tmp_I1(ind_AUC_I1(2,1):ind_AUC_I1(2,2),:),1) ; % P2
            mampl_VW(sub2ind([n_chan, n_stim_test, n_conds_max,n_subj, 3], ...
                repmat(i_chan, 1, nep_tmp), ind_stim, ind_cond, repmat(i_subj, 1, nep_tmp), repmat(2, 1, nep_tmp) )) = mampl_P2 ;  
            mampl_VW(sub2ind([n_chan, n_stim_test, n_conds_max,n_subj, 3], ...
                repmat(i_chan, 1, nep_tmp), ind_stim, ind_cond, repmat(i_subj, 1, nep_tmp), repmat(3, 1, nep_tmp) )) = mampl_P2+mampl_N2 ;
            
            %%% idem for peak amplitudes 
            y_tmp_I1 = y_tmp(:,s_tmp==1) ; 
            ampl_min = min(y_tmp_I1(ind_AUC_I1(1,1):ind_AUC_I1(1,2),:), [],1) ; % N2
            ampl_VW(sub2ind([n_chan, n_stim_test, n_conds_max,n_subj, 3], ...
                repmat(i_chan, 1, nep_tmp), ind_stim, ind_cond, repmat(i_subj, 1, nep_tmp), repmat(1, 1, nep_tmp) )) = ampl_min ;   
            ampl_max = max(y_tmp_I1(ind_AUC_I1(2,1):ind_AUC_I1(2,2),:), [],1) ; % P2
            ampl_VW(sub2ind([n_chan, n_stim_test, n_conds_max,n_subj, 3], ...
                repmat(i_chan, 1, nep_tmp), ind_stim, ind_cond, repmat(i_subj, 1, nep_tmp), repmat(2, 1, nep_tmp) )) = ampl_max ;  
            ampl_VW(sub2ind([n_chan, n_stim_test, n_conds_max,n_subj, 3], ...
                repmat(i_chan, 1, nep_tmp), ind_stim, ind_cond, repmat(i_subj, 1, nep_tmp), repmat(3, 1, nep_tmp) )) = ampl_max-ampl_min ; 
            
            % ======================================== I2
            y_tmp_I2 = abs(y_tmp(:,s_tmp==2)) ; 
            nep_tmp = size(y_tmp_I2,2) ; 
            [ind_stim, ind_cond] = ind2sub([n_stim_test, n_conds_max], lin_idx_I2) ; 
            
            % AUC of N2 -- fills AUC_VW(i_chan, :, :, i_subj, 1)
            AUC_tmp = trapz(y_tmp_I2(ind_AUC_I2(1,1):ind_AUC_I2(1,2),:)) ; 
            AUC_VW(sub2ind([n_chan, n_stim_test, n_conds_max,n_subj, 3], ...
                repmat(i_chan, 1, nep_tmp), ind_stim, ind_cond, repmat(i_subj, 1, nep_tmp), repmat(1, 1, nep_tmp) )) = AUC_tmp ; 
            % AUC of P2 -- AUC_VW(i_chan, :, :, i_subj, 2)
            AUC_tmp = trapz(y_tmp_I2(ind_AUC_I2(2,1):ind_AUC_I2(2,2),:)) ; 
            AUC_VW(sub2ind([n_chan, n_stim_test, n_conds_max,n_subj, 3], ...
                repmat(i_chan, 1, nep_tmp), ind_stim, ind_cond, repmat(i_subj, 1, nep_tmp), repmat(2, 1, nep_tmp) )) = AUC_tmp ; 
            % AUC of N2-P2 -- AUC_VW(i_chan, :, :, i_subj, 3)
            AUC_tmp = trapz(y_tmp_I2(ind_AUC_I2(3,1):ind_AUC_I2(3,2),:)) ; 
            AUC_VW(sub2ind([n_chan, n_stim_test, n_conds_max,n_subj, 3], ...
                repmat(i_chan, 1, nep_tmp), ind_stim, ind_cond, repmat(i_subj, 1, nep_tmp), repmat(3, 1, nep_tmp) )) = AUC_tmp ; 
            %%% idem for mean amplitudes 
            mampl_N2 = mean(y_tmp_I2(ind_AUC_I2(1,1):ind_AUC_I2(1,2),:),1) ; % N2
            mampl_VW(sub2ind([n_chan, n_stim_test, n_conds_max,n_subj, 3], ...
                repmat(i_chan, 1, nep_tmp), ind_stim, ind_cond, repmat(i_subj, 1, nep_tmp), repmat(1, 1, nep_tmp) )) = mampl_N2 ;   
            mampl_P2 = mean(y_tmp_I2(ind_AUC_I2(2,1):ind_AUC_I2(2,2),:),1) ; % P2
            mampl_VW(sub2ind([n_chan, n_stim_test, n_conds_max,n_subj, 3], ...
                repmat(i_chan, 1, nep_tmp), ind_stim, ind_cond, repmat(i_subj, 1, nep_tmp), repmat(2, 1, nep_tmp) )) = mampl_P2 ;  
            mampl_VW(sub2ind([n_chan, n_stim_test, n_conds_max,n_subj, 3], ...
                repmat(i_chan, 1, nep_tmp), ind_stim, ind_cond, repmat(i_subj, 1, nep_tmp), repmat(3, 1, nep_tmp) )) = mampl_P2+mampl_N2 ;
            %%% idem for peak amplitudes 
            y_tmp_I2 = y_tmp(:,s_tmp==2) ; 
            ampl_min = min(y_tmp_I2(ind_AUC_I2(1,1):ind_AUC_I2(1,2),:), [],1) ; % N2
            ampl_VW(sub2ind([n_chan, n_stim_test, n_conds_max,n_subj, 3], ...
                repmat(i_chan, 1, nep_tmp), ind_stim, ind_cond, repmat(i_subj, 1, nep_tmp), repmat(1, 1, nep_tmp) )) = ampl_min ;   
            ampl_max = max(y_tmp_I2(ind_AUC_I2(2,1):ind_AUC_I2(2,2),:), [],1) ; % P2
            ampl_VW(sub2ind([n_chan, n_stim_test, n_conds_max,n_subj, 3], ...
                repmat(i_chan, 1, nep_tmp), ind_stim, ind_cond, repmat(i_subj, 1, nep_tmp), repmat(2, 1, nep_tmp) )) = ampl_max ;  
            ampl_VW(sub2ind([n_chan, n_stim_test, n_conds_max,n_subj, 3], ...
                repmat(i_chan, 1, nep_tmp), ind_stim, ind_cond, repmat(i_subj, 1, nep_tmp), repmat(3, 1, nep_tmp) )) = ampl_max-ampl_min ; 
        end        
    end
end
if save_avg_eeg
    save(fn_avg_eeg, '-v7.3', 'indices_TP', 'time_vec', 'fs', 'xstep', ...
        'chanlocs', 'all_TPs', 'all_TPs_str', 'y_I1', 'y_I2', ...
        'y_I1_TP', 'y_I2_TP', 'y_I1_cond', 'y_I2_cond') ; 
end
if save_AUC_VW
    fn_AUC = [fn_data_tmp, 'AUC', str_params, str_ampl_max, '.mat'] ; 
    save(fn_AUC, '-v7.3', 'AUC_VW', 'ampl_VW', 'mampl_VW', 'time_vec', 'ind_AUC_I1', 'ind_AUC_I2', 'all_ep_kept') ; 
    disp('************************ AUC SAVED')
end


% ======================================================================= %
% ===== * ===== Save time-freq data (avg over all epochs)
% ======================================================================= %
if save_avg_TF
    % ======================================== filenames
    str_avg_TF = '' ; 
    if use_stft
       str_avg_TF = [str_avg_TF, '_w', num2str(round(100*win_width_sec))] ; % window width
    end
    allchans_TF = 1 ;     
    if allchans_TF
        all_chans_TF = 1:n_chan ; 
        str_chans = '' ; 
    else
        chan_name = 'FCZ' ; 
        [~, chan_idx] = define_chan_order_colors(chanlocs, {chan_name}) ;    
        all_chans_TF = chan_idx ; % 1:n_chan ; 
        str_chans = ['_', chan_name] ; 
    end
    fn_avg_TF = [fn_data_tmp, 'TF', str_params, str_ampl_max, str_avg_TF, str_chans, '.mat'] 
    % ========================================
    
    nchan_TF = length(all_chans_TF) ;     
    % Down-sampling
    df_bTF = 2 ;    % before TF --> 250Hz = 4ms step
    df_aTF = 10 ;   % after TF --> 25 Hz = 40ms step (500/20)
    t_min = -0.5 ; t_max = 1 ;
    time_bTF = downsample(time_vec, df_bTF) ; nT_bTF = length(time_bTF) ;
    %t_min = time_TF(1) ; t_max = time_TF(end) ;     
    time_TF = downsample(time_bTF, df_aTF) ; 
    % str_TF
    nT_TF = length(time_TF) ; 
    freqs_tmp = f_begin:f_step:f_stop ; 
    nF_TF = length(freqs_tmp) ;     
    
    y_TF_I1 = NaN(n_chan, nF_TF, nT_TF, n_subj) ;
    y_TF_I2 = NaN(n_chan, nF_TF, nT_TF, n_subj) ;
    % idem but using only data during first testing bloc
    y_TF_I1_c1 = NaN(n_chan, nF_TF, nT_TF, n_subj) ;
    y_TF_I2_c1 = NaN(n_chan, nF_TF, nT_TF, n_subj) ;
    
    %%%%%%%%%%%%%% Extract pre-alpha power
    time_alpha = [-0.42, -0.06] ;
    [~,alpha_tmin] = min(abs(time_alpha(1)-time_TF)) ;  
    [~,alpha_tmax] = min(abs(time_alpha(2)-time_TF)) ;    
    freq_alpha = [8, 12] ; 
    [~,alpha_fmin] = min(abs(freq_alpha(1)-freqs_tmp)) ;  
    [~,alpha_fmax] = min(abs(freq_alpha(2)-freqs_tmp)) ; 
    pre_alpha = NaN(n_chan, n_stim_test, n_conds_max,n_subj) ;
    fn_prealpha = [fn_data_tmp, 'PreA', str_params, str_ampl_max, str_avg_TF, str_chans, '.mat'] ; 
    
    
    ts_TF = tic ; 
    for i_subj=1:n_subj
        for i_chan=all_chans_TF           
            y_tmp = downsample(squeeze(all_EEG_data(i_chan, :, :, :, i_subj)), df_bTF) ; 
            % (n_times_TF, n_stim_test, n_conds_max)
            s_tmp = squeeze(all_s(:,:,i_subj)) ; 
            % (n_stim_test, n_conds_max)
            
            % ==*== Merging all conditions
            y_tmp = reshape(y_tmp,[nT_bTF,N_obs]) ; 
            s_tmp = reshape(s_tmp,[1,N_obs]) ; 
            cond_tmp = reshape(ones(n_stim_test,1)*[1:n_conds_max], [1, N_obs]) ; 

            if reject_ep
                ep_kept_tmp = reshape(squeeze(all_ep_kept(:,:,i_subj)), ...
                    [1, N_obs]) ; 
                % all_ep_kept = NaN(n_stim_test, n_conds_max, n_subj) ; 
                y_tmp = y_tmp(:,ep_kept_tmp) ; 
                
                lin_idx_I1 = find(ep_kept_tmp & s_tmp==1) ;
                lin_idx_I2 = find(ep_kept_tmp & s_tmp==2) ;
                s_tmp = s_tmp(ep_kept_tmp) ; 
                cond_tmp = cond_tmp(ep_kept_tmp) ; 
                cond_tmp = cond_tmp==1 ; 
            end
            
            %%%%% I1
            y_tmp_I1 = (y_tmp(:,s_tmp==1))' ; % n_ep, n_times_TF
            [all_STFT, freqs_TF, ~,~,~,~] = ...
                extract_time_freq_signal(y_tmp_I1, fs/df_bTF, time_bTF, ...
                'f_begin', f_begin, 'f_stop', f_stop, 'f_step', ...
                f_step,'win_width_sec', win_width_sec, 'stand_STFT', ...
                0, 'use_cwt', use_cwt, 'use_stft', use_stft, ...
                'max_filt_order', 100, 'f_half_width', f_step/2, ...
                'stft_with_conv', stft_with_conv) ; % [nF, nT, n_epoch]
            all_STFT = downsample(permute(abs(all_STFT),[2, 1, 3]),df_aTF) ; % (nT, nF, n_ep)
            all_STFT = permute(all_STFT, [3,2,1]) ; % (n_ep, nF, nT)
            y_TF_I1(i_chan, :,:, i_subj) = squeeze(mean(all_STFT,1)) ; 
            % idem using only epochs from boc 1
            y_TF_I1_c1(i_chan, :,:, i_subj) = squeeze(mean(all_STFT(cond_tmp(s_tmp==1),:,:),1)) ; 
            
            nep_tmp = size(all_STFT,1) ; 
            alpha_tmp = mean(mean(all_STFT(:,alpha_fmin:alpha_fmax,alpha_tmin:alpha_tmax), 2), 3) ; 
            [ind_stim, ind_cond] = ind2sub([n_stim_test, n_conds_max], lin_idx_I1) ; 
            pre_alpha(sub2ind([n_chan, n_stim_test, n_conds_max,n_subj], ...
                repmat(i_chan, 1, nep_tmp), ind_stim, ind_cond, repmat(i_subj, 1, nep_tmp))) = alpha_tmp ; 
            
            %%%%% I2
            y_tmp_I2 = (y_tmp(:,s_tmp==2))' ; % n_ep, n_times_TF
            [all_STFT, freqs_TF, ~,~,~,~] = ...
                extract_time_freq_signal(y_tmp_I2, fs/df_bTF, time_bTF, ...
                'f_begin', f_begin, 'f_stop', f_stop, 'f_step', ...
                f_step,'win_width_sec', win_width_sec, 'stand_STFT', ...
                0, 'use_cwt', use_cwt, 'use_stft', use_stft, ...
                'max_filt_order', 100, 'f_half_width', f_step/2, ...
                'stft_with_conv', stft_with_conv) ; % [nF, nT, n_epoch]
            all_STFT = downsample(permute(abs(all_STFT),[2, 1, 3]),df_aTF) ; 
            all_STFT = permute(all_STFT, [3,2,1]) ; % (n_ep, nF, nT)
            y_TF_I2(i_chan, :,:, i_subj) = squeeze(mean(all_STFT,1)) ; 
            % idem using only epochs from boc 1
            y_TF_I2_c1(i_chan, :,:, i_subj) = squeeze(mean(all_STFT(cond_tmp(s_tmp==2),:,:),1)) ; 
            
            nep_tmp = size(all_STFT,1) ; 
            alpha_tmp = mean(mean(all_STFT(:,alpha_fmin:alpha_fmax,alpha_tmin:alpha_tmax), 2), 3) ; 
            [ind_stim, ind_cond] = ind2sub([n_stim_test, n_conds_max], lin_idx_I2) ; 
            pre_alpha(sub2ind([n_chan, n_stim_test, n_conds_max,n_subj], ...
                repmat(i_chan, 1, nep_tmp), ind_stim, ind_cond, repmat(i_subj, 1, nep_tmp))) = alpha_tmp ; 
        end   
        disp([num2str(toc(ts_TF)), ' sec for ', num2str(nchan_TF), ' chans and ', ...
            num2str(i_subj), ' subjs'])
    end
    
    save(fn_avg_TF, '-v7.3', 'y_TF_I1', 'y_TF_I2', 'time_TF', 'freqs_TF', ...
        'f_begin', 'f_step', 'f_stop', 'y_TF_I1_c1', 'y_TF_I2_c1') ; 
    % Reload & plot in TSL_plot_avg_EEG.m
    fn_prealpha = fn_prealpha
    save(fn_prealpha, '-v7.3', 'all_trials_resp', 'all_s', 'pre_alpha', 'all_ep_kept', ...
        'indices_TP', 'all_TPs', 'all_TPs_str')
    % Reload and plot in TSL_corr_TF_conf.m
    return
end
if save_avg_eeg ; return ; end



if plot_avg_eeg
    % ======================================================================= %
    % ===== * ===== Plot the average EEG
    % ======================================================================= %
    chan_topo = logical([chanlocs.topo_enabled]) ;
    chanlocs_butterfly = chanlocs ; 

    y_I1_avg = squeeze(mean(y_I1,3)) ; 
    y_I2_avg = squeeze(mean(y_I2,3)) ; 

    if selected_chan
        % only consider chan_selected for the plots
        [~, chan_desired] = define_chan_order_colors(chanlocs, chan_selected) ; 
        y_I1_avg = y_I1_avg(chan_desired, :) ; 
        y_I2_avg = y_I2_avg(chan_desired, :) ; 
        chanlocs_butterfly = chanlocs(chan_desired) ; 

    elseif only_plot_topo_chan
        y_I1_avg = y_I1_avg(chan_topo, :) ; 
        y_I2_avg = y_I2_avg(chan_topo, :) ; 
        if topo_colors
            chanlocs_butterfly = chanlocs(chan_topo) ; 
        end
    end


    % in sec
    xlim_butterfly = [-0.25,0.75] ; % [time_vec(1), time_vec(end)] ; 
    time_vec_used = time_vec ; 

    if plot_ms
        % in seconds
        xlim_butterfly = 1000.*xlim_butterfly ;
        time_vec_used = 1000.*time_vec_used ; 
    end

    % ------------------------------------------------------------------- %
    % --*-- I_1
    % ------------------------------------------------------------------- %
    if selected_chan
        fn_res_I1 = [fn_res, num2str(n_chans_plot), 'chans_I1/'] ; 
        if ~exist(fn_res_I1, 'dir') ;  mkdir(fn_res_I1) ; end
    else
        fn_res_I1 = fn_res ; 
    end

    fn_I1 = [fn_res_I1, 'Butterfly_I1', str_params,str_plot] ; 
    plot_butterfly(time_vec_used, y_I1_avg, 'fn_save', fn_I1, 'xlim_val', ...
        xlim_butterfly, 'x_lab', 'time after stimulus (ms)', 'topo_colors', ...
        topo_colors, 'chanlocs', chanlocs_butterfly, 'eps_fig', eps_fig, ...
        'fig_fig', fig_fig, 'pdf_fig', pdf_fig) ; 

    % ------------------------------------------------------------------- %
    % --*-- I_2
    % ------------------------------------------------------------------- %
    if selected_chan
        fn_res_I2 = [fn_res, num2str(n_chans_plot), 'chans_I2/'] ; 
        if ~exist(fn_res_I2, 'dir') ;  mkdir(fn_res_I2) ; end
    else
        fn_res_I2 = fn_res ; 
    end

    fn_I2 = [fn_res_I2, 'Butterfly_I2', str_params,str_plot] ; 
    plot_butterfly(time_vec_used, y_I2_avg, 'fn_save', fn_I2, 'xlim_val', ...
        xlim_butterfly, 'x_lab', 'time after stimulus (ms)', 'topo_colors', ...
        topo_colors, 'chanlocs', chanlocs_butterfly, 'eps_fig', eps_fig, ...
        'fig_fig', fig_fig, 'pdf_fig', pdf_fig) ; 
end

% ======================================================================= %
% ===== * ===== IO fitting
% ======================================================================= %

%%%%%%%%% For all IO fittings
down_factor = 2 ; 
t_min = -0.5 ; 
t_max = 1 ; 

time_vec_IO = downsample(time_vec,down_factor) ;
[~, idx_m] = min(abs(time_vec_IO-t_min)) ; 
[~, idx_M] = min(abs(time_vec_IO-t_max)) ; 
time_vec_IO = time_vec_IO(idx_m:idx_M) ; 
n_times_IO = length(time_vec_IO) ;

%%%%%%%%% Time-frequency maps
downf_TF = 10 ; % in addition to down_factor, after TF
% --> 25 Hz = 40ms step
time_vec_TF = downsample(time_vec_IO,downf_TF) ; 


nT_TF = length(time_vec_TF) ; 
nF_TF = length(f_begin:f_step:f_stop) ; 

%%%%%%%%%%%%%%%%%%%%%%%%%% Principle:
% loop on all time steps
%   regress the EEG on the IO data (surprise, conf, ...)
%   save R^2, MSE, ... as a fct of time, channel, subjects, params
%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Do it separately when all params vs. when a single param
%%% Because 
%%% * can test variants in single param case (eg residual conf vs. conf)
%%% * should limit nb of computations in all params case

if single_elec 
    [~, chan_desired] = define_chan_order_colors(chanlocs, {elec_IO}) ;
    all_chan_idx = 1 ; 
    n_chan = 1 ; 
    all_EEG_data = all_EEG_data(chan_desired, :, :, :, :) ; 
else
    all_chan_idx = 1:n_chan ; 
end

if TF_regress
    nF_data =  nF_TF ;
    nT_data = nT_TF ;
    all_TFs = zeros(n_chan, nF_data, nT_data, 2) ; % average TF data for I1 and I2
else
    nF_data =  1 ; % baseband response
    nT_data = n_times_IO ;
end
alpha_level = 0.05 ; paired_tests = 0 ; mean_H0 = 0 ; always_t_test = 1 ;
two_sided = 1 ; all_names = {chanlocs.labels} ;  


switch IO_fit_opt
    case 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% Single parameter
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %all_IO_data: 
        % 6, n_stim_test, n_conds_max, n_subj
        % dim1: surprise - confidence - update - unpredictability - p1 -
        % pred error
        
        % Outcomes of the fitting (1 regression per IO_outcome)
        % - 1 regression with surprise and residual conf as predictors --> 2 betas
        % - 1 regression with surprise and conf as predictors --> 2 betas
        % /!\ separately for I_1 and I_2 --> x 2
        % if p_predictor --> add the predicted proba as predictor!
        if p_predictor
            n_betas = 5 ; 
        else
            n_betas = 4 ; 
        end
        all_beta = NaN(n_betas, n_chan, nF_data, nT_data, n_subj, 2) ; 
        all_t_beta = zeros(n_betas, n_chan, nF_data, nT_data, 2) ; 
        all_p_beta = ones(n_betas, n_chan, nF_data, nT_data, 2) ; 

        all_R_squared = NaN(2, n_chan, nF_data, nT_data, n_subj, 2) ; 
        all_Rs_VIF = NaN(2, n_chan, n_subj, 2) ; % btw (r)conf & PE
        % https://www.geeksforgeeks.org/detecting-multicollinearity-with-vif-python/
        all_Rs_IV = NaN(n_betas, n_chan, nF_data, nT_data, n_subj, 2) ; % for each indep variable separately
        all_MSE = NaN(2, n_chan, nF_data, nT_data, n_subj, 2) ; 
        
        n_blocs = 1 ; % repeat IO_fitting n_blocs times
        n_perm_bloc = n_perm ;
        
        % ----- avg EEG with median split on PE and rconf
        % Order:
        labels_split = {'low rc - high PE', 'low rc - low PE', ...
            'high rc - high PE', 'high rc - low PE', ...
            'low rc', 'high rc', 'high PE', 'low PE'} ; 
        y_I1_split = NaN(n_chan, nF_data, nT_data, n_subj, 8) ;
        y_I2_split = NaN(n_chan, nF_data, nT_data, n_subj, 8) ; 
        N_I1_split = zeros(n_subj, 8) ; % # of stim in each category 
        N_I2_split = zeros(n_subj, 8) ; 
        
        if CB_PT
            if ~single_elec                 
                % Compute n_perm_bloc at once, otherwise beta_CP too large
                n_blocs = 10 ;
                n_perm_bloc = n_perm/n_blocs
            end
            beta_CP = NaN(n_betas, n_chan, nF_data, nT_data, n_subj, 2, n_perm_bloc) ;
            perm_cluster_t = zeros(n_perm, n_betas, 2, 2) ; 
        end

        ts = tic ; 
        for i_bloc=1:n_blocs
            for i_subj=1:n_subj
                IO_data_tmp = squeeze(all_IO_data(:, :, :, i_subj)) ;
                % dim1: surprise - confidence - update - unpredictability - p1 -
                % prederror
                % x n_stim_test, n_conds_max, n_subj

                s_tmp = squeeze(all_s(:,:,i_subj)) ;    % (n_stim_test, n_conds_max)

                cond_vec = indices_TP(:,i_subj) ;       % (n_cond_max, 1)
                cond_vec = repmat(cond_vec', n_stim_test,1) ; 

                % ==*== Merging all conditions
                s_tmp = reshape(s_tmp,[1,N_obs]) ;
                IO_data_tmp = (reshape(IO_data_tmp, [n_IO_out, N_obs]))' ; 
                cond_vec = (reshape(cond_vec, [1, N_obs]))' ; 

                if reject_ep
                    ep_kept_tmp = reshape(squeeze(all_ep_kept(:,:,i_subj)), ...
                        [1, N_obs]) ; 
                    % all_ep_kept = NaN(n_stim_test, n_conds_max, n_subj) ; 
                    s_tmp = s_tmp(ep_kept_tmp) ; 
                    IO_data_tmp = IO_data_tmp(ep_kept_tmp,:) ; 
                    cond_vec = cond_vec(ep_kept_tmp) ;
                end


                for i_chan=all_chan_idx
                    y_tmp = downsample(squeeze(all_EEG_data(i_chan, :, :, :, i_subj)), down_factor) ; 
                    y_tmp = y_tmp(idx_m:idx_M,:,:) ; 
                    % (n_times, n_stim_test, n_conds_max)

                    % ==*== Merging all conditions
                    y_tmp = reshape(y_tmp,[n_times_IO,N_obs]) ;
                    if reject_ep
                        y_tmp = y_tmp(:,ep_kept_tmp) ; 
                    end

                    % =========================================================== %
                    % ===*=== I_1
                    % [beta, mse, Rs] = regress_Y_on_X(Y,X) ;
                    y_tmp_I1 = (y_tmp(:,s_tmp==1))' ;

                    if TF_regress
                        % y_tmp_I1: (n_epochs, n_times_IO)
                        % all_STFT: complex
                        [all_STFT, all_freq_Hz, all_t_STFT,~,f0,Fb] = ...
                            extract_time_freq_signal(y_tmp_I1, fs/down_factor, time_vec_IO, ...
                            'f_begin', f_begin, 'f_stop', f_stop, 'f_step', ...
                            f_step,'win_width_sec', win_width_sec, 'stand_STFT', ...
                            0, 'use_cwt', use_cwt, 'use_stft', use_stft, ...
                            'max_filt_order', 100, 'f_half_width', f_step/2, ...
                            'stft_with_conv', stft_with_conv) ; % [nF, nT, n_epoch]
                        all_STFT = downsample(permute(abs(all_STFT),[2, 1, 3]),downf_TF) ; % (nT, nF, n_ep)
                        all_STFT = permute(all_STFT, [3,2,1]) ; % (n_ep, nF, nT)
                        y_tmp_I1 = reshape(all_STFT,size(y_tmp_I1,1), nF_TF*nT_TF) ; 
                        % Store TF data
                        all_TFs(i_chan, :,:,1) = all_TFs(i_chan, :,:,1) + ...
                            (mean(all_STFT,1))./n_subj ; 
                        %zeros(n_chan, nF_data, nT_data, 2) ;
                    end


                    IO_data_I1 = IO_data_tmp(s_tmp==1,:) ; 
                    % columns: surprise - confidence - update - unpredictability - p1 -
                    % PE
                    IO_p1_I1 = IO_data_I1(:,5) ; % attention: need to extract p1_mean BEFORE z-score, because log2 later
                    IO_N_I1 = IO_data_I1(:,7) ; % N_used
                    IO_data_I1 = zscore(IO_data_I1) ; 
                    cond_vec_I1 = cond_vec(s_tmp==1) ; 

                    % -----> compute RESIDUAL confidence
                    if rconf_no_logN
                        [~, ~, ~, eps_c] = regress_Y_on_X(IO_data_I1(:,2), ...
                            [IO_p1_I1, IO_p1_I1.^2, log2(IO_p1_I1), log2(1-IO_p1_I1), ...
                            log2(IO_N_I1)], random_inter, cond_vec_I1) ;
                    else
                        [~, ~, ~, eps_c] = regress_Y_on_X(IO_data_I1(:,2), ...
                            [IO_p1_I1, IO_p1_I1.^2, log2(IO_p1_I1), log2(1-IO_p1_I1)], random_inter, cond_vec_I1) ;
                    end
                    
                    if i_bloc==1
                        % ----- avg EEG with median split on surprise and rconf
                        % labels_split = {'low rc - high s', 'low rc - low s', ...
                        % 'high rc - high s', 'high rc - low s', ...
                        % 'low rc', 'high rc', 'high s', 'low s'} ;                         
                        y_used = y_tmp_I1 - mean(y_tmp_I1(:)) ; % remove subject mean  
                        if rconf_split
                            conf_used = eps_c ; 
                        else
                            conf_used = IO_data_I1(:,2) ; % eps_c
                        end
                        %idx_lc = eps_c<median(eps_c) ; 

                        idx_lc = conf_used<median(conf_used) ;
                        idx_hc = conf_used>=median(conf_used) ; 
                        
                        idx_hs = IO_data_I1(:,1)>=median(IO_data_I1(:,1)) ; 
                        idx_ls = IO_data_I1(:,1)<median(IO_data_I1(:,1)) ; 
                        % ============== Split by conf OR surp (sep)
                        % ** low c
                        y_lc = y_used(idx_lc, :) ; 
                        y_I1_split(i_chan, :,:,i_subj,5) = ...
                            reshape(mean(y_lc,1), nF_data, nT_data);
                        N_I1_split(i_subj,5) = sum(idx_lc) ;
                        % ** high c
                        y_hc = y_used(idx_hc, :) ;
                        y_I1_split(i_chan, :,:,i_subj,6) = ...
                            reshape(mean(y_hc,1), nF_data, nT_data);
                        N_I1_split(i_subj,6) = sum(idx_hc) ;
                        % ** high s
                        y_I1_split(i_chan, :,:,i_subj,7) = ...
                            reshape(mean(y_used(idx_hs, :),1), ...
                            nF_data, nT_data); % high s
                        N_I1_split(i_subj,7) = sum(idx_hs) ;                        
                        % ** low s                        
                        y_I1_split(i_chan, :,:,i_subj,8) = ...
                            reshape(mean(y_used(idx_ls, :),1), ...
                            nF_data, nT_data); % low s
                        N_I1_split(i_subj,8) = sum(idx_ls) ;
                        
                        % ============== Split by conf AND surp                        
                        % ** low c, high s  
                        %s_lc = IO_data_I1(idx_lc,1) ; % to keep same #obs in 4 categories
                        %idx_hs_lc = s_lc>=median(s_lc) ; idx_ls_lc = s_lc<median(s_lc) ; 
                        idx_hs_lc = idx_hs & idx_lc ; 
                        y_I1_split(i_chan, :,:,i_subj,1) = ...
                            reshape(mean(y_used(idx_hs_lc, :),1), nF_data, nT_data);
                        N_I1_split(i_subj,1) = sum(idx_hs_lc) ;
                        % ** low c, low s  
                        idx_ls_lc = idx_ls & idx_lc ;
                        y_I1_split(i_chan, :,:,i_subj,2) = ...
                            reshape(mean(y_used(idx_ls_lc, :),1), nF_data, nT_data);
                        N_I1_split(i_subj,2) = sum(idx_ls_lc) ; 
                        % ** high c, high s 
                        %s_hc = IO_data_I1(idx_hc,1) ; 
                        %idx_hs_hc = s_hc>=median(s_hc) ; idx_ls_hc = s_hc<median(s_hc) ; 
                        idx_hs_hc = idx_hs & idx_hc ; 
                        y_I1_split(i_chan, :,:,i_subj,3) = ...
                            reshape(mean(y_used(idx_hs_hc, :),1), nF_data, nT_data);
                        N_I1_split(i_subj,3) = sum(idx_hs_hc) ;
                        % ** high c, low s 
                        idx_ls_hc = idx_ls & idx_hc ; 
                        y_I1_split(i_chan, :,:,i_subj,4) = ...
                            reshape(mean(y_used(idx_ls_hc, :),1), nF_data, nT_data);
                        N_I1_split(i_subj,4) = sum(idx_ls_hc) ;
                    end

                    % -----> surprise + RESIDUAL confidence (+ p1 if
                    % p_predictor)
                    % ------------------------------------------------------- %
                    if p_predictor
                        % select current predictors used
                        pred_tmp = [IO_data_I1(:,1),zscore(eps_c), IO_data_I1(:,5)] ; 
                    else
                        pred_tmp = [IO_data_I1(:,1), zscore(eps_c)] ; 
                    end
                    if i_bloc==1
                    [beta, mse, Rs] = regress_Y_on_X(y_tmp_I1, pred_tmp, random_inter, cond_vec_I1) ;

                    all_beta(1:(n_betas-2),i_chan, :, :, i_subj, 1) = reshape(beta, n_betas-2, nF_data, nT_data); 
                    all_R_squared(1,i_chan, :, :, i_subj, 1) = reshape(Rs, nF_data, nT_data); 
                    all_MSE(1, i_chan,:, :, i_subj, 1) = reshape(mse, nF_data, nT_data) ;
                    
                    % Regress rconf on PE to measure collinearity
                    [~, ~, Rs] = regress_Y_on_X(zscore(eps_c), IO_data_I1(:,1), random_inter, cond_vec_I1) ;
                    all_Rs_VIF(1, i_chan, i_subj, 1) = Rs ; % btw (r)conf & PE
                    
                    npred_first = size(pred_tmp,2) ; 
                    for ipred=1:npred_first
                    corr_tmp = corr(pred_tmp(:,ipred), y_tmp_I1) ; % (1, nF*NT) matrix
                    all_Rs_IV(ipred, i_chan,:,:, i_subj,1) = reshape(corr_tmp.^2, nF_data, nT_data) ;
                    %(n_betas, n_chan, nF_data, nT_data, n_subj, 2)
                    end
                    end
                    % ============================= %
                    if CB_PT
                        n_ep_tmp = size(IO_data_I1,1) ; 
                        for i_perm=1:n_perm_bloc
                            beta = regress_Y_on_X(y_tmp_I1, pred_tmp(randperm(n_ep_tmp),:), random_inter, cond_vec_I1) ;
                            beta_CP(1:(n_betas-2),i_chan, :, :, i_subj, 1, i_perm) = ...
                                reshape(beta, n_betas-2, nF_data, nT_data) ; 
                        end
                    end
                    % ============================= % 

                    % -----> surprise + conf
                    % ------------------------------------------------------- %
                    pred_tmp = IO_data_I1(:,[1:2]) ; 
                    if i_bloc==1
                    [beta, mse, Rs] = regress_Y_on_X(y_tmp_I1, pred_tmp, random_inter, cond_vec_I1) ;

                    all_beta(n_betas-1:n_betas,i_chan, :, :, i_subj, 1) = reshape(beta, 2, nF_data, nT_data); 
                    all_R_squared(2,i_chan, :, :, i_subj, 1) = reshape(Rs, nF_data, nT_data); 
                    all_MSE(2, i_chan,:, :, i_subj, 1) = reshape(mse, nF_data, nT_data) ; 
                    
                    % Regress conf on PE to measure collinearity
                    [~, ~, Rs] = regress_Y_on_X(IO_data_I1(:,2), IO_data_I1(:,1), random_inter, cond_vec_I1) ;
                    all_Rs_VIF(2, i_chan, i_subj, 1) = Rs ; % btw (r)conf & PE
                    for ipred=1:size(pred_tmp,2)
                    corr_tmp = corr(pred_tmp(:,ipred), y_tmp_I1) ; % (1, nF*NT) matrix
                    all_Rs_IV(npred_first+ipred, i_chan,:,:, i_subj,1) = reshape(corr_tmp.^2, nF_data, nT_data) ;
                    %(n_betas, n_chan, nF_data, nT_data, n_subj, 2)
                    end
                    end
                    % ============================= %
                    if CB_PT
                        for i_perm=1:n_perm_bloc
                            beta = regress_Y_on_X(y_tmp_I1, pred_tmp(randperm(n_ep_tmp),:), random_inter, cond_vec_I1) ;
                            beta_CP(n_betas-1:n_betas,i_chan, :, :, i_subj, 1, i_perm) = ...
                                reshape(beta, 2, nF_data, nT_data) ; 
                        end
                    end
                    % ============================= % 

                    % =========================================================== %
                    % ===*=== I_2
                    y_tmp_I2 = (y_tmp(:,s_tmp==2))' ;
                    if TF_regress
                        % y_tmp_I2: (n_epochs, n_times_IO)
                        all_STFT = ...
                            extract_time_freq_signal(y_tmp_I2, fs/down_factor, time_vec_IO, ...
                            'f_begin', f_begin, 'f_stop', f_stop, 'f_step', ...
                            f_step,'win_width_sec', win_width_sec, 'stand_STFT', ...
                            0, 'use_cwt', use_cwt, 'use_stft', use_stft, ...
                            'max_filt_order', 100, 'f_half_width', f_step/2, ...
                            'stft_with_conv', stft_with_conv) ; % [nF, nT, n_epoch]
                        all_STFT = downsample(permute(abs(all_STFT),[2, 1, 3]),downf_TF) ; % (nT, nF, n_ep)
                        all_STFT = permute(all_STFT, [3,2,1]) ; % (n_ep, nF, nT)
                        y_tmp_I2 = reshape(all_STFT,size(y_tmp_I2,1), nF_TF*nT_TF) ; 
                        all_TFs(i_chan, :,:,2) = all_TFs(i_chan, :,:,2) + ...
                            (mean(all_STFT,1))./n_subj ; 
                    end

                    IO_data_I2 = IO_data_tmp(s_tmp==2,:) ;
                    IO_p1_I2 = IO_data_I2(:,5) ;
                    IO_N_I2 = IO_data_I2(:,7) ; % N_used
                    IO_data_I2 = zscore(IO_data_I2) ; 
                    cond_vec_I2 = cond_vec(s_tmp==2) ; 

                    % -----> compute RESIDUAL confidence
                    if rconf_no_logN
                        [~, ~, ~, eps_c] = regress_Y_on_X(IO_data_I2(:,2), ...
                            [IO_p1_I2, IO_p1_I2.^2, log2(IO_p1_I2), log2(1-IO_p1_I2), ...
                            log2(IO_N_I2)], random_inter, cond_vec_I2) ;
                    else
                        [~, ~, ~, eps_c] = regress_Y_on_X(IO_data_I2(:,2), ...
                            [IO_p1_I2, IO_p1_I2.^2, log2(IO_p1_I2), log2(1-IO_p1_I2), ...
                            ], random_inter, cond_vec_I2) ;
                    end
                    
                    if i_bloc==1
                        % ----- avg EEG with median split on surprise and rconf
                        % labels_split = {'low rc - high s', 'low rc - low s', ...
                        % 'high rc - high s', 'high rc - low s', ...
                        % 'low rc', 'high rc', 'high s', 'low s'} ; 
                        
                        y_used = y_tmp_I2 - mean(y_tmp_I2(:)) ; % remove subject mean  
                        if rconf_split
                            conf_used = eps_c ; 
                        else
                            conf_used = IO_data_I2(:,2) ; % eps_c
                        end
                        %idx_lc = eps_c<median(eps_c) ; 

                        idx_lc = conf_used<median(conf_used) ;
                        idx_hc = conf_used>=median(conf_used) ; 
                        
                        idx_hs = IO_data_I2(:,1)>=median(IO_data_I2(:,1)) ; 
                        idx_ls = IO_data_I2(:,1)<median(IO_data_I2(:,1)) ; 
                        % ============== Split by conf OR surp (sep)
                        % ** low c
                        y_lc = y_used(idx_lc, :) ; 
                        y_I2_split(i_chan, :,:,i_subj,5) = ...
                            reshape(mean(y_lc,1), nF_data, nT_data);
                        N_I2_split(i_subj,5) = sum(idx_lc) ;
                        % ** high c
                        y_hc = y_used(idx_hc, :) ;
                        y_I2_split(i_chan, :,:,i_subj,6) = ...
                            reshape(mean(y_hc,1), nF_data, nT_data);
                        N_I2_split(i_subj,6) = sum(idx_hc) ;
                        % ** high s
                        y_I2_split(i_chan, :,:,i_subj,7) = ...
                            reshape(mean(y_used(idx_hs, :),1), ...
                            nF_data, nT_data); % high s
                        N_I2_split(i_subj,7) = sum(idx_hs) ;                        
                        % ** low s                        
                        y_I2_split(i_chan, :,:,i_subj,8) = ...
                            reshape(mean(y_used(idx_ls, :),1), ...
                            nF_data, nT_data); % low s
                        N_I2_split(i_subj,8) = sum(idx_ls) ;
                        
                        % ============== Split by conf AND surp (sep)                     
                        % ** low c, high s  
                        %s_lc = IO_data_I2(idx_lc,1) ; 
                        %idx_hs_lc = s_lc>=median(s_lc) ; idx_ls_lc = s_lc<median(s_lc) ; 
                        idx_hs_lc = idx_hs & idx_lc ; 
                        y_I2_split(i_chan, :,:,i_subj,1) = ...
                            reshape(mean(y_used(idx_hs_lc, :),1), nF_data, nT_data);
                        N_I2_split(i_subj,1) = sum(idx_hs_lc) ;
                        % ** low c, low s  
                        idx_ls_lc = idx_ls & idx_lc ;
                        y_I2_split(i_chan, :,:,i_subj,2) = ...
                            reshape(mean(y_used(idx_ls_lc, :),1), nF_data, nT_data);
                        N_I2_split(i_subj,2) = sum(idx_ls_lc) ; 
                        % ** high c, high s 
                        %s_hc = IO_data_I2(idx_hc,1) ; 
                        %idx_hs_hc = s_hc>=median(s_hc) ; idx_ls_hc = s_hc<median(s_hc) ; 
                        idx_hs_hc = idx_hs & idx_hc ; 
                        y_I2_split(i_chan, :,:,i_subj,3) = ...
                            reshape(mean(y_used(idx_hs_hc, :),1), nF_data, nT_data);
                        N_I2_split(i_subj,3) = sum(idx_hs_hc) ;
                        % ** high c, low s 
                        idx_ls_hc = idx_ls & idx_hc ; 
                        y_I2_split(i_chan, :,:,i_subj,4) = ...
                            reshape(mean(y_used(idx_ls_hc, :),1), nF_data, nT_data);
                        N_I2_split(i_subj,4) = sum(idx_ls_hc) ;
                       
                    end
                    
                    % -----> surprise + RESIDUAL confidence
                    % ------------------------------------------------------- %
                    if p_predictor
                        pred_tmp = [IO_data_I2(:,1), zscore(eps_c), zscore(1-IO_p1_I2)] ; 
                    else
                        pred_tmp = [IO_data_I2(:,1),zscore(eps_c)] ; 
                    end    
                    if i_bloc==1
                    [beta, mse, Rs] = regress_Y_on_X(y_tmp_I2, pred_tmp, random_inter, cond_vec_I2) ;
                    all_beta(1:(n_betas-2),i_chan, :, :, i_subj, 2) = reshape(beta, n_betas-2, nF_data, nT_data); 
                    all_R_squared(1,i_chan, :, :, i_subj, 2) = reshape(Rs, nF_data, nT_data); 
                    all_MSE(1, i_chan,:, :, i_subj, 2) = reshape(mse, nF_data, nT_data) ; 
                    
                    % Regress rconf on PE to measure collinearity
                    [~, ~, Rs] = regress_Y_on_X(zscore(eps_c), IO_data_I2(:,1), random_inter, cond_vec_I2) ;
                    all_Rs_VIF(1, i_chan, i_subj, 2) = Rs ;
                    for ipred=1:npred_first
                        corr_tmp = corr(pred_tmp(:,ipred), y_tmp_I2) ; % (1, nF*NT) matrix
                        all_Rs_IV(ipred, i_chan,:,:, i_subj,2) = reshape(corr_tmp.^2, nF_data, nT_data) ;
                    end
                    
                    end
                    % ============================= %
                    if CB_PT
                        n_ep_tmp = size(IO_data_I2,1) ; 
                        for i_perm=1:n_perm_bloc
                            beta = regress_Y_on_X(y_tmp_I2, pred_tmp(randperm(n_ep_tmp),:), random_inter, cond_vec_I2) ;
                            beta_CP(1:(n_betas-2),i_chan, :, :, i_subj, 2, i_perm) = ...
                                reshape(beta, n_betas-2, nF_data, nT_data) ; 
                        end
                    end
                    % ============================= % 

                    % -----> surprise + conf
                    % ------------------------------------------------------- %
                    pred_tmp = IO_data_I2(:,[1:2]) ; 
                    if i_bloc==1                        
                    [beta, mse, Rs] = regress_Y_on_X(y_tmp_I2, pred_tmp, random_inter, cond_vec_I2) ;
                    all_beta(n_betas-1:n_betas,i_chan, :, :, i_subj, 2) = reshape(beta, 2, nF_data, nT_data); 
                    all_R_squared(2,i_chan, :, :, i_subj, 2) = reshape(Rs, nF_data, nT_data); 
                    all_MSE(2, i_chan,:, :, i_subj, 2) = reshape(mse, nF_data, nT_data) ; 
                    
                    % Regress conf on PE to measure collinearity
                    [~, ~, Rs] = regress_Y_on_X(IO_data_I2(:,2), IO_data_I2(:,1), random_inter, cond_vec_I2) ;
                    all_Rs_VIF(2, i_chan, i_subj, 2) = Rs ;
                    for ipred=1:size(pred_tmp,2)
                        corr_tmp = corr(pred_tmp(:,ipred), y_tmp_I2) ; % (1, nF*NT) matrix
                        all_Rs_IV(npred_first+ipred, i_chan,:,:, i_subj,2) = reshape(corr_tmp.^2, nF_data, nT_data) ;
                    end
                    end
                    % ============================= %
                    if CB_PT
                        for i_perm=1:n_perm_bloc
                            beta = regress_Y_on_X(y_tmp_I2, pred_tmp(randperm(n_ep_tmp),:), random_inter, cond_vec_I2) ;
                            beta_CP(n_betas-1:n_betas, i_chan, :, :, i_subj, 2, i_perm) = ...
                                reshape(beta, 2, nF_data, nT_data) ; 
                        end
                    end
                    % ============================= % 
                end
                if CB_PT
                   disp(['-- IO fitting with ', num2str(i_bloc*n_perm_bloc), ' permutations ', ...
                       'done for subject ', num2str(i_subj), '/', num2str(n_subj), ...
                       ' (',num2str(toc(ts)),' sec)']) 
                end
            end
            
            if CB_PT
                % =============================================================== %
                % ======== Compute shuffled cluster statistics for i_bloc              
                % =============================================================== %
                perm_cluster_t_i = compute_cluster_perm_stats(beta_CP, ...
                    cf_alpha, all_names, n_perm_bloc, n_betas, n_chan, nF_data, nT_data, n_subj) ; 
                % perm_cluster_t = zeros(n_perm, n_betas, 2, 2) ; 
                perm_cluster_t((1+n_perm_bloc*(i_bloc-1)):(n_perm_bloc*i_bloc),:,:,:) = ...
                    perm_cluster_t_i ; 
            end
        
        end
        
        % =============================================================== %
        % ======== Compute significance of beta coefficients (real)
        % =============================================================== %
        % all_beta: (n_betas, n_chan, nF_data, nT_data, n_subj, 2)
              
        %all_t_beta,all_p_beta: (n_betas, n_chan, nF_data, nT_data, 2) ; 
        for i_beta = 1:n_betas
            for i_I=1:2    
                curr_betas = shiftdim(all_beta(i_beta, :, :, :, :,i_I),1) ;
                % (n_chan, nF_data, nT_data, n_subj)
                curr_betas = permute(curr_betas, [4, 1, 2, 3]) ;
                curr_betas = reshape(curr_betas, n_subj, n_chan*nF_data*nT_data) ;
                
                [p_values, h_mask, ~, ~, ~, t_stats] = ...
                    compute_significance_matrix(curr_betas, cf_alpha, 10, ...
                    paired_tests, mean_H0, always_t_test, two_sided) ;%compute_significance_matrix
                
                t_stats = reshape(t_stats, n_chan, nF_data, nT_data) ; 
                p_values = reshape(p_values, n_chan, nF_data, nT_data) ; 
                all_t_beta(i_beta, :, :, :, i_I) = t_stats ; 
                all_p_beta(i_beta, :, :,:, i_I) = p_values ; 
            end
        end
        
        if CB_PT
            % =============================================================== %
            % ======== Compute shuffled cluster statistics
            % =============================================================== %
            criticals = compute_criticals(perm_cluster_t, n_betas) ; 
            % =============================================================== %
            % ======== Compute CB significance for real data
            % =============================================================== %
            [p_beta_CP, p_beta_CP_max] = compute_cluster_pv(all_t_beta, all_p_beta, ...
                perm_cluster_t, criticals, all_names, n_perm, ...
                n_betas, n_chan, nF_data, nT_data, cf_alpha) ; 
            % p_beta_CP, p_beta_CP_max: (n_betas, n_chan, nF_data, nT_data, 2) ; 
        end
        
        
        disp(['=== Time to fit all regressions (1 param): ', num2str(toc(ts))]) ; 
        if CB_PT
            save(fn_IO_fit, '-v7.3', 'indices_TP', 'time_vec', 'fs', 'xstep', ...
                'chanlocs', 'y_I1', 'y_I2', 'time_vec_IO', 'all_beta', 'all_R_squared', ...
                'all_MSE', 'all_t_beta', 'all_p_beta', 'time_vec_TF', 'f_begin', ...
                'f_step', 'f_stop', 'p_beta_CP', 'p_beta_CP_max', 'all_Rs_VIF', 'all_Rs_IV') ; 
        else
            save(fn_IO_fit, '-v7.3', 'indices_TP', 'time_vec', 'fs', 'xstep', ...
                'chanlocs', 'y_I1', 'y_I2', 'time_vec_IO', 'all_beta', 'all_R_squared', ...
                'all_MSE', 'all_t_beta', 'all_p_beta', 'time_vec_TF', 'f_begin', ...
                'f_step', 'f_stop', 'all_Rs_VIF', 'all_Rs_IV') ; 
        end
        % Save avg EEG with median splittings
        save(fn_avg_split, '-v7.3', 'labels_split', 'y_I1_split', 'y_I2_split', ...
            'time_vec_IO', 'f_begin', 'f_stop', 'f_step', 'N_I1_split', ...
            'N_I2_split', 'chanlocs') ; 
        if TF_regress
           save(fn_IO_fit,'all_TFs',  '-append') ;  
        end

    case 2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% All parameters (for model comparisons)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %all_IO_data: 
        % 3, n_stim_test, n_conds_max, n_subj, n_params
        % dim1: surprise - confidence - p1
        
        % Outcomes of the fitting (1 regression per param)
        % - 1 regression with surprise and residual conf as predictors --> 2 betas
        % /!\ separately for I_1 and I_2 --> x 2
        

        
        all_beta = NaN(2, n_chan, nF_data, nT_data, n_subj, 2, n_params) ; 
        all_R_squared = NaN(n_chan, nF_data, nT_data, n_subj, 2, n_params) ; 
        all_MSE = NaN(n_chan, nF_data, nT_data, n_subj, 2, n_params) ; 

        ts = tic ; 
        for i_subj=1:n_subj
            s_tmp = squeeze(all_s(:,:,i_subj)) ;    % (n_stim_test, n_conds_max)
            
            cond_vec = indices_TP(:,i_subj) ;       % (n_cond_max, 1)
            cond_vec = repmat(cond_vec', n_stim_test,1) ; 
            % ==*== Merging all conditions
            s_tmp = reshape(s_tmp,[1,N_obs]) ;            
            cond_vec = (reshape(cond_vec, [1, N_obs]))' ; 

            if reject_ep
                ep_kept_tmp = reshape(squeeze(all_ep_kept(:,:,i_subj)), ...
                    [1, N_obs]) ; 
                % all_ep_kept = NaN(n_stim_test, n_conds_max, n_subj) ; 
                s_tmp = s_tmp(ep_kept_tmp) ; 
                cond_vec = cond_vec(ep_kept_tmp) ; 
            end
            cond_vec_I1 = cond_vec(s_tmp==1) ;
            cond_vec_I2 = cond_vec(s_tmp==2) ;
            
            for i_chan=all_chan_idx
                y_tmp = downsample(squeeze(all_EEG_data(i_chan, :, :, :, i_subj)), down_factor) ; 
                y_tmp = y_tmp(idx_m:idx_M,:,:) ; 
                % (n_times, n_stim_test, n_conds_max)

                % ==*== Merging all conditions
                y_tmp = reshape(y_tmp,[n_times_IO,N_obs]) ;
                if reject_ep
                    y_tmp = y_tmp(:,ep_kept_tmp) ; 
                end
                y_tmp_I1 = (y_tmp(:,s_tmp==1))' ;
                y_tmp_I2 = (y_tmp(:,s_tmp==2))' ;
                
                if TF_regress
                    % y_tmp_I1: (n_epochs, n_times_IO)
                    [all_STFT, all_freq_Hz, all_t_STFT,~,f0,Fb] = ...
                        extract_time_freq_signal(y_tmp_I1, fs/down_factor, time_vec_IO, ...
                        'f_begin', f_begin, 'f_stop', f_stop, 'f_step', ...
                        f_step,'win_width_sec', win_width_sec, 'stand_STFT', ...
                        0, 'use_cwt', use_cwt, 'use_stft', use_stft, ...
                        'max_filt_order', 100, 'f_half_width', f_step/2, ...
                        'stft_with_conv', stft_with_conv) ; % [nF, nT, n_epoch]
                    all_STFT = downsample(permute(abs(all_STFT),[2, 1, 3]),downf_TF) ; % (nT, nF, n_ep)
                    all_STFT = permute(all_STFT, [3,2,1]) ; % (n_ep, nF, nT)
                    y_tmp_I1 = reshape(all_STFT,size(y_tmp_I1,1), nF_TF*nT_TF) ;
                    
                    [all_STFT, all_freq_Hz, all_t_STFT,~,f0,Fb] = ...
                        extract_time_freq_signal(y_tmp_I2, fs/down_factor, time_vec_IO, ...
                        'f_begin', f_begin, 'f_stop', f_stop, 'f_step', ...
                        f_step,'win_width_sec', win_width_sec, 'stand_STFT', ...
                        0, 'use_cwt', use_cwt, 'use_stft', use_stft, ...
                        'max_filt_order', 100, 'f_half_width', f_step/2, ...
                        'stft_with_conv', stft_with_conv) ; % [nF, nT, n_epoch]
                    all_STFT = downsample(permute(abs(all_STFT),[2, 1, 3]),downf_TF) ; % (nT, nF, n_ep)
                    all_STFT = permute(all_STFT, [3,2,1]) ; % (n_ep, nF, nT)
                    y_tmp_I2 = reshape(all_STFT,size(y_tmp_I2,1), nF_TF*nT_TF) ;
                end                
                
                for j_param=1:n_params
                    IO_data_tmp = squeeze(all_IO_data(:, :, :, i_subj, j_param)) ;
                    % dim1: surprise - confidence  p1
                    % x n_stim_test, n_conds_max, n_subj

                    % ==*== Merging all conditions
                    IO_data_tmp = (reshape(IO_data_tmp, [n_IO_out, N_obs]))' ;
                    if reject_ep
                        IO_data_tmp = IO_data_tmp(ep_kept_tmp,:) ; 
                    end
                    
                    % =========================================================== %
                    % ===*=== I_1
                    IO_data_I1 = IO_data_tmp(s_tmp==1,:) ; 
                    % columns: surprise - confidence - p1 - N_used
                    IO_p1_I1 = IO_data_I1(:,3) ; % attention: need to extract p1_mean BEFORE z-score, because log2 later
                    IO_N_I1 = IO_data_I1(:,4) ;
                    IO_data_I1 = zscore(IO_data_I1) ; 
                    
                    % -----> compute RESIDUAL confidence
                    if rconf_BMC
                        if rconf_no_logN
                            [~, ~, ~, eps_c] = regress_Y_on_X(IO_data_I1(:,2), ...
                                [IO_p1_I1, IO_p1_I1.^2, log2(IO_p1_I1), log2(1-IO_p1_I1), ...
                                log2(IO_N_I1)], random_inter, cond_vec_I1) ;
                        else
                            [~, ~, ~, eps_c] = regress_Y_on_X(IO_data_I1(:,2), ...
                                [IO_p1_I1, IO_p1_I1.^2, log2(IO_p1_I1), log2(1-IO_p1_I1)], random_inter, cond_vec_I1) ;
                        end

                        % -----> surprise + (RESIDUAL) confidence
                        [beta, mse, Rs] = regress_Y_on_X(y_tmp_I1, [IO_data_I1(:,1),zscore(eps_c)], random_inter, cond_vec_I1) ;
                    else
                        [beta, mse, Rs] = regress_Y_on_X(y_tmp_I1, IO_data_I1(:,1:2), random_inter, cond_vec_I1) ;
                    end

                    all_beta(:,i_chan, :, :, i_subj, 1, j_param) = reshape(beta, 2, nF_data, nT_data) ; %(2, n_chan, n_times_IO, n_subj, 2) ;
                    all_R_squared(i_chan, :, :, i_subj, 1, j_param) = reshape(Rs, nF_data, nT_data) ;  %(n_chan, n_times_IO, n_subj, 2) ;
                    all_MSE(i_chan, :, :, i_subj, 1, j_param) = reshape(mse, nF_data, nT_data) ; %(n_chan, n_times_IO, n_subj, 2) ;
                    
                    % =========================================================== %
                    % ===*=== I_2
                    IO_data_I2 = IO_data_tmp(s_tmp==2,:) ;
                    IO_p1_I2 = IO_data_I2(:,3) ;
                    IO_N_I2 = IO_data_I2(:,4) ;
                    IO_data_I2 = zscore(IO_data_I2) ; 

                    if rconf_BMC
                        % -----> compute RESIDUAL confidence
                        if rconf_no_logN
                            [~, ~, ~, eps_c] = regress_Y_on_X(IO_data_I2(:,2), ...
                                [IO_p1_I2, IO_p1_I2.^2, log2(IO_p1_I2), log2(1-IO_p1_I2), ...
                                log2(IO_N_I2)], random_inter, cond_vec_I2) ;
                        else
                            [~, ~, ~, eps_c] = regress_Y_on_X(IO_data_I2(:,2), ...
                               [IO_p1_I2, IO_p1_I2.^2, log2(IO_p1_I2), log2(1-IO_p1_I2)], random_inter, cond_vec_I2) ;
                        end

                        % -----> surprise + RESIDUAL confidence
                        [beta, mse, Rs] = regress_Y_on_X(y_tmp_I2, [IO_data_I2(:,1),zscore(eps_c)], ...
                            random_inter, cond_vec_I2) ;
                    else
                        [beta, mse, Rs] = regress_Y_on_X(y_tmp_I2, IO_data_I2(:,1:2), ...
                            random_inter, cond_vec_I2) ;
                    end

                    all_beta(:,i_chan, :, :, i_subj, 2, j_param) = reshape(beta, 2, nF_data, nT_data) ; %(2, n_chan, n_times_IO, n_subj, 2, n_params) ;
                    all_R_squared(i_chan, :, :, i_subj, 2, j_param) = reshape(Rs, nF_data, nT_data) ; %(n_chan, n_times_IO, n_subj, 2, n_params) ;
                    all_MSE(i_chan, :, :, i_subj, 2, j_param) = reshape(mse, nF_data, nT_data) ; %(n_chan, n_times_IO, n_subj, 2, n_params) ;
                end
            end
        end
        
        
        % Remove dimensions with single entries:
        %   if single_elec --> n_chan = 1
        %   if  TF_regress --> nF_data = 1
        all_beta = squeeze(all_beta) ; 
        all_R_squared = squeeze(all_R_squared) ;
        all_MSE = squeeze(all_MSE) ; 
        disp(['=== Time to fit all regressions (',num2str(n_params),...
            ' params): ', num2str(toc(ts))]) ; 
        xstep_IO = time_vec_IO(2) - time_vec_IO(1) ; 
        xstart_IO = time_vec_IO(1) ; 
        save(fn_IO_fit, '-v7.3', 'indices_TP', 'xstart', 'xstep', 'fs', ...
            'chanlocs', 'xstep_IO', 'xstart_IO','all_beta', 'all_R_squared', ...
            'all_MSE', 'all_params', 'all_ep_kept', 'time_vec_TF', 'f_begin', ...
            'f_step', 'f_stop') ; 
        % for y_I1 and y_I2 (eg to show GFP along with the IO fitting) -->
        % can load fn_avg_eeg!
end

if save_EEG_model && IO_fit_opt==1
    chan_names = {chanlocs.labels} ; 
    
    % all_s: NaN(n_stim_test, n_conds_max, n_subj) ;
    % all_ep_kept: true(n_stim_test, n_conds_max, n_subj) ; 
    % all_EEG_data: NaN(n_chan, n_times, n_stim_test, n_conds_max, n_subj) ; --> downsample(..., down_factor)
    
    % surprise - confidence - update - unpredictability - p1 -
        % prediction error - N_used (= N observations for current
        % inference)
    % all_IO_data:  NaN(7, n_stim_test, n_conds_max, n_subj) ; % PE - confidence - update - unpredictability - p1 - surprise - N_used
    
    
    % *** A *** EEG data
    % ========== I1
    all_EEG_I1_model = all_EEG_data(:,1:down_factor:end,:,:,:) ; 
    all_EEG_I1_model = all_EEG_I1_model(:,idx_m:idx_M,:,:,:) ; 
    all_EEG_I1_model(:,:, all_s==2 & ~all_ep_kept) = NaN ; 
    
    fn_EEG_I1 = [fn_data_tmp, 'EEG_I1', str_params, str_ampl_max, '.mat'] ; 
    save(fn_EEG_I1, '-v7.3', 'all_EEG_I1_model', 'chan_names', 'indices_TP', 'time_vec_IO') ; 
    clear all_EEG_I1_model ;     
    
    % ========== I2
    all_EEG_I2_model = all_EEG_data(:,1:down_factor:end,:,:,:) ; 
    all_EEG_I2_model = all_EEG_I2_model(:,idx_m:idx_M,:,:,:) ; 
    all_EEG_I2_model(:,:, all_s==1 & ~all_ep_kept) = NaN ; 
    
    fn_EEG_I2 = [fn_data_tmp, 'EEG_I2', str_params, str_ampl_max, '.mat'] ; 
    save(fn_EEG_I2, '-v7.3', 'all_EEG_I2_model', 'chan_names', 'indices_TP', 'time_vec_IO') ; 
    clear all_EEG_I2_model ; 
    
    % *** B *** regressors
    % ========== I1
    PE_conf_I1 = all_IO_data(1:2,:,:,:) ;
    PE_conf_I1(:,all_s==2 & ~all_ep_kept) = NaN ; 
    fn_reg_I1 = [fn_data_tmp, 'PE_conf_I1', str_params, str_ampl_max, str_model, '.mat'] ; 
    save(fn_reg_I1, '-v7.3', 'PE_conf_I1', 'indices_TP') ; 
    
    
    % ========== I2
    PE_conf_I2 = all_IO_data(1:2,:,:,:) ;
    PE_conf_I2(:,all_s==1 & ~all_ep_kept) = NaN ; 
    fn_reg_I2 = [fn_data_tmp, 'PE_conf_I2', str_params, str_ampl_max, str_model, '.mat'] ; 
    save(fn_reg_I2, '-v7.3', 'PE_conf_I1', 'indices_TP') ; 
    
    
    disp('************************ EEG and regressors SAVED')
end



end

