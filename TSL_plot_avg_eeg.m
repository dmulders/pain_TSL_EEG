function [chanlocs,chan_desired, n_chan] = TSL_plot_avg_eeg()
% Reload the data allowing to display the average EEG responses
%   - merging all conditions: y_I1 & y_I2
%   - per condition: y_I1_cond & y_I2_cond
%   - per parameter (TPs): y_I1_TP & y_I2_TP

% required data are saved in TSL_analyze_EEG.m w/ save_avg_eeg = 1 ; 

% ======================================================================= %
% ===== * ===== General parameters
% ======================================================================= %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
all_subj_idx    = setdiff(1:36,[1, 11, 15, 28, 33]) ; 
LP_filter       = 1 ;       % Low-pass filter the epochs
LP_freq         = 30 ;      % frequency cutoff for the LP filtering
reject_ampl     = 1 ;       % reject epochs based on amplitude
max_ampl_rej    = 80 ;      % amplitude criterion to reject epochs
                            % applied after LP AND after bc...
                            % EEG data of all subjects. 
show_merged_avg = 1 ;       % only show the avg merging all conditions
moving_avg_ep   = 0 ;       % reload data where moving avg was applied on each epoch
moving_avg      = 0 ;       % moving avg on EEG from each subject
movmu_ms        = 80 ;      % msec for the moving avg
reref           = 0 ; 
elecs_ref = {'T7','T8'} ; %'PZ' ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Time-frequency params
plot_TF_maps = 0 ;          % 
use_stft = 1 ;              % TF parameter; otherwise: CWTs
chan_TF = {'FCZ'} ;           % FCZ CZ CPZ C3
%chan_TF = {'FCZ', 'CZ', 'CPZ', 'PZ'} ; 
%chan_TF = {'FCZ', 'CZ', 'CPZ', 'C3', 'C4'} ; 
win_width_sec = 0.1 ; %0.1 ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
show_topos      = 0 ; 
% Topoplot options
% times at which to show topo, in ms
t_topo_I1 = [200, 330, 470] ;
% N1-N2: 190-200ms (mixed, N2 but lateralized)
% P2: 330
t_topo_I2 = [355, 508, 682] ; 
% N1-N2: 355 (cfr TSL_plot_IO_fit: 208ms is not lateralized here)
% P2: 508

% TF
str_avg_TF = '' ; 
if use_stft
   str_avg_TF = [str_avg_TF, '_w', num2str(round(100*win_width_sec))] ; % window width
end

% to find latencies of peaks, based on grand mean plots
N2_I1_win = [170, 236] ; P2_I1_win = [240, 410] ; 
N2_I2_win = [316, 430] ; P2_I2_win = [450, 600 ] ; %668] ; 

taille_axis = 24 ; 
taille_axis_cbar = 15 ; 
avg_cmap_levels = 0 ; 
numC_topo = 8 ;                 % if avg topoplots use filled levels 
topo_style = 'both' ;           % 'both', 'map', cfr test_scalp_plot
n_colors_def = 64 ;             % default nb of colors for the colormaps
if avg_cmap_levels
    cmap_topo = parula(numC_topo) ; %parula ; %jet ; %gnuplot2 ; %jet ; %plasma ; % viridis
else
    cmap_topo = parula(n_colors_def) ; %jet(n_colors_def) ; %parula ; 
    cmap_topo = jet(n_colors_def) ;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% General plot params --------------------------------------------------- %
eps_fig = 0 ; fig_fig = 0 ; pdf_fig = 0 ; % otherwise: .png
% Specific plot params -------------------------------------------------- %
plot_ms = 1 ;               % msec time unit

% Define channels to show on plot
only_plot_topo_chan = 0 ;   % don't show channels out of topoplot in butterfly plot
topo_colors = 1 ;           % show color codes for the avg signals on topo
selected_chan = 1 ;         % only plot chan_selected on the butterfly plot
chan_selected = {'FCZ', 'C3', 'C4', 'CPz', 'Cz'} ;
% ----------------------------------------------------------------------- %


fn_res = ['./results_EEG/'] ; 

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
str_movmu_ep = '' ; 
if moving_avg_ep
    str_movmu_ep = ['_mmu', num2str(movmu_ms)] ;
end
str_ref = '' ; 
if reref 
    for ie=1:length(elecs_ref)
        str_ref = [str_ref,'_', elecs_ref{ie}] ;
    end
end

str_params = [str_subj, str_LP, str_movmu_ep] ;

n_chans_plot = length(chan_selected) ;
str_plot = '' ; 
if selected_chan
    str_plot = ['_',num2str(n_chans_plot),'chans'] ; 
elseif only_plot_topo_chan
    str_plot = '_topo_chans' ; 
end
str_plot = [str_plot, str_ref] ; 

fn_avg_eeg = [fn_res, 'data/Avg_eeg', str_params, str_ampl_max, '.mat'] 

saved_data = load(fn_avg_eeg) ;
indices_TP = saved_data.indices_TP ; 
time_vec = saved_data.time_vec;
fs = saved_data.fs;
xstep = saved_data.xstep;
chanlocs = saved_data.chanlocs;
all_TPs = saved_data.all_TPs;
all_TPs_str = saved_data.all_TPs_str;
y_I1 = saved_data.y_I1; % (n_chan, n_times, n_subj) 
y_I2 = saved_data.y_I2 ; 
y_I1_TP = saved_data.y_I1_TP; % (n_chan, n_times, n_subj, n_TPs) ; 
y_I2_TP = saved_data.y_I2_TP ; 
y_I1_cond = saved_data.y_I1_cond; % (n_chan, n_times, n_subj, n_conds_max)
y_I2_cond = saved_data.y_I2_cond ; 

[n_chan, n_times, ~] = size(y_I1) ; 
n_TPs = length(all_TPs_str) ; 

chan_topo = logical([chanlocs.topo_enabled]) ;
chanlocs_butterfly = chanlocs ; 
% only consider chan_selected for the plots
[~, chan_desired] = define_chan_order_colors(chanlocs, chan_selected) ; 
if selected_chan
    chanlocs_butterfly = chanlocs(chan_desired) ; 
elseif only_plot_topo_chan && topo_colors
    chanlocs_butterfly = chanlocs(chan_topo) ; 
end


nbins_mmu = 2*floor(fs*movmu_ms/1000/2)+1 ;

% in sec
xlim_butterfly = [-0.25, 1] ; %[-0.25,0.75] ; % [time_vec(1), time_vec(end)] ; 
time_vec_used = time_vec ; 
time_ms = 1000.*time_vec_used ; 
xlab = 'time after stimulus (s)' ; 
if plot_ms
    % in seconds
    xlim_butterfly = 1000.*xlim_butterfly ;
    time_vec_used = time_ms ; 
    xlab = 'time after stimulus (ms)' ; 
end

if selected_chan
    n_curves = n_chans_plot ; 
elseif only_plot_topo_chan
    n_curves = length(chan_topo) ; 
end
if n_curves==1
    cmap_but = [0, 128, 255]./255 ; 
else
    cmap_but = viridis(n_curves+1) ; %jet(n_curves) ;   % 
    % perceptually uniform & colorblind friendly cmaps: 
    % viridis, inferno, plasma, magma
    cmap_but = cmap_but(2:end,:) ; 
end


% ------------------------------------------------------------------- %
% --*-- Find latencies of N2 and P2
% ------------------------------------------------------------------- %
[~,istart] = min(abs(time_ms-N2_I1_win(1))) ; [~,istop] = min(abs(time_ms-N2_I1_win(2))) ;  
idx_N2_I1 = istart:istop ;
[~,istart] = min(abs(time_ms-P2_I1_win(1))) ; [~,istop] = min(abs(time_ms-P2_I1_win(2))) ;  
idx_P2_I1 = istart:istop ;
[~,istart] = min(abs(time_ms-N2_I2_win(1))) ; [~,istop] = min(abs(time_ms-N2_I2_win(2))) ;  
idx_N2_I2 = istart:istop ;
[~,istart] = min(abs(time_ms-P2_I2_win(1))) ; [~,istop] = min(abs(time_ms-P2_I2_win(2))) ;  
idx_P2_I2 = istart:istop ; 

all_lat = NaN(4, n_subj) ; 
chan_lat = 'FCZ' ; 
[~, idx_chan_lat] = define_chan_order_colors(chanlocs, {chan_lat}) ; 
fn_lat = [fn_res, 'lat', str_subj, '_', chan_lat, '.txt'] ;

% y_I1, y_I2: (n_chan, n_times, n_subj) 

%%%%%%%%%%%%%%%%%% I1 
y_tmp = squeeze(y_I1(idx_chan_lat,:,:)) ; % n_times, n_subj
% ==== N2
[N2_I1,i_N2_I1] = min(y_tmp(idx_N2_I1,:),[],1) ; %(1, n_subj)
all_lat(1, :) = time_ms(idx_N2_I1(1) + i_N2_I1-1) ; 
% ==== P2
[P2_I1, i_P2_I1] = max(y_tmp(idx_P2_I1,:),[],1) ; 
all_lat(2, :) = time_ms(idx_P2_I1(1) + i_P2_I1-1) ;     

%%%%%%%%%%%%%%%%%% I2
y_tmp = squeeze(y_I2(idx_chan_lat, :,:)) ; 
% ==== N2
[N2_I2, i_N2_I2] = min(y_tmp(idx_N2_I2,:),[],1) ;
all_lat(3, :) = time_ms(idx_N2_I2(1) + i_N2_I2-1) ; 
% ==== P2
[P2_I2, i_P2_I2] = max(y_tmp(idx_P2_I2,:),[],1) ; 
all_lat(4, :) = time_ms(idx_P2_I2(1) + i_P2_I2-1) ; 

str_lat_I1 = ['[I1] Latencies of N2: ', num2str(round(mean(all_lat(1,:)),4)),'(', num2str(round(std(all_lat(1,:)),4)),')', ...
    '; P2: ', num2str(mean(all_lat(2,:))),'(', num2str(std(all_lat(2,:))),')'] ; 
str_lat_I2 = ['[I2] Latencies of N2: ', num2str(round(mean(all_lat(3,:)),4)),'(', num2str(round(std(all_lat(3,:)),4)),')', ...
    '; P2: ', num2str(mean(all_lat(4,:))),'(', num2str(std(all_lat(4,:))),')'] ; 
disp(str_lat_I1) ; disp(str_lat_I2) ; 
write_empty_lines(fn_lat, 0,'w') ;
write_to_txt(fn_lat, str_lat_I1,'a',2) ;
write_to_txt(fn_lat, str_lat_I2,'a',2) ;


% ------------------------------------------------------------------- %
% --*-- Merging all TPS (idem in TSL_analyze_EEG.m)
% ------------------------------------------------------------------- %
if show_merged_avg
    if reref
       [~, idx_elec] = define_chan_order_colors(chanlocs, elecs_ref) ; 
       %y_I1: (n_chan, n_times, n_subj) 
       y_I1 = y_I1 - repmat(mean(y_I1(idx_elec,:,:),1), n_chan, 1, 1) ;
       y_I2 = y_I2 - repmat(mean(y_I2(idx_elec,:,:), 1), n_chan, 1, 1) ;
    end
    
    if moving_avg
        disp('Doing movmean!')
        y_I1_avg = squeeze(mean(movmean(y_I1,nbins_mmu,2),3)) ; % n_chan, n_times
        y_I2_avg = squeeze(mean(movmean(y_I2,nbins_mmu,2),3)) ; 
    else
        y_I1_avg = squeeze(mean(y_I1,3)) ; % n_chan, n_times
        y_I2_avg = squeeze(mean(y_I2,3)) ; 
    end
    
    if selected_chan
        y_I1_avg = y_I1_avg(chan_desired, :) ; 
        y_I2_avg = y_I2_avg(chan_desired, :) ; 
    elseif only_plot_topo_chan
        y_I1_avg = y_I1_avg(chan_topo, :) ; 
        y_I2_avg = y_I2_avg(chan_topo, :) ; 
    end
    
    % non overlapping means:
    %  mean(reshape(data, [], 5),2);
    %  movmean(A,k,dim)
    
    % ------------------------------------------------------------------- %
    % --*-- Compute GFP
    % ------------------------------------------------------------------- %
    %y_I1, y_I2: (n_chan, n_times, n_subj) 
    GFP_I1 = squeeze(compute_GFP(y_I1)) ; % (n_times, n_subj)
    GFP_I2 = squeeze(compute_GFP(y_I2)) ; 
    mean_GFP_I1 = (squeeze(mean(GFP_I1,2)))' ;  mean_GFP_I2 = (squeeze(mean(GFP_I2,2)))' ; 
    std_GFP_I1 = (squeeze(std(GFP_I1,[],2)))' ; std_GFP_I2 = (squeeze(std(GFP_I2,[],2)))' ; 

    % ------------------------------------------------------------------- %
    % --*-- I_1
    % ------------------------------------------------------------------- %
    if selected_chan
        fn_res_I1 = [fn_res, num2str(n_chans_plot), 'chans_I1/'] ; 
        if ~exist(fn_res_I1, 'dir') ;  mkdir(fn_res_I1) ; end
    else
        fn_res_I1 = fn_res ; 
    end
    gray_col = [200,200,200]./255 ; 
    %%%%%%%%%%%%%%% ADD GFP
    fn_I1 = [fn_res_I1, 'Butterfly_I1', str_params,str_movmu, str_plot] ; 
    fig_h = figure('units','centimeters','outerposition',[1 2 18 12] ,...
        'Name', 'I_1') ; 
    hold on ; 
    % shaded std
    fill([time_vec_used, fliplr(time_vec_used)], ...
        [mean_GFP_I1-std_GFP_I1, fliplr(mean_GFP_I1+std_GFP_I1)], gray_col, ...
        'EdgeColor','none', 'FaceAlpha', 0.2) ; hold on ;      
    plot(time_vec_used,mean_GFP_I1,'color', gray_col, 'linewidth', 2) ;     
    hold on ; 
    
    plot_butterfly(time_vec_used, y_I1_avg, 'fn_save', fn_I1, 'xlim_val', ...
        xlim_butterfly, 'x_lab', xlab, 'topo_colors', ...
        topo_colors, 'chanlocs', chanlocs_butterfly, 'eps_fig', eps_fig, ...
        'fig_fig', fig_fig, 'pdf_fig', pdf_fig, 'create_fig', 0, 'save_fig', 0, ...
        'cmap', cmap_but) ; 
    hold on ; 
    
    % Indicate times at which show topo
    yL = get(gca, 'ylim') ; 
    for t_I1=t_topo_I1 
        [~,idx_topo] = min(abs(time_ms-t_I1));  
        plot([t_I1, t_I1],[yL(1),mean_GFP_I1(idx_topo)],':','color', 'k', 'linewidth', 1) ; 
        hold on ; 
    end
    
    if topo_colors ; fn_I1 = [fn_I1, '_topo'] ; end
    my_save_fig(fig_h, fn_I1, eps_fig, fig_fig, pdf_fig) ;

    % ------------------------------------------------------------------- %
    % --*-- I_2
    % ------------------------------------------------------------------- %
    if selected_chan
        fn_res_I2 = [fn_res, num2str(n_chans_plot), 'chans_I2/'] ; 
        if ~exist(fn_res_I2, 'dir') ;  mkdir(fn_res_I2) ; end
    else
        fn_res_I2 = fn_res ; 
    end
    
    %%%%%%%%%%%%%%% ADD GFP
    fig_h = figure('units','centimeters','outerposition',[1 2 18 12] ,...
        'Name', 'I_2') ; 
    hold on ; 
    % shaded std
    fill([time_vec_used, fliplr(time_vec_used)], ...
        [mean_GFP_I2-std_GFP_I2, fliplr(mean_GFP_I2+std_GFP_I2)], gray_col, ...
        'EdgeColor','none', 'FaceAlpha', 0.2) ; hold on ;      
    plot(time_vec_used,mean_GFP_I2,'color', gray_col, 'linewidth', 2) ;     
    hold on ; 
    %%%%%%%%%%%%%%% other chans
    fn_I2 = [fn_res_I2, 'Butterfly_I2', str_params,str_movmu,str_plot] ; 
    plot_butterfly(time_vec_used, y_I2_avg, 'fn_save', fn_I2, 'xlim_val', ...
        xlim_butterfly, 'x_lab', xlab, 'topo_colors', ...
        topo_colors, 'chanlocs', chanlocs_butterfly, 'eps_fig', eps_fig, ...
        'fig_fig', fig_fig, 'pdf_fig', pdf_fig, 'create_fig', 0, 'save_fig', 0, ...
        'cmap', cmap_but) ; 
    % Indicate times at which show topo
    yL = get(gca, 'ylim') ; 
    for t_I2=t_topo_I2 
        [~,idx_topo] = min(abs(time_ms-t_I2));  
        plot([t_I2, t_I2],[yL(1),mean_GFP_I2(idx_topo)],':','color', 'k', 'linewidth', 1) ; 
        hold on ; 
    end
    
    if topo_colors ; fn_I2 = [fn_I2, '_topo'] ; end
    my_save_fig(fig_h, fn_I2, eps_fig, fig_fig, pdf_fig) ;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Topoplots at selected times
    if show_topos
        for t_I1=t_topo_I1

            % Extract topo at time t_I1
            [~,idx_topo] = min(abs(time_ms-t_I1));  

            % 1. standardize topo of each subject
            y_t = squeeze(y_I1(:,idx_topo,:)) ; % (n_chan, n_subj) 
            y_t = (y_t - repmat(mean(y_t,1),n_chan,1))./...
                (repmat(sqrt(sum(y_t.^2,1)), n_chan, 1)) ; 
            % 2. average over subjects
            y_t = mean(y_t,2) ; 

            fn_tmp = [fn_res_I1, 'Topo_I1_', num2str(t_I1), 'ms', str_params,str_movmu] ;
            plot_single_topo(y_t, chanlocs, fn_tmp, ['I_1 - t = ', num2str(t_I1), ' ms'], ...
                avg_cmap_levels, eps_fig, fig_fig, pdf_fig, topo_style, numC_topo, cmap_topo) ; 
        end
        % Idem I2
        for t_I2=t_topo_I2

            % Extract topo at time t_I1
            [~,idx_topo] = min(abs(time_ms-t_I2));  

            % 1. standardize topo of each subject
            y_t = squeeze(y_I2(:,idx_topo,:)) ; % (n_chan, n_subj) 
            y_t = (y_t - repmat(mean(y_t,1),n_chan,1))./...
                (repmat(sqrt(sum(y_t.^2,1)), n_chan, 1)) ; 
            % 2. average over subjects
            y_t = mean(y_t,2) ; 

            fn_tmp = [fn_res_I2, 'Topo_I2_', num2str(t_I2), 'ms', str_params,str_movmu] ;
            plot_single_topo(y_t, chanlocs, fn_tmp, ['I_2 - t = ', num2str(t_I2), ' ms'], ...
                avg_cmap_levels, eps_fig, fig_fig, pdf_fig, topo_style, numC_topo, cmap_topo) ; 
        end
    end
    
    
    % ------------------------------------------------------------------- %
    % --*-- Time-frequency maps (avg of TFs of epochs)
    % ------------------------------------------------------------------- %
    str_chans = '' ; %'' ; %'_FCZ' ; %'_FCZ' ; %'' ; 
    fn_res_TF = [fn_res, 'avg_TF/'] ; if ~exist(fn_res_TF, 'dir') ;  mkdir(fn_res_TF) ; end
    fn_avg_TF = [fn_res, 'data/TF', str_subj, '_maxA150', str_avg_TF,str_chans,'.mat'] 
    if plot_TF_maps && exist(fn_avg_TF, 'file')
        data_TF = load(fn_avg_TF) ; 
        y_TF_I1 = data_TF.y_TF_I1 ; % (n_chan, nF_TF, nT_TF, n_subj)
        y_TF_I2 = data_TF.y_TF_I2 ; 

        time_TF = data_TF.time_TF ; 
        freqs_TF = data_TF.freqs_TF ;  
        f_begin = data_TF.f_begin ; 
        f_step = data_TF.f_step ; 
        f_stop = data_TF.f_stop ; 
        ylab = 'frequency (Hz) ' ;

        time_TF = 1000.*time_TF ; 
        gapx = diff(time_TF(1:2)) ; 
        xlim_v = [time_TF(2),time_TF(end-1)+gapx/2] ;%[time_TF(2)-gapx/2,time_TF(end-1)+gapx/2] % [time_TF(1), time_TF(end)] ;

        xticks = 1000.*[-0.5:0.5:1] ; 
        yticks = 10:20:f_stop ; 
        cmap = jet ; %redblue ; jet ; % parula

        % ----- select data 
        cond1_TF = 0 ; str_norm = '' ; 
        if cond1_TF
            str_norm = '_C1' ; 
            y_TF_I1_c1 = data_TF.y_TF_I1_c1 ; % NEW (cond 1 only)
            y_TF_I2_c1 = data_TF.y_TF_I2_c1 ; % NEW (cond 1 only)
            y_TF_I1 = y_TF_I1_c1 ;
            y_TF_I2 = y_TF_I2_c1 ; 
        end

        % ----- Normalization, per participant
        bcorr_STFT = 1 ; zscore_STFT = 0 ; 
        if bcorr_STFT
            % pre-stimulus baseline subtraction
            % ref interval: [-0.5,-0.1]sec (idem Liberati2017)
            nT_TF = length(time_TF) ; 
            [~, idx_first] = min(abs(time_TF+500)) ;
            [~, idx_last] = min(abs(time_TF+100)) ;
            ref = repmat(mean(y_TF_I1(:, :, idx_first:idx_last,:),3), 1,1, nT_TF, 1) ;
            %y_TF_I1 = y_TF_I1 - ref ; 
            y_TF_I1 = (y_TF_I1 - ref)./ref.*100 ; % test w & w/o divide (dividing: idem Liberati2017)
            ref = repmat(mean(y_TF_I2(:, :,idx_first:idx_last,:),3), 1, 1,nT_TF, 1) ;
            %y_TF_I2 = y_TF_I2 - ref ; 
            y_TF_I2 = (y_TF_I2 - ref)./ref.*100 ; 
            str_norm = [str_norm, '_bc'] ; 
        elseif zscore_STFT
            y_TF_I1 = (y_TF_I1 - repmat(mean(y_TF_I1,3),1, 1,nT_TF, 1))./( repmat(std(y_TF_I1,0,3),1, 1,nT_TF, 1) ) ;
            y_TF_I2 = (y_TF_I2 - repmat(mean(y_TF_I2,3),1, 1,nT_TF, 1))./( repmat(std(y_TF_I2,0,3),1, 1,nT_TF, 1) ) ;
            str_norm = [str_norm, '_zs'] ; 
        end

        % ----- Select chan
        y_TF_I1_all = y_TF_I1 ; % keep all chans for topoplots
        y_TF_I2_all = y_TF_I2 ; % (n_chan, nF_TF, nT_TF, n_subj)
        [~, idx_chan_TF] = define_chan_order_colors(chanlocs, chan_TF) ; 
        str_chan_TF = '' ; 
        for icTF=1:length(chan_TF)
            str_chan_TF = [str_chan_TF, '_', chan_TF{icTF}] ; 
        end

        y_TF_I1 = squeeze(mean(y_TF_I1(idx_chan_TF, :, :, :),1)) ; % (nF_TF, nT_TF, n_subj)
        y_TF_I2 = squeeze(mean(y_TF_I2(idx_chan_TF, :, :, :), 1)) ; 


        % compare STFT amplitudes along different intervals
        all_freqs_inter = [60, 90; 60, 90] ;
        all_time_inter_I1 = [150, 1500; 170, 240] ;
        all_time_inter_I2 = [150, 1500; 320,430] ;
        % N2_I1_win = [170, 236] ; 
        % N2_I2_win = [316, 430] ;  
        n_inter = size(all_freqs_inter,1) ; 
        fn_gamma = [fn_res_TF, 'gamma_ER', str_avg_TF, str_subj, str_chan_TF, str_norm, '_',num2str(n_inter), 'TF_int.txt'] ; 
        write_to_txt(fn_gamma, '','w',1) ;
        show_topos_TF = 1 ; 

        for iinter = 1:n_inter
            freqs_inter = all_freqs_inter(iinter,:) ;
            time_inter_I1 = all_time_inter_I1(iinter,:) ;      
            time_inter_I2 = all_time_inter_I2(iinter,:) ;    
            % ----- extract amplitude of stimulus-evoked changes
            indf_gamma = msec_to_ind(freqs_TF, freqs_inter) ;      % frequency  
            indt_gamma_I1 = msec_to_ind(time_TF, time_inter_I1) ; 
            indt_gamma_I2 = msec_to_ind(time_TF, time_inter_I2) ; 
            gamma_I1 = squeeze(mean(mean(y_TF_I1(indf_gamma(1):indf_gamma(2), indt_gamma_I1(1):indt_gamma_I1(2),:),1),2)) ; % (n_subj,1)
            gamma_I2 = squeeze(mean(mean(y_TF_I2(indf_gamma(1):indf_gamma(2), indt_gamma_I2(1):indt_gamma_I2(2),:),1),2)) ; 

            gamma_ER = [gamma_I1, gamma_I2] ; 
            [p_values, h_mask, ~, ~, conf_intervals, test_statistic, cohend] = ...
                compute_significance_matrix(gamma_ER, 0.05, 10, 1, 0, 0, 1) ; 
            str_tmp = ['Diff. of evoked ampl. for I2-I1 in [',num2str(freqs_inter(1)),',',num2str(freqs_inter(2)),'] Hz and [',num2str(time_inter_I1(1)),',',num2str(time_inter_I1(2)),'] ms (I1) -- [',num2str(time_inter_I2(1)),',',num2str(time_inter_I2(2)),'] ms (I2): ', ...
                num2str(round(diff(nanmean(gamma_ER)),3)), ' (',...
                num2str(round(test_statistic(2,1),3)), ', p = ',...
                num2str(round(p_values(2,1),5)),', Cohen''s d = ',num2str(round(cohend(2,1),5)),')'] ; 
            disp(str_tmp) ; write_to_txt(fn_gamma, str_tmp,'a',1) ;

            if show_topos_TF && strcmpi(str_chans, '')
                % ----- plot topographies of TF amplitudes
                % 1. average along selected time-frequency interval
                % (n_chan, n_subj)
                topo_TF_I1 = squeeze(mean(mean(y_TF_I1_all(:,indf_gamma(1):indf_gamma(2), indt_gamma_I1(1):indt_gamma_I1(2),:),2),3)) ; 
                topo_TF_I2 = squeeze(mean(mean(y_TF_I2_all(:,indf_gamma(1):indf_gamma(2), indt_gamma_I2(1):indt_gamma_I2(2),:),2),3)) ; 

                % 2. standardize topo of each subject
                topo_TF_I1 = (topo_TF_I1 - repmat(mean(topo_TF_I1,1),n_chan,1))./...
                    (repmat(sqrt(sum(topo_TF_I1.^2,1)), n_chan, 1)) ; 
                topo_TF_I2 = (topo_TF_I2 - repmat(mean(topo_TF_I2,1),n_chan,1))./...
                    (repmat(sqrt(sum(topo_TF_I2.^2,1)), n_chan, 1)) ; 

                % 3. average over subjects
                topo_TF_I1 = mean(topo_TF_I1,2) ; topo_TF_I2 = mean(topo_TF_I2,2) ; 
                str_inter_I1 = ['_', num2str(freqs_inter(1)),num2str(freqs_inter(2)),'Hz_',...
                    num2str(time_inter_I1(1)),num2str(time_inter_I1(2)),'ms'] ; 
                str_inter_I2 = ['_', num2str(freqs_inter(1)),num2str(freqs_inter(2)),'Hz_',...
                    num2str(time_inter_I2(1)),'',num2str(time_inter_I2(2)),'ms'] ; 

                fn_tmp = [fn_res_TF, 'Topo_TF_I1', str_avg_TF, str_subj, str_norm, str_inter_I1] ;
                plot_single_topo(topo_TF_I1, chanlocs, fn_tmp, ['I_1 - ', str_inter_I1], ...
                    avg_cmap_levels, eps_fig, fig_fig, pdf_fig, topo_style, numC_topo, cmap_topo, 1, 'a.u.') ; 
                fn_tmp = [fn_res_TF, 'Topo_TF_I2', str_avg_TF, str_subj, str_norm, str_inter_I2] ;
                plot_single_topo(topo_TF_I2, chanlocs, fn_tmp, ['I_2 - ', str_inter_I2], ...
                    avg_cmap_levels, eps_fig, fig_fig, pdf_fig, topo_style, numC_topo, cmap_topo, 1, 'a.u.') ;    
            end
        end

        fig_sz = [1 2 32 12] ; taille = 20 ;     
        fmax_plot = 50 ; % limit y axis
        ylim_map = [f_begin-f_step/2,fmax_plot] ; 
        ylim_map = [40,f_stop] ; 
        %ylim_map = [15,f_stop] ; 
        ind_ylim = msec_to_ind(freqs_TF, ylim_map) ; 
        ind_xlim = msec_to_ind(time_TF, xlim_v) ;

        avg_subj = 1 ; % avg data of all participants (otherwise: show indiv maps)

        if avg_subj
            % ----- Avg data of all participants 
            y_TF_I1 = mean(y_TF_I1, 3) ; 
            y_TF_I2 = mean(y_TF_I2, 3) ;

            idx_map = 1 ; nmaps = 1 ; 
        else
            % ----- plot individual maps
            idx_map = 1:n_subj ; nmaps = n_subj ; 
        end

        fn_TF = [fn_res_TF, 'TF', str_avg_TF, str_subj, str_chan_TF, str_norm] ;

        for imap=1:nmaps
            fn_TF_tmp = fn_TF ; 
            y_tmp_I1 = squeeze(y_TF_I1(:,:,imap)) ; 
            y_tmp_I2 = squeeze(y_TF_I2(:,:,imap)) ; 

            if ~avg_subj
                fn_TF_tmp = [fn_TF, '_subj', num2str(all_subj_idx(imap))] ; 
            end
            fig_tmp = figure('units','centimeters','outerposition',fig_sz,...
                'Name', 'TF') ;     
            % --*-- ----> I_1    
            h_tmp(1) = subplot(1, 2, 1) ; hold on ;         
            I1_tmp = y_tmp_I1(ind_ylim(1):ind_ylim(2),ind_xlim(1):ind_xlim(2)) ; 
            I2_tmp = y_tmp_I2(ind_ylim(1):ind_ylim(2),ind_xlim(1):ind_xlim(2)) ; 

            ampl = max(max(I1_tmp(:)), max(I2_tmp(:))) ; 
            mampl = min(min(I1_tmp(:)), min(I2_tmp(:))) ; 
            curr_min = mampl ; curr_max = ampl ;  
            plot_heatmap(' ', time_TF, freqs_TF, y_tmp_I1, cmap, ...
                'ampl. (%)', ylab, xlab, 'map_opt', 1,  ...
                'eps_fig', eps_fig, 'fig_fig', fig_fig, 'pdf_fig', pdf_fig, ...
                'xticks', xticks, 'yticks', yticks, 'xtick_labs', xticks, ...
                'ytick_labs', yticks, 'clim_val',[curr_min, curr_max], 'create_fig', 0, ...
                'vline', 1, 'ylim_v', ylim_map, 'xlim_v', xlim_v) ; 
            title('I_1', 'FontSize', taille, 'FontName', 'Arial') ;     
            % --*-- ----> I_2
            h_tmp(2) = subplot(1, 2, 2) ; hold on ; 
            %curr_min = min(B_I2(:)) ; curr_max = max(B_I2(:)) ; 
            plot_heatmap('', time_TF, freqs_TF, y_tmp_I2, cmap, ...
                'ampl. (%)', ylab, xlab, 'map_opt', 1,  ...
                'eps_fig', eps_fig, 'fig_fig', fig_fig, 'pdf_fig', pdf_fig, ...
                'xticks', xticks, 'yticks', yticks, 'xtick_labs', xticks, ...
                'ytick_labs', yticks, 'clim_val',[curr_min, curr_max], 'create_fig', 0, ...
                'vline', 1, 'ylim_v', ylim_map, 'xlim_v', xlim_v) ; 
            title('I_2', 'FontSize', taille, 'FontName', 'Arial') ;
            %linkaxes(h_tmp,'y') ; 
            my_save_fig(fig_tmp, fn_TF_tmp, eps_fig, fig_fig, pdf_fig) ; 
        end
    end
    
    return
end






fn_res_tp = [fn_res, 'per_TP/'] ; 
if ~exist(fn_res_tp, 'dir') ;  mkdir(fn_res_tp) ; end
% ------------------------------------------------------------------- %
% --*-- Compare 2 TPs
% ------------------------------------------------------------------- %
i_TPs_compare = [2, 3] ; 
% all_TPs = [0.5, 0.5; 0.3, 0.7; 0.7, 0.3; 0.3, 0.3; 0.7, 0.7] ; 
str_TPs = [all_TPs_str{i_TPs_compare(1)}, all_TPs_str{i_TPs_compare(2)}] ; 
str_TPs = str_TPs(2:end) ; % remove 1st _

fn_comp_tp = [fn_res_tp, str_TPs, '/'] ; 
if ~exist(fn_comp_tp, 'dir') ;  mkdir(fn_comp_tp) ; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Global Field Power (GFP) in rare vs. frequent conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
show_sem = 1 ; 

x_to_plot = [time_vec_used; time_vec_used] ;           
col_mean = [0, 128, 255; 102, 204, 0]./255 ; col_std = col_mean ; 

GFP_I1 = squeeze(compute_GFP(y_I1_TP(:,:,:,i_TPs_compare))) ; 
mean_GFP_I1 = (squeeze(mean(GFP_I1,2)))' ;  mean_GFP_I1 = mean_GFP_I1([2,1],:) ; 
std_GFP_I1 = (squeeze(std(GFP_I1,[],2)))' ; std_GFP_I1 = std_GFP_I1([2,1],:) ; 
if show_sem ; std_GFP_I1 = std_GFP_I1./sqrt(n_subj) ; end
caption_lgd = {'P(I_1)=0.7', 'P(I_1)=0.3'} ;
fn_GFP_I1 = [fn_comp_tp, 'GFP_I1_', str_TPs, str_params] ; 
plot_mean_std(x_to_plot, mean_GFP_I1, std_GFP_I1, 'col_mean', ...
    col_mean, 'x_lab', xlab, 'y_lab', 'GFP','fn_save',...
    fn_GFP_I1, 'fig_fig', fig_fig, 'pdf_fig',pdf_fig, 'eps_fig', eps_fig, ...
    'add_prior', 0, 'shaded_std', 1, 'show_max', 0, 'col_std', col_std, ...
    'caption_lgd', caption_lgd, 'logx_curve', 0, 'xlim_val', xlim_butterfly) ;


GFP_I2 = squeeze(compute_GFP(y_I2_TP(:,:,:,i_TPs_compare))) ; 
mean_GFP_I2 = (squeeze(mean(GFP_I2,2)))' ; % 2 x n_times
std_GFP_I2 = (squeeze(std(GFP_I2,[],2)))' ;
if show_sem ; std_GFP_I2 = std_GFP_I2./sqrt(n_subj) ; end
caption_lgd = {'P(I_2)=0.7', 'P(I_2)=0.3'} ;
fn_GFP_I2 = [fn_comp_tp, 'GFP_I2_', str_TPs, str_params] ; 
plot_mean_std(x_to_plot, mean_GFP_I2, std_GFP_I2, 'col_mean', ...
    col_mean, 'x_lab', xlab, 'y_lab', 'GFP','fn_save',...
    fn_GFP_I2, 'fig_fig', fig_fig, 'pdf_fig',pdf_fig, 'eps_fig', eps_fig, ...
    'add_prior', 0, 'shaded_std', 1, 'show_max', 0, 'col_std', col_std, ...
    'caption_lgd', caption_lgd, 'logx_curve', 0, 'xlim_val', xlim_butterfly) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Avg EEG curves (subfigs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig_sz = [1 2 32 12] ;
taille = 20 ; 
fig_I1 = figure('units','centimeters','outerposition',fig_sz,...
    'Name', 'I1') ; 
fig_I2 = figure('units','centimeters','outerposition',fig_sz,...
    'Name', 'I2') ; 
idx_subplot = 1 ; 
for i_tp = i_TPs_compare
    y_I1_avg = squeeze(mean(y_I1_TP(:,:,:,i_tp),3)) ; 
    y_I2_avg = squeeze(mean(y_I2_TP(:,:,:,i_tp),3)) ; 
    
    if selected_chan
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
    
    % --*-- ----> I_1
    figure(fig_I1) ; %hold on ; 
    h_I1(idx_subplot) = subplot(1, 2, idx_subplot) ; hold on ; 
    plot_butterfly(time_vec_used, y_I1_avg, 'xlim_val', ...
        xlim_butterfly, 'x_lab', xlab, 'topo_colors', ...
        topo_colors, 'chanlocs', chanlocs_butterfly, 'eps_fig', eps_fig, ...
        'fig_fig', fig_fig, 'pdf_fig', pdf_fig, 'create_fig', 0, 'save_fig', 0) ;
    title(strrep(all_TPs_str{i_tp},'_', '\_'), 'FontSize', taille,'Interpreter','Latex') ; 
    
    % --*-- ----> I_2
    figure(fig_I2) ; %hold on ; 
    h_I2(idx_subplot) = subplot(1, 2, idx_subplot) ; hold on ; 
    plot_butterfly(time_vec_used, y_I2_avg, 'xlim_val', ...
        xlim_butterfly, 'x_lab', xlab, 'topo_colors', ...
        topo_colors, 'chanlocs', chanlocs_butterfly, 'eps_fig', eps_fig, ...
        'fig_fig', fig_fig, 'pdf_fig', pdf_fig, 'create_fig', 0, 'save_fig', 0) ; 
    title(strrep(all_TPs_str{i_tp},'_', '\_'), 'FontSize', taille,'Interpreter','Latex') ; 
    
    idx_subplot = idx_subplot + 1 ; 
end
linkaxes(h_I1,'y') ;
linkaxes(h_I2,'y') ;
fn_I1 = [fn_comp_tp, 'Butterfly_I1_', str_TPs, str_params,str_plot] ; 
fn_I2 = [fn_comp_tp, 'Butterfly_I2_', str_TPs, str_params,str_plot] ; 
if topo_colors
    fn_I1 = [fn_I1, '_topo'] ; fn_I2 = [fn_I2, '_topo'] ; 
end

%reduce_hspace_subplot(h_I1, 0.6, 1) ; 
% annotation('textbox', [0 0.07 1 0.06], ... % [0 0.02 1 0.05]
%     'String', xlab,'EdgeColor', 'none', ...
%     'HorizontalAlignment', 'center','FontSize',20, 'Interpreter', 'latex') ;
my_save_fig(fig_I1, fn_I1, eps_fig, fig_fig, pdf_fig) ; 

%reduce_hspace_subplot(h_I2, 0.6, 1) ; 
% annotation('textbox', [0 0.07 1 0.06], ... % [0 0.02 1 0.05]
%     'String', xlab,'EdgeColor', 'none', ...
%     'HorizontalAlignment', 'center','FontSize',20, 'Interpreter', 'latex') ;
my_save_fig(fig_I2, fn_I2, eps_fig, fig_fig, pdf_fig) ;

%set(gca,'YLim', ylim_val) ;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO: for GFP
% ** can use plot_mean_std(x_data, y_data, y_std, varargin) for a single
% elec and to find the peaks!!

return
% ------------------------------------------------------------------- %
% --*-- Plot all avg per param (TP)
% ------------------------------------------------------------------- %

for i_tp = 1:n_TPs
    y_I1_avg = squeeze(mean(y_I1_TP(:,:,:,i_tp),3)) ; 
    y_I2_avg = squeeze(mean(y_I2_TP(:,:,:,i_tp),3)) ; 
    
    if selected_chan
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
    
    % --*-- ----> I_1
    if selected_chan
        fn_res_I1 = [fn_res_tp, num2str(n_chans_plot), 'chans_I1/'] ; 
        if ~exist(fn_res_I1, 'dir') ;  mkdir(fn_res_I1) ; end
    else
        fn_res_I1 = fn_res_tp ; 
    end
    fn_I1 = [fn_res_I1, 'Butterfly_I1', all_TPs_str{i_tp}, str_params,str_plot] ; 
    plot_butterfly(time_vec_used, y_I1_avg, 'fn_save', fn_I1, 'xlim_val', ...
        xlim_butterfly, 'x_lab', xlab, 'topo_colors', ...
        topo_colors, 'chanlocs', chanlocs_butterfly, 'eps_fig', eps_fig, ...
        'fig_fig', fig_fig, 'pdf_fig', pdf_fig) ; 
    
    % --*-- ----> I_2
    if selected_chan
        fn_res_I2 = [fn_res_tp, num2str(n_chans_plot), 'chans_I2/'] ; 
        if ~exist(fn_res_I2, 'dir') ;  mkdir(fn_res_I2) ; end
    else
        fn_res_I2 = fn_res_tp ; 
    end
    fn_I2 = [fn_res_I2, 'Butterfly_I2', all_TPs_str{i_tp}, str_params,str_plot] ; 
    plot_butterfly(time_vec_used, y_I2_avg, 'fn_save', fn_I2, 'xlim_val', ...
        xlim_butterfly, 'x_lab', xlab, 'topo_colors', ...
        topo_colors, 'chanlocs', chanlocs_butterfly, 'eps_fig', eps_fig, ...
        'fig_fig', fig_fig, 'pdf_fig', pdf_fig) ;
end





end

