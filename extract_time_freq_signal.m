function [all_STFT, all_freq_Hz, all_t_STFT, mean_all_ER_percents,f0,Fb] = ...
    extract_time_freq_signal(in_data, f_samp, x_vec, varargin) 

% Compute the short-time fourier transform along x_vec in time and from 
% f_begin to f_stop Hz in frequency.
% If f_stop < f_begin+f_step, only one frequency is considered. 
%
% In: 
%   in_data:    dataset with the time in 2nd dimension, and possibly
%               epochs or channels in 1st dimension. 
%   f_samp:     sampling frequency (Hz).
%   x_vec:      time vector for in_data.
%   varagin:
%       use_cwt:    extract the phase using CWT.    
%       use_stft:   extract the phase using STFT.
%                   if ~use_cwt && ~use_stft: use bandpass filters and then 
%                   analytic signals.
%       f_begin     default 5Hz.
%       f_stop      default 90Hz.
%       f_step:     default 1Hz. frequencies all_freq_Hz = 
%                   f_begin:f_step:f_stop are considered when use_stft or  
%                   use_cwt. Otherwise, when analytic signals (AS) are  
%                   used, AS are computed in the frequency bands freq_bands.
%       freq_bands: (n_bands,2) matrix giving the frequency bands 
%                   considered when ~use_stft && ~use_cwt.
%                   Default: [all_freq_Hz'-1, all_freq_Hz'+1].  
%                   If it is specified, all_freq_Hz = mean(freq_bands, 2),
%                   and f_begin, f_stop, f_step are NOT used.
%       max_filt_order: maximum order of the filter used for the AS, ie 
%                       when ~use_stft && ~use_cwt. 
%       win_width_sec:  duration (in seconds) of the hann window used if
%                       use_stft.
%   
% Out: 
%   all_STFT: STFT in the considered time-frequency region, across 
%           epochs (one entry per epoch). Dim: [nF, nT, n_epoch]. 

% 3 ways of extracting the phase information: 
% 1) using STFT
% 2) using wavelets
% 3) using a narrow bandpass filter and then computing the analytic signal.
% (after bandass: standardize, cfr PAC_analysis)
%
% adapted from extract_phase_time_freq_plane which is used in 
% PAC_link_phase_ERP.m.
%
% to be used in TCS2_analyze_SSEP.
show_spectrum_AS = 0 ;          % plot the bandpass filtered signals used when AS is used

use_sub_tfa_stft = 0 ;          % true to use the LW function sub_tfa_stft
                                % to compute the PSD (when use_stft!)

% Default varargin %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stand_STFT = 0 ;                % true to standardize the time-frequency
                                % signals in each freq band.
stft_with_conv = 0 ;            % implement the STFT (when use_stft) with a
                                % convolution, as for the CWT. 
use_cwt = 0 ; use_stft = 1 ; 
f_begin = 5 ; f_stop = 90 ; t_begin = -1 ; t_stop = 1.5 ; f_step = 1 ;
f_half_width = 1 ; 
win_width_sec = 0.5 ; % 0.5 sec hann window used 
% hann window used if use_stft 
% --> can choose a smaller one at higher frequencies
% for small frequencies (eg 0.2Hz): have to choose a long window
% ie >=5sec (eg). 
max_filt_order = 100 ;  % only use for the AS, ie when ~use_stft && ~use_cwt
freq_bands_spec = 0 ;   % freq_bands specified?
nargs = nargin ;
% Parse varargin %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargs > 1            %
    nargs_left = nargs - 1 ; 
    if ~(round(nargs_left/2) == nargs_left/2)
        error('-- There should be an even number of varargin arguments')
    end
    for i = 1:2:length(varargin)
        Name_arg = varargin{i} ; 
        Val_arg = varargin{i+1};
        if ~isstr(Name_arg)
            error('Flag arguments must be strings')
        end
        Name_arg = lower(Name_arg);%%%!!!!!!
        switch Name_arg
            case 'use_cwt'
                use_cwt = Val_arg ;
            case 'use_stft'
                use_stft = Val_arg ; 
            case 'f_begin' 
                f_begin = Val_arg ; 
            case 'f_stop' 
                f_stop = Val_arg ; 
            case 't_begin'
                t_begin = Val_arg ;
            case 't_stop'
                t_stop = Val_arg ; 
            case 'f_step'
                f_step = Val_arg ; 
            case 'f_half_width'
                f_half_width = Val_arg ; % for freq_bands
            case 'win_width_sec'
                win_width_sec = Val_arg ; 
            case 'freq_bands'
                freq_bands_spec = 1 ; 
                freq_bands = Val_arg ; 
                all_freq_Hz = mean(freq_bands, 2) ; % mean frequency in each band
            case 'max_filt_order'
                max_filt_order = Val_arg ; 
            case 'stand_stft'
                stand_STFT = Val_arg ;
            case 'stft_with_conv'
                stft_with_conv = Val_arg ; 
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if t_stop<=t_begin
    error('Final time should be higher than initial one')
end
if f_stop<f_begin
    error('Final frequency should be higher than initial one')
end
if ~freq_bands_spec
    all_freq_Hz = f_begin:f_step:f_stop ;
    %freq_bands = [all_freq_Hz'-1, all_freq_Hz'+1] ; 
    freq_bands = [all_freq_Hz'-f_half_width, all_freq_Hz'+f_half_width] ;
end
nF = length(all_freq_Hz) ;              % can be = 1!

ref_I_sec = [-0.4, -0.1] ;              % Reference time interval to 
                                        % compute the ER%.    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default parameters (modified only when use_stft)
[n_epoch, nT] = size(in_data) ; % nT = length(x_vec) ; 
all_t_STFT = x_vec ; 
f0 = NaN ; Fb = NaN ; % set to their real values of use_cwt
time_start = tic ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if use_stft
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameters for the STFT           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % gamma paper GL: 
    % STFT with a fixed 200-ms width Hanning window, chosen to achieve a 
    % good tradeoff between time and frequency resolution in the 
    % range of gamma-band frequencies. The STFT yielded, for each trial, a 
    % complex time–frequency spectral estimate F(t, f) at each point (t, f) 
    % of the time–frequency plane extending from -0.5 to 1.0 s in the time 
    % domain, and from 20 to 150 Hz (in steps of 1 Hz) in the frequency domain
    %
    % Pinprick Bart:
    % 1) STFFT with a fixed Hanning window of 500 ms was used to decompose 
    % the band-pass filtered (0.3–30 Hz) EEG signals. 
    % The explored frequencies ranged from 1 to 30 Hz.    
    % 2) STFFT with a fixed Hanning window of 200 ms to explore the 
    % frequencies ranging from 30 to 100 Hz.
    % 3) to visualize changes in post-stimulus EEG activity relative to the 
    % pre-stimulus baseline period, a baseline correction was applied:
    % average EEG signal between -500 and -100 ms relative to stimulus 
    % onset was subtracted from each post-stimulus time point. This was 
    % performed separately for each estimated frequency (Hu et al., 2014).
    n_samples_win = ceil(win_width_sec*f_samp) ; % nb of samples in the time window
    
    if stft_with_conv
        Fb = 0.05 ;
        % f0 depends on frequency, based on n_samples_win
        f_N = all_freq_Hz./f_samp ; 
        all_STFT = zeros(nF, nT, n_epoch) ; 

        for idx_epoch = 1:n_epoch
            all_STFT(:,:, idx_epoch) = my_stft_with_conv(...
                in_data(idx_epoch, :), f_N, n_samples_win, Fb, f_samp) ; 
        end

        all_PSD = all_STFT.*conj(all_STFT) ; 
        all_PSD = all_PSD./f_samp ;
        
    
    else
        % ==*== ~stft_with_conv ==*==
        %win_width_sec = 0.2 ; % 0.2 sec hann window used 
        % 0.2 seconds = 1 oscillation period at 5Hz 
        % (Useless to consider frequencies way below 5Hz!, or should increase the window width!)
        
        hann_win = gausswin(n_samples_win) ; %hanning(n_samples_win) ; 
        noverlap = min(round(n_samples_win/1.05), n_samples_win-1) ;  % must be < n_samples_win !!
        % larger overlap means more time pts at which the STFT is evaluated
        % if noverlap = n_samples_win-1: each original time pt is considered,
        % except first and last n_samples_win/2 pts!
        %[]: default: to obtain 50% overlap, ie using noverlap = floor(n_samples_win/2)
        % round(n_samples_win/1.05)leads to 237 pts in [-1, 1.5] sec
        % -- take largest nb of samples bec when cwt or AS: all time pts kept

        [~,~, curr_t] = spectrogram(in_data(1, :), ...
                hann_win, noverlap, [all_freq_Hz(1),all_freq_Hz(1)], f_samp) ; 
        % of even if nF=1

        all_t_STFT = curr_t ; % time vector corresponsing to the time 
                              % steps from the STFT.
                              % ATTENTION: should be x_start + curr_t if
                              % x_start neq 0!!!

        % Number of time steps in the considered time window: depends on noverlap ! 
        % Number of frequency steps: fixed by Nfft (or specifying a required freq vector)
        nT = length(all_t_STFT) ;   
     
        all_STFT = zeros(nF, nT, n_epoch) ; 
        if ~use_sub_tfa_stft
            for idx_epoch = 1:n_epoch
                % use spectrogram OR sub_tfa_stft from external in LW
                if nF==1
                    % at least 2 freqs are needed, otherwise 4th arg = NFFT!
                    [tmp_stft, ~, ~] = spectrogram(in_data(idx_epoch, :), ...
                        hann_win, noverlap, [all_freq_Hz,all_freq_Hz], f_samp) ; % [s,f,t]
                    all_STFT(:,:, idx_epoch) = tmp_stft(1,:) ;
                else
                    [all_STFT(:,:, idx_epoch), ~, ~] = spectrogram(in_data(idx_epoch, :), ...
                        hann_win, noverlap, all_freq_Hz, f_samp) ; % [s,f,t]
                    % in fact not the spectrogram but rather the complex STFT...
                end
            end
            all_PSD = all_STFT.*conj(all_STFT) ;
            all_PSD = all_PSD./f_samp./sum(hann_win.^2) ; 
        else
            % using sub_tfa_stft --> gives similar plots!
            [all_STFT, all_PSD] = sub_tfa_stft(in_data', x_vec*1000, all_t_STFT*1000, all_freq_Hz, ...
                f_samp, win_width_sec*1000, 'hann') ; 
        end
    end
    
elseif use_cwt  
    % ------------------------------- %
    % === DEFAULT: leads to 5.69 oscillations
    %f0 = 6 ; 
    %Fb = 0.05 ; 
    % === 3 oscillations:
    %f0 = 3.16 ;   
    %Fb = 0.05 ;
    % === 4 oscillations
    %f0 = 4.22 ; 
    %Fb = 0.05 ;
    % === 5 oscillations
    f0 = 5.27 ; 
    Fb = 0.05 ;
    % === 6 oscillations
    %f0 = 6.32 ; 
    %Fb = 0.05 ;
    % ------------------------------- % 
    f_N = all_freq_Hz./f_samp ; 
    all_STFT = zeros(nF, nT, n_epoch) ; 
    
    for idx_epoch = 1:n_epoch
        all_STFT(:,:, idx_epoch) = my_sub_sep_cwt(...
            in_data(idx_epoch, :), f_N,'normal', f0, Fb, f_samp) ; 
        % [P, hw, Fc, Fb] = my_sub_sep_cwt(x, f, option, f0, Fb)
    end
    
    all_PSD = all_STFT.*conj(all_STFT) ; 
    all_PSD = all_PSD./f_samp ; 
else
    % use narrow bandpass filters and then analytic signals to construct
    % time-frequency signal
    all_STFT = zeros(nF, nT, n_epoch) ; 
    % cfr PAC_analysis (or extract_power, or extract_PTPA)
    
    % Bandpass filtering !
    filt_order = min(max_filt_order, floor(length(x_vec)/3)) ;
    butter_order = 4 ; 
    use_butter = 1 ; % good to isolate thin frequency bands (fir1 ko, cfr show_spectrum_AS)
    
    for idx_freq = 1:nF
        % Consider 2*f_step (default 2Hz) bandwidth if nargin<13
        if use_butter
            [z,p,k] = butter(butter_order,freq_bands(idx_freq,:)./(f_samp/2),...
                'bandpass') ;
            [sos,G] = zp2sos(z,p,k) ;
        else
            disp(['Order used for the bandpass filter: ', num2str(filt_order)])
            [b,a] = fir1(filt_order,...
                freq_bands(idx_freq,:)./(f_samp/2),'bandpass') ;
        end
            
        for idx_epoch=1:n_epoch              
            curr_epoch_signal = in_data(idx_epoch, :) ; 
            if use_butter
                curr_filtsig = filtfilt(sos, G, double(curr_epoch_signal)); 
            else
                curr_filtsig = filtfilt(b, a, double(curr_epoch_signal)); % Zero-phase digital filtering (a=1)
                % consider the analytic signal (otherwise, no phase
                % information).
            end
            
            %curr_filtsig = (curr_filtsig-mean(curr_filtsig))./std(curr_filtsig) ; % zscore
            analytic_sig = hilbert(curr_filtsig) ; 
            all_STFT(idx_freq, :, idx_epoch) = analytic_sig ; 
        end
    end
    all_PSD = all_STFT.*conj(all_STFT) ; 
    all_PSD = all_PSD./f_samp ;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if show_spectrum_AS
        % show the spectrum of the filtered signals
        all_bls = struct() ; names_yticks = cell(nF,1) ; 
        for idx_freq = 1:nF
            curr_freq = all_freq_Hz(idx_freq) ;
            tmp_T_cell = get_cellstr_freq_band([curr_freq, 0], 1000) ;

            curr_label = ['Sig', tmp_T_cell{1}] ;
            all_bls.(curr_label) = (real(squeeze(all_STFT(idx_freq, :,:))))' ;
            names_yticks{idx_freq} = curr_label ; 
        end
        xlim_spec = [0, 1] ;
        remove_dc = 0 ; 
        % only display frequencies after 0.5Hz, to reduce the yLim of each spectrum
        % plot --> cfr plot_signals_spectrum
        scale_range = 30 ; norm_sig = 1 ; add_spectrum = 1 ;  
        fig_tmp = figure('units','normalized','outerposition',...
            [0.1 0.1 0.6 0.6], 'Name', 'AS') ; hold on ; 
        plot_signals_spectrum(all_bls, ...
            x_vec, f_samp, n_epoch, norm_sig, add_spectrum, 0, ...
            [0, x_vec(end)], xlim_spec, remove_dc, scale_range, names_yticks) ;     
    % disp('=== Close figs and press a key to continue')
    % pause ;
    
%     for idx_freq = 1:nF
%         sig_tmp = squeeze(all_STFT(idx_freq, :,:)) ; 
%         [tf_tmp, freqs] = get_spectrum(sig_tmp, fs, 1,1) ;
%         plot(freqs, tf_tmp, '-','LineWidth', 1) ; hold on;
%     end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
time_stop = toc(time_start) ; 
%disp(['Time needed to compute time-freq maps (use_stft=',num2str(use_stft),...
%    ', use_cwt=',num2str(use_cwt),'): ', num2str(time_stop)])

% can standardize each row of all_STFT (for each epoch and frequency
% considered). 
if stand_STFT
    disp('Standardize STFT')
    all_STFT = (all_STFT-repmat(mean(real(all_STFT), 2), 1, nT, 1))./...
        (repmat(std(real(all_STFT), 0,2), 1, nT, 1)) ;  
end

if nargout>3
    mean_all_ER_percents = compute_ER_values(all_PSD, all_freq_Hz, all_t_STFT,...
        f_begin, f_stop,ref_I_sec, t_begin, t_stop) ; 
end

end

function mean_all_ER_percents = compute_ER_values(all_PSD, ...
    all_freq_Hz, all_t_STFT, f_begin, f_stop, ref_I_sec,t_begin, t_stop)
% From an STFT, compute the ER%, cfr gamma paper GL, for each point in the 
% time-freq domain and for each epoch. Data are then averaged over the
% epochs. 

nT = length(all_t_STFT) ; 
nF = length(all_freq_Hz) ; 

% Reference interval
idx_begin_ref = find(all_t_STFT>=ref_I_sec(1), 1) ;  
idx_stop_ref = find(all_t_STFT>=ref_I_sec(2), 1) ; 
if isempty(idx_stop_ref) ; idx_stop_ref = nT ; end

% Time-frequency interval in which the power should be computed
idx_begin_time_I = find(all_t_STFT>= t_begin,1) ; 
idx_stop_time_I = find(all_t_STFT>= t_stop,1) ;
if isempty(idx_stop_time_I) ; idx_stop_time_I = nT ; end

idx_begin_freq_I = find(all_freq_Hz>= f_begin,1) ; 
idx_stop_freq_I = find(all_freq_Hz>= f_stop,1) ;
if isempty(idx_stop_freq_I) ; idx_stop_freq_I = nF ; end


% Refrence values R(f) (mean along time in reference time window)
all_Rf = mean(all_PSD(:,idx_begin_ref:idx_stop_ref,:), 2) ;% all_PSD: (nF, nT, k)
% DO NOT USE SQUEEZE, to be able to use repmat afterwards!

rep_all_Rf = repmat(all_Rf, 1, nT, 1) ; 
all_ER_percents =  ((all_PSD - rep_all_Rf).*100)./rep_all_Rf ; 
% ER%: corresponds to PSD normalized by frequency band to account for the 
% power decay as f increases

all_ER_TF_win = all_ER_percents(idx_begin_freq_I:idx_stop_freq_I, ...
    idx_begin_time_I:idx_stop_time_I,:) ; 

% Plot the time-frequency power map (across-trials average of the ER%).
mean_all_ER_percents = squeeze(mean(all_ER_TF_win, 3)) ;  % mean over epochs


end

