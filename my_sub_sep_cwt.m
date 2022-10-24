function [P, hw, Fc, Fb] = my_sub_sep_cwt(x, f, option, f0, Fb, fs, plot_cwt)
% this sub-function is used to estimate PSD using wavelet transform (the 
% wavelet basis is selected as the type proposed By Tognola)
% 'P' is used for illustration, and 'phi' is used for signal reconstruction
% option: 'normal' means slow but useful for icwt, 
%         'fast' means fast but not usable for icwt, 
%         'phi' only computes the wavelet basis for icwt
                                        
% fprintf('Calculating CWT ... ')

% 'x' - the original data samples
% 'f' - the discrete frequency samples (normalized)
%
% 'fs':     only required for the plots (default=1000hz).

% Copyright (C) 2009 University of Oxford & The University of Hong Kong
% Authors:   Zhiguo Zhang & Li Hu
%            huli@hku.hk

% pre-processing
if size(x,2)==1
    x = x.';
end

% (the wavelet basis was proposed By Tognola)
if nargin<4 % f0 not specified 
    f0 = 6 ; % central frequency
end
Fc = f0 ;
if nargin<5
    Fb = 0.05 ;
end

N_F = length(f);
N_T = length(x);
t = [1:N_T];

P = zeros(N_F ,N_T) ; %[];
phi = [];

% continuous wavelet transform

L_hw = N_T; % filter length
tau_vec = [-L_hw:L_hw] ;
full_length = length(tau_vec) ; 
hw = zeros(N_F, full_length) ; 
P_full = zeros(N_F, N_T+ (2*L_hw+1) -1 ) ; 
% length of u: (2*L_hw+1)

% CWT: implemented by convolution --> not
% straightforward to evaluate it to a reduced
% nb of samples
% Can subsample the original signal...
if (strcmp(option,'normal')==1)||(strcmp(option,'fast')==1)
    for fi=1:N_F
        u = (f(fi)/f0) * -(tau_vec);
        hw(fi,:) = sqrt((f(fi)/f0)) *((pi*Fb)^(-0.5))*exp(2*1i*pi*Fc*u).*exp(-(u.*u)/Fb);
        
        %P_full(fi,:) = conv(x,conj(hw(fi,:)));
        %P(fi,:) = P_full(fi,L_hw+1:L_hw+N_T);   % tau = 1:N_T
        P(fi,:) = conv(x,conj(hw(fi,:)), 'same') ;
    end
end
if nargin<7
    plot_cwt = 0 ; 
end
if plot_cwt
    if nargin<6
       fs = 1000 ;  
    end
    taille = 20 ; 
    xvec_hw = [1:full_length]./fs ; 
    xvec_hw = xvec_hw - mean(xvec_hw) ; % center the vector around 0
    idx_to_keep = ceil(N_T/2):(full_length-ceil(N_T/2)) ; 
    xvec_hw = xvec_hw(idx_to_keep) ; 
    figure() ; lwdt = 1 ; 
    n_freqs_plot = min(N_F,4) ; %floor(N_F/2) ; 
    % ----> Wavelet -----------
    % -------------------------   
    %f_used = f(1:n_freqs_plot)*fs
    subplot(311)
    for idx_freq=1:n_freqs_plot
        hw_tmp = hw(idx_freq,idx_to_keep) ; f_tmp = f(idx_freq)*fs ;
        plot(xvec_hw ,abs(hw_tmp), '-', 'LineWidth', lwdt) ; 
        hold on ; 
        ax = gca ; ax.ColorOrderIndex = max(1,ax.ColorOrderIndex - 1) ;
        plot( xvec_hw ,real(hw_tmp), '--', 'LineWidth', lwdt) ; hold on

        ax = gca ; ax.ColorOrderIndex = ax.ColorOrderIndex - 1 ;
        plot( xvec_hw ,imag(hw_tmp), ':', 'LineWidth', lwdt, 'MarkerSize', 2) ; hold on
               
        yL = get(gca,'YLim');
        line([sqrt(4.5)*sqrt(Fb)*Fc/f_tmp sqrt(4.5)*sqrt(Fb)*Fc/f_tmp],...
            yL,'lineStyle', ':', 'Color',[0 0 0], 'linewidth',1) ; hold on
        line([-sqrt(4.5)*sqrt(Fb)*Fc/f_tmp -sqrt(4.5)*sqrt(Fb)*Fc/f_tmp],...
            yL,'lineStyle', ':', 'Color',[0 0 0], 'linewidth',1) ; hold on    
    end
    title('Wavelet', 'FontSize', taille, 'Interpreter', 'Latex')
    my_lgd = legend('Norm', 'Real part', 'Imaginary part') ; 
    xlim([xvec_hw(1), xvec_hw(end)]) ; 
    %xlim([-0.5, 0.5]) 
    set(ax, 'fontsize', taille-2)
    set(my_lgd,'FontSize',taille-2);      
    
    % ----> Original signal -------
    % -----------------------------
    t_vec = [1:N_T]./fs ;
    subplot(312)
    plot(t_vec, x, '-', 'LineWidth', 1) ; hold on ; xlim([t_vec(1), t_vec(end)])
    yL = get(gca,'YLim');
    title('Signal', 'FontSize', taille, 'Interpreter', 'Latex')
    set(gca, 'fontsize', taille-2)
    % ----> Filtered signal (convolution (wavelet, original signal). --
    % -----------------------------------------------------------------
    subplot(313)
    tmp_ampl = abs(P(1:n_freqs_plot,:)) ; 
    plot(t_vec, tmp_ampl./repmat(max(tmp_ampl,[],2),1,N_T), '-', 'LineWidth', 2) ; hold on
    yL = get(gca,'YLim');
    xlim([t_vec(1), t_vec(end)])
    set(gca, 'fontsize', taille-2)
    %ax = gca ; ax.ColorOrderIndex = ax.ColorOrderIndex - 1 ; 
    %plot(t_vec, phase_band_vec, '--', 'LineWidth', 1) 
    
    disp('OK, pause')
    pause ; 

end
