function [P, hw, win_width, Fb] = my_stft_with_conv(x, f, win_width, Fb, fs, plot_cwt)
% Estimate the PSD of x using short-time Fourier transform (STFT). 
% 
% Inputs
%   'x'     the original data samples.
%   'f'     the discrete frequency samples (normalized).
%   'win_width'     window width, in number of samples. 
%   'Fb'    bandwidth parameter of the Gaussian window. 
%   'fs'    only required for the plots (default=1000hz).
% Outputs
%   'P'     STFT matrix (freqs along the rows and time along cols).
%   'hw'    wavelet employed at each frequency for the convolution. 
%   'win_width'     window width, in number of samples.
%   'Fb'    bandwidth parameter of the Gaussian window.
%
% Written based on sub_sep_cwt from STEP toolbox. 
% Dounia Mulders. 

% pre-processing
if size(x,2)==1
    x = x.';
end

if nargin<3 % f0 not specified 
    win_width = 500 ; % window width, in number of samples
    % assuming f_s = 1000Hz, use window of 0.5 seconds
end
if nargin<4
    Fb = 0.05 ;
end

N_F = length(f);
N_T = length(x);
% STFT matrix
P = zeros(N_F ,N_T) ;


L_hw = N_T; % filter length
tau_vec = [-L_hw:L_hw] ; % sufficient length s.t. multiply all input signal for all tau in [0, L_hw]
full_length = length(tau_vec) ; 
hw = zeros(N_F, full_length) ; 
%P_full = zeros(N_F, N_T+ (2*L_hw+1) -1 ) ; 
% length of u: (2*L_hw+1)

for fi=1:N_F
    f0 = 0.5*win_width*f(fi)/(sqrt(Fb*4.5)) ; 
    f_f0 = sqrt(Fb*4.5)/(0.5*win_width) ; % f/f0

    u = f_f0 * -(tau_vec);
    hw(fi,:) = sqrt(f_f0) *((pi*Fb)^(-0.5))*exp(2*1i*pi*f0*u).*exp(-(u.*u)/Fb);

    %P_full(fi,:) = conv(x,conj(hw(fi,:))); % length = N_T + full_length - 1
    %P(fi,:) = P_full(fi,L_hw+1:L_hw+N_T);   % tau = 1:N_T        
    P(fi, :) = conv(x,conj(hw(fi,:)), 'same');
end

if nargin<6
    plot_cwt = 0 ; 
end
if plot_cwt
    if nargin<5
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
        line([sqrt(4.5)*sqrt(Fb)*f0/f_tmp sqrt(4.5)*sqrt(Fb)*f0/f_tmp],...
            yL,'lineStyle', ':', 'Color',[0 0 0], 'linewidth',1) ; hold on
        line([-sqrt(4.5)*sqrt(Fb)*f0/f_tmp -sqrt(4.5)*sqrt(Fb)*f0/f_tmp],...
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
    % ----> Filtered signal (convolution (wavelet, original signal)). --
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
