
function plot_butterfly(x_data, y_data, varargin)
% Plot the signals in y_data (each row = 1 signal).
% used in TSL_analyze_EEG.
%
% Inputs
%   x_data: (1, n_pts) matrix of x values
%   y_data: (n_curves, n_pts) matrix

[n_curves,n_data] = size(y_data) ;

% === Default arguments =========== %
eps_fig = 1 ; fig_fig = 0 ; pdf_fig = 0 ; 
fn_save         = 'test' ;          %
save_fig        = 1 ;               %
fig_name        = '' ;              %
x_lab           = 'time after stimulus (s)' ; 
y_lab           = 'amplitude (\muV)' ; 
add_HL          = 0 ;               % horizontal line of color col_HL
neg_pos_HL      = 0 ;               % add at pm HL_val
col_HL          = 0.6*ones(1,3) ;   %
lwdt_HL         = 2 ;               %
lwdt            = 1.5 ;             %
taille_tick     = 20 ;              %
taille_stick    = 16 ;              % small ticks
fig_sz          = [1 2 18 12] ; %[1 2 12 12] ;
logx_curve      = 0 ; 
topo_colors     = 1 ;               % Show the color codes on a scalp topography inset
order_top_colors= 1 ;               % order the channels as if they were to be plotted on topo
cmap            = jet(n_curves) ;   % 
create_fig      = 1 ;               % 
auto_order      = 1 ;               % automate the channel ordering 
% ================================= %


xlim_val = [min(x_data(:)), max(x_data(:))] ; 
HL_val = 1/(n_data*n_curves) ; 

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
            case 'fn_save'
                fn_save = Val_arg ;
            case 'save_fig'
                save_fig = Val_arg ;
            case 'eps_fig'
                eps_fig = Val_arg ; 
            case 'fig_fig'
                fig_fig = Val_arg ; 
            case 'pdf_fig'
                pdf_fig = Val_arg ; 
            case 'fig_name'
                fig_name = Val_arg ; 
            case 'x_lab'
                x_lab = Val_arg ; 
            case 'y_lab'
                y_lab = Val_arg ; 
            case 'fig_sz'
                fig_sz = Val_arg ; 
            case 'hl_val'
                HL_val = Val_arg ; 
            case 'neg_pos_hl'
                neg_pos_HL = Val_arg ; 
            case 'taille_tick'
                taille_tick = Val_arg ; 
            case 'taille_stick'
                taille_stick = Val_arg ; 
            case 'add_hl'
                add_HL = Val_arg ; 
            case 'col_hl'
                col_HL = Val_arg ; 
            case 'cmap'
                cmap = Val_arg ; 
            case 'logx_curve'
                logx_curve = Val_arg ; 
            case 'topo_colors'
                topo_colors = Val_arg ; 
            case 'xlim_val'
                xlim_val = Val_arg ; 
            case 'chanlocs'
                chanlocs = Val_arg ; % required if topo_colors
            case 'create_fig'
                create_fig = Val_arg ; 
            case 'fig_h'
                fig_h = Val_arg ; 
            case 'order_top_colors'
                order_top_colors = Val_arg ;
            case 'auto_order'
                auto_order = Val_arg ; 
        end
    end
end

if create_fig
    fig_h = figure('units','centimeters','outerposition',fig_sz,...
        'Name', fig_name) ; 
end
hold on ; 
main_ax = gca ; 
if auto_order
    set(main_ax, 'ColorOrder', flipud(cmap));
else
    set(main_ax, 'ColorOrder', cmap) ; % 'NextPlot', 'replacechildren'); 
end

if add_HL
    x_tmp = xlim_val ;
    y_tmp = HL_val*ones(1,2) ; 
    plot(x_tmp, y_tmp, 'linewidth', lwdt_HL, 'Color',col_HL) ; hold on ; 
    if neg_pos_HL
        plot(x_tmp, -y_tmp, 'linewidth', lwdt_HL, 'Color',col_HL) ; hold on ; 
    end
end

% ====*==== Indicate signals
if topo_colors || order_top_colors
    % Order desired for the topo colormap
   desired_chans = {'Fpz', 'Fz', 'FCZ', 'CZ', 'CPZ', 'PZ', 'POz', 'Oz', ...
        'O1', 'PO3', 'PO7', 'P7', 'P5', 'P3', 'P1', 'CP1', 'CP3', 'CP5', 'TP7', ...
        'T7', 'C5', 'C3', 'C1', 'FC1', 'FC3', 'FC5', 'FT7', 'F7', 'F5', 'F3', ...
        'F1', 'AF3', 'AF7', 'Fp1', 'Fp2', 'AF8', 'AF4', 'F2', 'F4', 'F6', 'F8',...
        'FT8', 'FC6', 'FC4', 'FC2', 'C2', 'C4', 'C6', 'T8', 'TP8', 'CP6', ...
        'CP4', 'CP2', 'P2', 'P4', 'P6', 'P8', 'PO8', 'PO4', 'O2'} ; 
    chan_order = define_chan_order_colors(chanlocs, desired_chans) ;
    
    if auto_order
        %y_data = y_data(chan_order,:) ; 
        y_data = y_data(fliplr(chan_order),:) ; 
        % plot FCZ on top of other curves; but keep same cmap, cfr flipud(cmap)
    end
end
hold on ; 

% Ensure that colormap starts at first color
% to match colors indicated on topoplot
set(gca,'ColorOrderIndex',1) ;

plot(x_data,y_data, '-', 'Linewidth', lwdt) ; hold on ; 
  
if logx_curve
   set(gca,'XScale','log') ;  
end
set(main_ax,'YGrid', 'off', 'XGrid','off','XMinorGrid','off', ...
    'FontSize', taille_stick, 'XLim', xlim_val) ;
xlabel(x_lab, 'FontSize', taille_tick,'FontName','Arial') ; %'Interpreter','Latex'); 
ylabel(y_lab, 'FontSize', taille_tick,'FontName','Arial') ; %'Interpreter','Latex'); 

if topo_colors
    % ==*== Positioning axes for topoplot (cfr TCS2_plot_fft.m)
    set(main_ax, 'Units', 'normalized');
    pos = get(main_ax,'position') ;
    x_full = pos(1) ; w_full = pos(3) ; y_full = pos(2) ; h_full = pos(4) ; 
    
    pos_for_topo = [x_full+0.8*w_full, y_full+0.7*h_full, 0.2*w_full, 0.3*h_full] ; 
    %pos_for_topo = [x_full+0.5*w_full, y_full+0.5*h_full, 0.5*w_full, 0.5*h_full] ; 
    h_topo = axes('position',pos_for_topo,'visible','off') ; hold on ; 
    
    % ==*== Do the topoplot with the colors (cfr test_scalp_plot). 
    n_chan = length(chanlocs) ;
    marker_def = {'.','none', 16, 1} ; % {markerchar color size linewidth}
    
    cmap_cell = cell(1,n_chan) ;
    for idx_chan=1:n_chan
        cmap_cell{idx_chan} = cmap(idx_chan,:) ;
    end
    
    % plotted values: used to select color within the colormap
    if auto_order
        chanlocs_plot = chanlocs(chan_order) ;
    else
        chanlocs_plot = chanlocs ;
    end
    
    my_topoplot(1:n_chan, chanlocs_plot, 'electrodes', 'on', 'style','blank',...
        'emarker', marker_def, 'emarkercolors', cmap_cell, 'emarkersizemark',...
        16) ;
    %  cb = colorbar ;
    %  caxis([1 n_chan+1]) ;% sets colorbar limits
    %  set(cb,'FontSize',18)
    
    axis tight ;
end

set(gcf,'CurrentAxes',main_ax) ; 


if save_fig
    if topo_colors
        fn_save = [fn_save, '_topo'] ; 
    end
    my_save_fig(fig_h, fn_save, eps_fig, fig_fig, pdf_fig) ;
end

end

