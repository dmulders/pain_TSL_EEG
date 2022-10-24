
function plot_bars(data_mat, varargin)
% Plot bars of input data.
% used in TSL_fit_on_ratings (exceedance proba).
%
% data_mat
%   Matrix (M, N) with M x-positions and N bars per x-position. 

% === Default arguments
eps_fig = 1 ; fig_fig = 0 ; pdf_fig = 0 ; 
fn_save         = 'test' ; 
save_fig        = 1 ; 
fig_name        = '' ; 
x_lab           = 'x' ; 
y_lab           = 'y' ; 
cbar_lab        = 'Index' ; 
add_prior       = 0 ;
col_prior       = 0.6*ones(1,3) ; 
lwdt            = 2 ; 
taille_tick     = 20 ;  %
taille_stick    = 16 ;  % small ticks
fig_sz          = [1 2 18 12] ; %[1 2 12 12] ;
use_colormap    = 1 ;   % colormap encoding heights of bars; 
                        % otherwise, one color per bar (i=1...N)
c_bars          = NaN ; 
add_colbar      = 0 ;   % only if use_colormap  
view_3D         = 1 ;   % if 2D and more than one group, cannot use colormap for height
                        % because we use the 3D plotting for the colormap
caxis_val       = [0,1] ; 

[n_pos,n_gr] = size(data_mat) ; 
prior_val = 1/n_pos ; 
x_descr = cell(1, n_pos) ; 
for i_pos = 1:n_pos    
    x_descr{i_pos} = [num2str(i_pos)] ; 
end

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
            case 'c_bars'
                c_bars = Val_arg ; 
            case 'x_descr'
                x_descr = Val_arg ; 
            case 'x_lab'
                x_lab = Val_arg ; 
            case 'y_lab'
                y_lab = Val_arg ; 
            case 'fig_sz'
                fig_sz = Val_arg ; 
            case 'use_colormap'
                use_colormap = Val_arg ; 
            case 'add_colbar'
                add_colbar = Val_arg ; 
            case 'prior_val'
                prior_val = Val_arg ; 
            case 'taille_tick'
                taille_tick = Val_arg ; 
            case 'taille_stick'
                taille_stick = Val_arg ; 
            case 'view_3d'
                view_3D = Val_arg ; 
            case 'add_prior'
                add_prior = Val_arg ; 
            case 'col_prior'
                col_prior = Val_arg ; 
            case 'caxis_val'
                caxis_val = Val_arg ; 
        end
    end
end

if sum(isnan(c_bars))
    if use_colormap
        c_bars = jet(150) ; %parula ; magma
    else
        %c_bars = jet(n_gr) ; %parula ; magma
        [ ~, c_bars] = get_some_nice_colors(1) ; 
        c_bars = c_bars(1:n_gr,:) ;
    end
end

fig_h = figure('units','centimeters','outerposition',fig_sz,...
    'Name', fig_name) ;

if use_colormap && n_gr==1
    if add_prior
        %plot3([0.5,n_pos+0.5],[0.5,n_pos-0.5], [prior_val, prior_val], ...
        %    '-', 'color', col_prior, 'linewidth', 2) ; 
        x_tmp = [0.5,n_gr+0.5] ; 
        y_tmp = [0.5,n_pos+0.5] ; % ATTENTION: positions are along the y-axis!!
        z_tmp = meshgrid(x_tmp,y_tmp) ; z_tmp(:) = prior_val ; 
        surf(x_tmp, y_tmp, z_tmp, 'linewidth', 2, 'EdgeColor',col_prior, ...
            'faceColor', col_prior, 'faceAlpha', 0.4); 
        hold on ; 
    end
    
    h = bar3(1:n_pos, data_mat) ; 
    if add_colbar
        colorbar ; 
    end
    caxis(caxis_val) ; 
    colormap(c_bars) ; 
    
    if n_gr==1 && ~view_3D
        set(gca,'YDir','reverse') ;
        view(-90,0) ; 
    end
    
    for k = 1:length(h)
        zdata = h(k).ZData ; 
        h(k).CData = zdata ; 
        h(k).FaceColor = 'interp' ; 
    end
    
else
    if view_3D
        if add_prior
            x_tmp = [0.5,n_gr+0.5] ; 
            y_tmp = [0.5,n_pos+0.5] ; % ATTENTION: positions are along the y-axis!!
            z_tmp = meshgrid(x,y) ; z_tmp(:) = prior_val ; 
            surf(x_tmp, y_tmp, z_tmp, 'linewidth', 2, 'EdgeColor',col_prior, ...
                'faceColor', col_prior, 'faceAlpha', 0.4); 
            hold on ; 
        end

        h = bar3(1:n_pos, data_mat) ; 
        caxis(caxis_val) ; 
        colormap(c_bars) ;
        
    else
        if add_prior
            x_tmp = [0.5,n_pos+0.5] ;
            y_tmp = prior_val*ones(1,2) ; 
            plot(x_tmp, y_tmp, 'linewidth', 2, 'Color',col_prior) ; 
            hold on ; 
        end

        h = bar(1:n_pos, data_mat) ;
        caxis(caxis_val) ; 
        colormap(c_bars) ;
        
    end
    
end

if view_3D || use_colormap
    % x_descr are along the y axis
    set(gca,'YGrid', 'off', 'XGrid','off','XMinorGrid','off', ...
        'FontSize', taille_stick, 'YLim', [0.5, n_pos+0.5], 'XLim', ...
        [0.5,n_gr+0.5], 'Ytick',1:n_pos, 'YtickLabel', x_descr) ;
    ylabel(x_lab, 'FontSize', taille_tick,'FontName', 'Arial') ; %'Interpreter','Latex'); 
    zlabel(y_lab, 'FontSize', taille_tick,'FontName', 'Arial') ; % 'Interpreter','Latex'); 
    
else
    % x_descr are along the x axis
    set(gca,'YGrid', 'on', 'XGrid','off','XMinorGrid','off', ...
        'FontSize', taille_stick, 'XLim', [0.5, n_pos+0.5], ...
        'Xtick',1:n_pos, 'XtickLabel', x_descr) ; %, 'FontName', 'Arial') ;
    xlabel(x_lab, 'FontSize', taille_tick,'FontName', 'Arial') ; %'Interpreter','Latex'); 
    ylabel(y_lab, 'FontSize', taille_tick,'FontName', 'Arial') ; % 'Interpreter','Latex'); 
end

if save_fig
    my_save_fig(fig_h, fn_save, eps_fig, fig_fig, pdf_fig) ;
end

end

