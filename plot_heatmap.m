function [] = plot_heatmap(fig_name, x_vec, y_vec, all_data, cmap, ...
    label_str, y_label, x_label, varargin)
% Display a colormap. 

% used in TSL_fit_on_ratings.

clim_val    = NaN ; % for the heatmap (+ hc if hc_lim not specified)
hc_lim      = NaN ; % for the hc only
filename    = NaN ; 
map_opt     = 1 ;   % 1: imagesc
                    % 2: contourf
                    % 3: contour
levels      = 10 ;  % when map_opt in {2,3}
                    % can be an integer OR a vector of increasing level 
                    % values.
two_maps    = 0 ;   % show two datasets --> only with opt_map=3                    
all_data2   = NaN ; 
create_fig  = 1 ;   % create figure 
fig_size    = [5 10 15 12] ;
log_xscale  = 0 ;   %
log_yscale  = 0 ;   % 
vline       = 0 ; % vertical white dotted line at x=0 

% Colorbar opts
hc_ticks = NaN ; 
hc_ticklabels = NaN ; 
ylim_v = NaN ; 
xlim_v = NaN ; 
xticks = [] ; 
yticks = [] ; 
xtick_labs = [] ; 
ytick_labs = [] ; 

 

eps_fig = 1 ; fig_fig = 0 ; pdf_fig = 0 ; 

% Parse varargin
nargs_left = length(varargin) ; 
if nargs_left>0 
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
            case 'clim_val'
                clim_val = Val_arg ;
            case 'filename'
                filename = Val_arg ; 
            case 'map_opt'
                map_opt = Val_arg ; 
            case 'levels'
                levels = Val_arg ; 
            case 'hc_ticks'
                hc_ticks = Val_arg ; 
            case 'hc_ticklabels'
                hc_ticklabels = Val_arg ; 
            case 'two_maps'
                two_maps = Val_arg ; 
            case 'all_data2'
                all_data2 = Val_arg ; 
            case 'create_fig'
                create_fig = Val_arg ; 
            case 'log_xscale'
                log_xscale = Val_arg ; 
            case 'log_yscale'
                log_yscale = Val_arg ; 
            case 'xticks'
                xticks = Val_arg ; 
            case 'yticks'
                yticks = Val_arg ; 
            case 'xtick_labs'
                xtick_labs = Val_arg ; 
            case 'ytick_labs'
                ytick_labs = Val_arg ; 
            case 'eps_fig'
                eps_fig = Val_arg ; 
            case 'fig_fig'
                fig_fig = Val_arg ; 
            case 'pdf_fig'
                pdf_fig = Val_arg ; 
            case 'fig_size'
                fig_size = Val_arg ; 
            case 'vline'
                vline = Val_arg ;
            case 'hc_lim'
                hc_lim = Val_arg ; 
            case 'ylim_v'
                ylim_v = Val_arg ; 
            case 'xlim_v'
                xlim_v = Val_arg ; 
        end
    end
end

if sum(isnan(cmap(:)))
    if length(levels)==1
        n_levels = levels ; 
    else
        n_levels = length(levels) ; 
    end
    switch map_opt
        case 1
            cmap = parula ;
        case 2
            cmap = magma(n_levels-1) ; 
        case 3
            cmap = magma(n_levels) ; 
    end
end

opt_xscale = 'linear' ; opt_yscale = 'linear' ; 
if log_xscale ; opt_xscale = 'log' ; end
if log_yscale ; opt_yscale = 'log' ; end

if isnan(ylim_v) ; ylim_v = [y_vec(1), y_vec(end)] ; end
if isnan(xlim_v) 
    %xlim_v = [x_vec(1),x_vec(end)] ; 
    gap = diff(x_vec(1:2)) ; 
    xlim_v = [x_vec(1)-gap/2,x_vec(end)+gap/2] ; 
end

if sum(isnan(all_data2))
    two_maps = 0 ; 
end

taille_tit = 28 ; %20 ; 
taille = 20 ; 
taille_lab = 24 ; %28 ; %26 ; % when latex char, smaller
if create_fig
    curr_fig = figure('Name', fig_name, 'units','centimeters',...
         'outerposition',fig_size) ; 
end
% when limiting caxis (with imagesc or caxis), we set the min/max color 
% limits so that values of clims(1) OR LESS map to the first color in the 
% colormap and values of clims(2) OR MORE map to the last color
switch map_opt
    case 1
        if log_xscale || log_yscale
            % cannot use imagesc with an axis set to log scale (this would
            % change the values on the axes). 
            surf(double(x_vec), double(y_vec), all_data); 
            if ~sum(isnan(clim_val))
                caxis(clim_val) ; 
            end
            set(gca, 'XLim', xlim_v, 'YLim', ylim_v) ; 
            view(2) ; 
        else        
            if ~sum(isnan(clim_val))
                imagesc(double(x_vec), double(y_vec), all_data, clim_val) ; 
            else
                imagesc(double(x_vec), double(y_vec), all_data) ; 
                % if with a mask: [clims_PC(1)-dmap, clims_PC(2)]
            end
            
            set(gca, 'XLim', xlim_v, 'YLim', ylim_v) ; % [xlim_v(1)-gap/2,xlim_v(end)+gap/2]
        end
    case 2
        contourf(double(x_vec), double(y_vec), all_data, levels) ; 
        if ~sum(isnan(clim_val))
            caxis(clim_val) ; 
        end
    case 3
        if two_maps
            [C,h(1)] = contour(double(x_vec), double(y_vec), all_data,levels,...
                'color',cmap(1,:),'LineWidth', 2);
        else
            [C,h] = contour(double(x_vec), double(y_vec), all_data, levels,...
                'LineWidth', 2);
            if ~sum(isnan(clim_val))
                caxis(clim_val) ; 
            end
        end
        
        
        if ~max(isnan(hc_ticklabels)) && isfloat(hc_ticklabels)
            clabel(C,h, hc_ticklabels,'LabelSpacing',1000,'FontSize',14) ;
        else
            clabel(C,h, 'LabelSpacing',1000) ; %,'LabelSpacing',300) % define space between labels, specified as a scalar value in point units
        end
%         contour(double(x_vec), double(y_vec), all_data, levels,'ShowText','on',...
%             'LineWidth', 2,'LabelSpacing',400)% default LabelSpacing = 144
end
hold on ; 
if vline
   plot([0,0], ylim_v, '--', 'color', [1, 1, 1],'LineWidth', 2) ; 
end

if two_maps
    [C,h(2)] = contour(double(x_vec), double(y_vec), all_data2,levels,...
        'color', cmap(2,:), 'LineWidth', 2) ;
    if ~max(isnan(hc_ticklabels)) && isfloat(hc_ticklabels)
        clabel(C,h(2), hc_ticklabels,'LabelSpacing',1000,'FontSize',14) ;
    else
        clabel(C,h(2), 'LabelSpacing',1000) ;
    end
else
    colormap(cmap) ; hold on ; % 10*log10()
end

if ~isempty(xticks)
    labs_used = xticks ; 
    if ~isempty(xtick_labs)
        labs_used = xtick_labs ;  
    end
    set(gca, 'xtick',xticks, 'xticklabel', labs_used) ; 
end
if ~isempty(yticks)
    labs_used = yticks ; 
    if ~isempty(ytick_labs)
        labs_used = ytick_labs ;  
    end
    set(gca, 'ytick',yticks, 'yticklabel', labs_used) ; 
end

shading flat ; set(gca,'YDir','Normal', 'XScale', opt_xscale, 'YScale', opt_yscale) ;
% YDir to 'reverse'. Values along the y-axis increase from top to bottom. 
% To decrease the values from top to bottom, set YDir to 'normal'. 
% This setting reverses both the y-axis and the image.

if two_maps
%     lgd = legend(label_str,'location','best') ;
%     set(lgd, 'FontSize', taille) ; 
    lgd = legendflex(h, label_str, 'ncol',2,...
        'anchor',{'ne','ne'},'box','on', 'fontsize',taille, 'FontName', 'Arial') ; %'Interpreter','LateX',) ;
    set(lgd,'color','w'); %'none'
else
    hc = colorbar('peer', gca, 'Ticklength', 0.05) ;
    % OR
    % pos_fig = get(gca, 'Position') ;
    % mid_ypos = pos_fig(2)+pos_fig(4)*0.5 ; 
    % hc = colorbar('peer', gca, 'Position', ...
    %     [pos_fig(1)+pos_fig(3)+0.01,  ...
    %     pos_fig(2), 0.02, pos_fig(4)],'Ticklength', 0.05) ;

    hc.Label.Interpreter = 'tex' ; %'latex';
    %hc.Label.String = label_str ;
    title(hc, label_str, 'FontSize', taille_tit) ;%  {label_str; ' '}
    hc.Label.FontSize = taille_tit ;
    if ~max(isnan(hc_ticks))
        hc.Ticks = hc_ticks ;
        if ~max(isnan(hc_ticklabels))
            hc.TickLabels = hc_ticklabels ;
        else
            hc.TickLabels = hc_ticks ;
        end
    end
    if ~isnan(hc_lim)
        ylim(hc,hc_lim) ; 
    end
end

%hc.TickLabelInterpreter = 'tex' ; 
if ~max(isnan(clim_val)) 
    caxis(clim_val) ;
    % Caxis rescales the colormap so anything below caxis(1) is the minimal 
    % color and anything above caxis(2) is the top of the colormap
    if ~two_maps
        ylim(hc,clim_val) ;      
        % ylim only affects what is shown on the colorbar but not the range
        % of colors displayed on the figure
    end
end
% if add_signific
%     ylim(hc,clims_PC)
% else
%     caxis(clims_PC)
%
set(gca,'fontsize',taille,'XGrid','off','YGrid','off','XMinorGrid','off',...
    'YMinorGrid','off','layer','bottom'); % xlim([xlim_down, xlim_up])
ylabel(y_label, 'FontSize', taille_lab, 'FontName', 'Arial') ; %'Interpreter', 'Latex')

xlabel(x_label, 'FontSize', taille_lab, 'FontName', 'Arial') ; 

% shift figure 
pos_fig = get(gca, 'Position') ;
% towards the left to have the hc label
pos_fig(3) = 0.94*pos_fig(3) ; 
% upwards
pos_fig(2) = pos_fig(2) + 0.15*pos_fig(4); %0.05*pos_fig(4); 
pos_fig(4) = 0.8*pos_fig(4) ; 
set(gca,'Position',pos_fig) ; 

% == * == Saving
if create_fig && sum(~isnan(filename))
    my_save_fig(curr_fig, filename, eps_fig, fig_fig, pdf_fig) ;
end


end

