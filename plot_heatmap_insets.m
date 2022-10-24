function [] = plot_heatmap_insets(fig_name, x_vec, y_vec, all_data, cmap, ...
    label_str, y_label, x_label, x_vec_R, mean_R, std_R, x_vec_T, mean_T, std_T, ...
    varargin)
% Display a colormap with Right and Top insets.
%
% x_vec_R, mean_R, std_R: data for the Right inset. 
% x_vec_T, mean_T, std_T: data for the Top inset.

% used in TSL_fit_on_ratings.

clim_val    = NaN ; 
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
log_xscale  = 0 ;   %
log_yscale  = 0 ;   % 

% Colorbar opts
hc_ticks = NaN ; 
hc_ticklabels = NaN ; 
xticks = [] ; 
yticks = [] ; 
xtick_labs = [] ; 
ytick_labs = [] ; 

% insets params
col_std     = 0.6*ones(1,3) ; 
lwdt_mean   = 2 ; 
lwdt        = 2 ; % prior
taille_stick= 16 ; 
add_prior   = 0 ; 
col_prior   = 0.6*ones(1,3) ; 
sz_pt       = 49 ; % 25 
show_max    = 0 ; 
prior_val   = NaN ; 

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
            case 'show_max'
                show_max = Val_arg ; 
            case 'add_prior'
                add_prior = Val_arg ; 
            case 'prior_val'
                prior_val = Val_arg ; 
        end
    end
end

opt_xscale = 'linear' ; opt_yscale = 'linear' ; 
if log_xscale ; opt_xscale = 'log' ; end
if log_yscale ; opt_yscale = 'log' ; end

if isnan(all_data2)
    two_maps = 0 ; 
end

taille_tit = 28 ; %20 ; 
taille = 20 ; 
taille_lab = 28 ; %26 ; % when latex char, smaller
fig_size = [5 5 18 14] ; 
if create_fig
    curr_fig = figure('Name', fig_name, 'units','centimeters',...
         'outerposition',[5 5 18 14]) ; hold on ; 
end

n_rows = 3 ; 
n_cols = 4 ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Heatmap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idx_map = [(n_cols+1):(n_cols+2), (2*n_cols+1):(2*n_cols+2)] ; 
h_m = subplot(n_rows, n_cols, idx_map) ; hold on ; 

% when limiting caxis (with imagesc or caxis), we set the min/max color 
% limits so that values of clims(1) OR LESS map to the first color in the 
% colormap and values of clims(2) OR MORE map to the last color
switch map_opt
    case 1
        if log_xscale || log_yscale
            % cannot use imagesc with an axis set to log scale (this would
            % change the values on the axes). 
            surf(double(x_vec), double(y_vec), all_data); 
            if ~isnan(clim_val)
                caxis(clim_val) ; 
            end
            set(gca, 'XLim', [x_vec(1),x_vec(end)], 'YLim', [y_vec(1), y_vec(end)]) ; 
            view(2) ; 
        else        
            if ~isnan(clim_val)
                imagesc(double(x_vec), double(y_vec), all_data, clim_val) ; 
            else
                imagesc(double(x_vec), double(y_vec), all_data) ; 
                % if with a mask: [clims_PC(1)-dmap, clims_PC(2)]
            end
        end
    case 2
        contourf(double(x_vec), double(y_vec), all_data, levels) ; 
        if ~isnan(clim_val)
            caxis(clim_val) ; 
        end
    case 3
        if two_maps
            [C,h(1)] = contour(double(x_vec), double(y_vec), all_data,levels,...
                'color',cmap(1,:),'LineWidth', 2);
        else
            [C,h] = contour(double(x_vec), double(y_vec), all_data, levels,...
                'LineWidth', 2);
            if ~isnan(clim_val)
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

shading flat ; 
set(gca,'YDir','Normal', 'XScale', opt_xscale, 'YScale', opt_yscale, ...
    'fontsize',taille) ;
% YDir to 'reverse'. Values along the y-axis increase from top to bottom. 
% To decrease the values from top to bottom, set YDir to 'normal'. 
% This setting reverses both the y-axis and the image.

if ~isnan(clim_val) 
    caxis(clim_val) ;   
    % Caxis rescales the colormap so anything below caxis(1) is the minimal 
    % color and anything above caxis(2) is the top of the colormap
end

% set(gca,'fontsize',taille,'XGrid','on','YGrid','on','XMinorGrid','on',...
%     'YMinorGrid','on','layer','bottom');
ylabel(y_label, 'FontSize', taille_lab, 'Interpreter', 'Latex') ;
xlabel(x_label, 'FontSize', taille_lab, 'Interpreter', 'Latex') ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Top inset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idx_top = 1:2 ; 
h_T = subplot(n_rows, n_cols, idx_top) ; hold on ; 

if add_prior
    if isnan(prior_val)
        y_tmp = ones(1,2)./length(mean_T) ; 
    else
        y_tmp = ones(1,2).*prior_val ; 
    end
    plot([x_vec(1),x_vec(end)], y_tmp, 'linewidth', lwdt, 'Color',col_prior) ;  hold on ; 
    
    yL_T = [min(y_tmp(1), min(mean_T-std_T)), max(y_tmp(1), max(mean_T+std_T))] ; 
else
    yL_T = [min(mean_T-std_T), max(mean_T+std_T)] ; 
end

fill([x_vec_T, fliplr(x_vec_T)], [mean_T-std_T, fliplr(mean_T+std_T)], col_std, ...
    'EdgeColor','none', 'FaceAlpha', 0.2) ; hold on ;
plot(x_vec_T, mean_T, '-','color', 'k', 'Linewidth', lwdt_mean) ;

if show_max    
    [y_max, max_idx] = max(mean_T) ;
    scatter(x_vec_T(max_idx), y_max, sz_pt, 'o', 'LineWidth', 1.5, ...
        'MarkerEdgeColor', [0,0,0], 'MarkerFaceColor', col_std) ; hold('on') ;
end

set(gca,'YGrid', 'on', 'XGrid','off','XMinorGrid','off', ...
    'FontSize', taille_stick, 'XLim', [x_vec(1),x_vec(end)], ...
    'YLim', yL_T,'XScale', opt_xscale) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Right inset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ATTENTION: plot x data along the y-axis (in normal dir, same as map)
idx_right = [n_cols+3, 2*n_cols+3] ; 
h_R = subplot(n_rows, n_cols, idx_right) ; hold on ; 

if add_prior
    if isnan(prior_val)
        y_tmp = ones(1,2)./length(mean_R) ; 
    else
        y_tmp = ones(1,2).*prior_val ; 
    end
    plot(y_tmp, [y_vec(1), y_vec(end)], 'linewidth', lwdt, 'Color',col_prior) ; hold on ; 
    xL_R = [min(y_tmp(1), min(mean_R-std_R)), max(y_tmp(1), max(mean_R+std_R))] ; 
else
    xL_R = [min(mean_R-std_R), max(mean_R+std_R)] ; 
end

fill([mean_R-std_R, fliplr(mean_R+std_R)], [x_vec_R, fliplr(x_vec_R)], col_std, ...
    'EdgeColor','none', 'FaceAlpha', 0.2) ; hold on ;
plot(mean_R, x_vec_R, '-','color', 'k', 'Linewidth', lwdt_mean) ;

if show_max    
    [y_max, max_idx] = max(mean_R) ;
    scatter(y_max, x_vec_R(max_idx), sz_pt, 'o', 'LineWidth', 1.5, ...
        'MarkerEdgeColor', [0,0,0], 'MarkerFaceColor', col_std) ; hold('on') ;
end

set(gca,'YGrid', 'off', 'XGrid','on','XMinorGrid','off', ...
    'FontSize', taille_stick, 'YLim', [y_vec(1), y_vec(end)], ...
    'XLim', xL_R, 'YScale', opt_yscale) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Adjust the position of the subplots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% TO DO: cfr TCS2_plot_TF_maps
% shift figure 
pos_m = get(h_m, 'Position') ;
pos_R = get(h_R, 'Position') ; 
pos_T = get(h_T, 'Position') ; 

% Shrink map to have the hc label
pos_m(3) = 0.94*pos_m(3) ; 
% upwards
pos_m(2) = pos_m(2) + 0.14*pos_m(4); 
% to right
pos_m(1) = pos_m(1) + 0.08*pos_m(3) ; 

% ===*=== Top inset
pos_T(3) = pos_m(3) ; pos_T(1) = pos_m(1) ; 
pos_T(2) = (pos_m(2)+pos_m(4)) + 0.1*pos_T(4) ;

ticks_tmp = get(h_T, 'YTickLabel') ; 
ticks_tmp{1} = '' ; 
set(h_T,'XTickLabel','') ; %, 'YTickLabel', ticks_tmp) ; 

% ===*=== Right inset
pos_R(4) = pos_m(4) ; pos_R(2) = pos_m(2) ; 
pos_R(1) = pos_m(1)+pos_m(3) + 0.1*pos_R(3) ; 
pos_R(3) = pos_T(4)*fig_size(4)/fig_size(3) ; 

ticks_tmp = get(h_R, 'XTickLabel') ; 
ticks_tmp{1} = '' ; 
set(h_R,'YTickLabel','') ; %, 'XTickLabel', ticks_tmp) ; 


set(h_m,'Position',pos_m) ; 
set(h_R,'Position',pos_R) ; 
set(h_T,'Position',pos_T) ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if two_maps
    lgd = legendflex(h, label_str, 'ncol',2,...
        'anchor',{'ne','ne'},'box','on', 'Interpreter','LateX','fontsize',taille) ;
    set(lgd,'color','w'); %'none'
else
    hc = colorbar('peer', h_m, 'Ticklength', 0.05, 'Position', ...
        [pos_R(1)+1.1*pos_R(3), pos_R(2), 0.2*pos_R(3), pos_R(4)]) ;
    % OR
    % pos_fig = get(gca, 'Position') ;
    % mid_ypos = pos_fig(2)+pos_fig(4)*0.5 ; 
    % hc = colorbar('peer', gca, 'Position', ...
    %     [pos_fig(1)+pos_fig(3)+0.01,  ...
    %     pos_fig(2), 0.02, pos_fig(4)],'Ticklength', 0.05) ;

    hc.Label.Interpreter = 'latex';
    hc.Label.String = label_str ;
    %title(hc, label_str, 'FontSize', taille_tit) ;%  {label_str; ' '}
    hc.Label.FontSize = taille_tit ;
    if ~isnan(hc_ticks)
        hc.Ticks = hc_ticks ;
        if ~isnan(hc_ticklabels)
            hc.TickLabels = hc_ticklabels ;
        else
            hc.TickLabels = hc_ticks ;
        end
    end
end

%hc.TickLabelInterpreter = 'tex' ; 
if ~sum(isnan(clim_val)) && ~two_maps
    ylim(hc,clim_val) ;      
    % ylim only affects what is shown on the colorbar but not the range
    % of colors displayed on the figure
end

% == * == Saving
if create_fig && sum(~isnan(filename))
    my_save_fig(curr_fig, filename, eps_fig, fig_fig, pdf_fig) ;
end


end

