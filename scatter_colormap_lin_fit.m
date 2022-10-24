
function scatter_colormap_lin_fit(data_struct, varargin)
% Plot data points from different groups. 
% used in TSL_analyze_ratings (lin fits).
%
% data_struct
%   Structure array with one entry per group of data, which will be
%   associated with a marker. Required fields: 
%       'x' and 'y' for the data to plot along the x and y axis respectively
%       'coeffs': coefficients of the linear regression

% === Default arguments
eps_fig = 1 ; fig_fig = 0 ; pdf_fig = 0 ; 
fn_save         = 'test' ; 
save_fig        = 1 ; 
fig_name        = '' ; 
x_lab           = 'x' ; 
y_lab           = 'y' ; 
cbar_lab        = 'Index' ; 
sz_pts          = 25 ; 
markers         = {'o' , '+', 'p', '*', '^', 'x', 's', 'd'} ; 
lwdt            = 2 ; 
taille_tick     = 20 ;  %
taille_stick    = 16 ;  % small ticks
fig_sz          = [1 2 18 12] ; %[1 2 12 12] ;
use_colormap    = 1 ; % Along the *pts* within group. Otherwise, one color per group
c_pts           = NaN ; 
add_lin_fit     = 0 ; % require the 'coeffs' field to data_struct 
single_fit      = 0 ; % if add_lin_fit, only one
all_corr_tit    = 0 ; % show each group correaltion in the plot title
aggregated_corr = 0 ; % if single_fit; compute corr by merging data from all groups!
axis_equal      = 0 ; 
add_identity    = 0 ; 
add_lgd         = 1 ;   
xL              = NaN ; 
create_fig      = 1 ; 
show_data       = 0 ; % data cloud or not
show_box        = 0 ; % boxplot (to bin data according to x_box)
x_box           = [1, 2] ; 

n_gr = length(data_struct) ; 
all_n_pts = zeros(n_gr,1) ; 
caption_lgd = cell(1, n_gr) ; 
for i_gr = 1:n_gr
    tmp_x = data_struct(i_gr).x ; 
    tmp_y = data_struct(i_gr).y ; 
    tmp_n = length(tmp_x) ; 
    
    % Ensure col vectors (to compute corr after)
    data_struct(i_gr).x = reshape(tmp_x,tmp_n,1) ; 
    data_struct(i_gr).y = reshape(tmp_y,tmp_n,1) ; 
    all_n_pts(i_gr) = tmp_n ; 
    caption_lgd{i_gr} = ['Gr',num2str(i_gr)] ; 
end
n_pts = max(all_n_pts) ; 

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
            case 'c_pts'
                c_pts = Val_arg ; 
            case 'caption_lgd'
                caption_lgd = Val_arg ; 
            case 'x_lab'
                x_lab = Val_arg ; 
            case 'y_lab'
                y_lab = Val_arg ; 
            case 'cbar_lab'
                cbar_lab = Val_arg ; 
            case 'fig_sz'
                fig_sz = Val_arg ; 
            case 'use_colormap'
                use_colormap = Val_arg ; 
            case 'add_lin_fit'
                add_lin_fit = Val_arg ; 
            case 'single_fit'
                single_fit = Val_arg ; 
            case 'taille_tick'
                taille_tick = Val_arg ; 
            case 'taille_stick'
                taille_stick = Val_arg ;     
            case 'axis_equal'
                axis_equal = Val_arg ; 
            case 'add_identity'
                add_identity = Val_arg ; 
            case 'sz_pts'
                sz_pts = Val_arg ; 
            case 'add_lgd'
                add_lgd = Val_arg ; 
            case 'markers'
                markers = Val_arg ; 
            case 'aggregated_corr'
                aggregated_corr = Val_arg ; 
            case 'xl'
                xL = Val_arg ; 
            case 'create_fig'
                create_fig = Val_arg ; 
            case 'show_data'
                show_data = Val_arg ; 
            case 'show_box'
                show_box = Val_arg ; 
            case 'x_box'
                x_box = Val_arg ; 
        end
    end
end

if show_box
    % Assign the data to box positions
    nbox = length(x_box) ; 
    data_box = NaN(sum(all_n_pts), nbox) ; % NaN are ignored in boxplots
    inext = 1 ; 
    for i_gr = 1:n_gr
        tmp_x = data_struct(i_gr).x ; 
        tmp_y = data_struct(i_gr).y ;     

        % assign all curr_x to box positions        
        ind = dsearchn(x_box, delaunayn(x_box), tmp_x) ; 
        % interp1(x_box, 1:nbox, tmp_x) ; % only ok if range of tmp_x is covered in x_box! 
        n_boxi_max = 0 ; 
        for ibox = 1:nbox
            ind_kept = ind==ibox ; 
            
            n_boxi = sum(ind_kept) ; 
            n_boxi_max = max(n_boxi_max, n_boxi) ; 
            
            data_box(inext:(inext+n_boxi-1), ibox) = tmp_y(ind_kept) ;
        end
        inext = inext + n_boxi_max ; 
    end
    data_box(inext:end,:) = [] ; % empty useless NaN values 
end

n_mark = length(markers) ; 
if sum(isnan(c_pts))
    if use_colormap
        c_pts = jet(n_pts) ; %parula ; magma
    else
        c_pts = viridis(n_gr) ; %parula ; magma
        %[ ~, c_pts] = get_some_nice_colors(1) ; 
        % c_pts = c_pts(1:n_gr,:) ;
        % perceptually uniform & colorblind friendly: viridis, inferno, (plasma, magma)
        % diverging: my_PiYG, my_bwr, my_seismic, my_coolwarm
        % best: seismic (white around 0) & coolwarm (no white, from blue to red)
        if ~show_data
            c_pts = zeros(n_gr,3) ; 
        end
    end
end

if create_fig
    fig_h = figure('units','centimeters','outerposition',fig_sz,...
        'Name', fig_name) ;
end

for i_gr = 1:n_gr
    tmp_x = data_struct(i_gr).x ; 
    tmp_y = data_struct(i_gr).y ; 
    tmp_n = length(tmp_x) ; 
    idx_mark = mod(i_gr-1,n_mark) + 1 ; 
    
    if use_colormap
        col_used = c_pts(1:tmp_n,:) ; 
    else
        col_used = c_pts(i_gr,:) ; 
    end
    
    % == * == Plot i_gr th data group
    %sz_pts = 0.01 ;  lwdt = 0.01 ; 
    if show_data
        scatter(tmp_x, tmp_y, sz_pts, col_used, markers{idx_mark}, 'LineWidth', lwdt) ;
    else
        scatter(tmp_x, tmp_y, 1, [1,1,1], markers{idx_mark}, 'LineWidth', lwdt) ;
    end
    hold('on') ;
    tmp_corr = corr(tmp_x, tmp_y, 'rows','complete') ;
    
    tmp_S_corr = corr(tmp_x, tmp_y, 'type','Spearman', 'rows','complete') ; % cfr comm_correlation.txt
    
    % plot(surprise(idx_prev1), distUpdate(idx_prev1), '.', 'Color', colors(1,:), ...
    %    'MarkerSize', 10) ; hold('on');
    
    if ~all_corr_tit
        if aggregated_corr
            if i_gr==1
                full_x = vertcat(data_struct(:).x) ; 
                full_y = vertcat(data_struct(:).y) ; 
                full_corr = corr(full_x(:), full_y(:), 'rows','complete') ;
                full_S_corr = corr(full_x(:), full_y(:), 'type','Spearman', 'rows','complete') ;
                title_str = ['r: ', num2str(round(full_corr,3)), ', S: ',...
                    num2str(round(full_S_corr,3))] ; 
            end
        else
            % Average all the group-correlations. 
            if i_gr==1
                mean_tmp_corr = 0 ;
                mean_tmp_S_corr = 0 ; 
            end
            mean_tmp_corr = mean_tmp_corr + tmp_corr ; 
            mean_tmp_S_corr = mean_tmp_S_corr + tmp_S_corr ;
            
            if i_gr==n_gr
                mean_tmp_corr = mean_tmp_corr/n_gr ; 
                mean_tmp_S_corr = mean_tmp_S_corr/n_gr ; 
                %title_str = ['$\bar{r}$: ', num2str(round(mean_tmp_corr,3)), ', $\bar{S}$: ',...
                %    num2str(round(mean_tmp_S_corr,3))] ; 
                title_str = ['$\bar{r}$: ', num2str(round(mean_tmp_corr,3))] ; 
            end
        end
    else
        if i_gr==1
            title_str = ['r1: ', num2str(round(tmp_corr,3)), ', S1: ',...
                num2str(round(tmp_S_corr,3))] ; 
        else
            title_str = [title_str, '; r',num2str(i_gr),': ', num2str(round(tmp_corr,3)),...
                ', S', num2str(i_gr), ': ', num2str(round(tmp_S_corr,3))] ; 
        end
    end
    
    
end
if show_box
    hb = boxplot(data_box, 'colors', 'k', 'Labels', x_box, 'positions', x_box) ; 
    set(hb,{'linew'},{2}) ; hold on ; box off ; 
    %gca.XAxis.TickLabelInterpreter = 'tex' ;
    % add the means
    plot(x_box, nanmean(data_box), 's', 'color', 'b', 'MarkerSize', 10, 'LineWidth', 2) ;
    if nbox>=2
        gap_tmp = x_box(2) - x_box(1) ; 
    else
        gap_tmp = 1 ; 
    end
    xlim([x_box(1)-gap_tmp/2, x_box(end)+gap_tmp/2]) ; 
end

if ~sum(isnan(xL))
   xlim(xL) ;  
   if axis_equal
       ylim(xL) ;
   end
end

% == * == Legend
if n_gr>1 && add_lgd
    xL = get(gca,'xlim') ; 
    for i_gr = 1:n_gr
        idx_mark = mod(i_gr-1,n_mark) + 1 ; 
        h(i_gr) = scatter(xL(1)-1,0,sz_pts, 'k', markers{idx_mark},'LineWidth', 2) ; hold('on') ; 
    end
    set(gca,'xLim', xL) ; 
    lgd = legend(h, caption_lgd) ;  set(lgd,'FontSize',taille_stick) ; 
end
set(gca,'YGrid', 'off', 'XGrid','off','XMinorGrid','off', ...
    'FontSize', taille_stick) ; 
xlabel(x_lab, 'FontSize', taille_tick,'FontName', 'Arial') ; %'Interpreter','Latex'); 
ylabel(y_lab, 'FontSize', taille_tick,'FontName', 'Arial') ; %'Interpreter','Latex'); 
title(title_str, 'FontSize', taille_tick,'Interpreter','Latex')

if axis_equal
    xL = get(gca,'XLim') ; yL = get(gca,'YLim') ; 
    xL(1) = min(xL(1),yL(1)) ; xL(2) = max(xL(2),yL(2)) ; 
    yL = xL ; 
    xlim(xL) ; ylim(yL) ; 
end

if add_lin_fit
    xL = get(gca,'XLim') ;
    hold on ; 
    yL = get(gca,'YLim') ; 
    xval_poly = linspace(xL(1), xL(2), 30) ; 
    poly_mean = zeros(1,length(xval_poly)) ; 
    if single_fit
        gr_to_test = 1 ; 
    else
        gr_to_test = 1:n_gr ;         
    end
    
    for i_gr=gr_to_test       
        poly_tmp = polyval(data_struct(i_gr).coeffs, xval_poly) ; 
        if use_colormap
            col_used = 'k' ; 
        else
            col_used = c_pts(i_gr,:) ; 
        end
        
        plot(xval_poly, poly_tmp, '-', 'Color',col_used,'Linewidth',1.5) ;% '--'
        hold on ;
        xlim(xL) ;ylim(yL) ; 
        
        poly_mean = poly_mean + poly_tmp./n_gr ;
    end
    if ~single_fit
        % show mean fit on top
        plot(xval_poly, poly_mean, '-', 'Color',[0,0,0],'Linewidth',4) ; hold on ;xlim(xL) ;ylim(yL) ; 
    end
end
if add_identity
    xL = get(gca,'XLim') ; hold on ;
    yL = get(gca,'YLim') ;
    plot(xL, xL, ':','LineWidth',2, 'Color', 0.5*ones(1,3)) ; hold on
    xlim(xL) ;
end

% == * == Colorbar
if use_colormap
    colormap(c_pts) ; 
    step = round(n_pts/5,-1) ; 
    real_ticks = [1, step:step:n_pts] ; 
    ticks_hc = (real_ticks-1)/(n_pts-1) ; % colorbar: from 0 to 1
    names_ticks = num2cell(real_ticks) ; 
    hc = colorbar('eastoutside', 'Ticks',ticks_hc, 'TickLabels',names_ticks,...
        'FontSize',taille_stick) ;    
    hc.Label.String = cbar_lab ;
end



if create_fig && save_fig
    %disp([' === Saving fig ', fn_save])
    my_save_fig(fig_h, fn_save, eps_fig, fig_fig, pdf_fig) ;
end

end

