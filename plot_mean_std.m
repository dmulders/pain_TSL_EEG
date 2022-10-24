
function plot_mean_std(x_data, y_data, y_std, varargin)
% Plot y_data \pm y_std (inputs: row vectors, or each row is plotted 
% separately).
% used in TSL_analyze_ratings.
%
% Inputs
%   x_data: (n_curves, n_pts) matrix of x values
%   y_data: (n_curves, n_pts) matrix

% === Default arguments
eps_fig = 1 ; fig_fig = 0 ; pdf_fig = 0 ; 
fn_save         = 'test' ; 
save_fig        = 1 ; 
create_fig      = 1 ; 
fig_name        = '' ; 
x_lab           = 'x' ; 
y_lab           = 'y' ; 
styles          = {'-'} ; 
add_prior       = 0 ;
neg_pos_prior   = 0 ; 
col_prior       = 0.6*ones(1,3) ; 
lwdt            = 2 ; 
lwdt_mean       = 2.5 ; 
taille_tick     = 20 ;  %
taille_stick    = 16 ;  % small ticks
fig_sz          = [1 2 18 12] ; %[1 2 12 12] ;

col_mean        = [0,0,0] ; 
col_std         = 0.6*ones(1,3) ; 
col_mean        = [0, 128, 255; 102, 204, 0]./255 ; 
col_std         = col_mean ; 
caption_lgd     = {} ; 

shaded_std      = 1 ; % otherwise: dotted
show_max        = 0 ; % indicate max value with a dot

[n_curves,n_data] = size(y_data) ;

logx_curve      = 1 ; 

xlim_val = [min(x_data(:)), max(x_data(:))] ; 
prior_val = 1/(n_data*n_curves) ; 

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
            case 'prior_val'
                prior_val = Val_arg ; 
            case 'taille_tick'
                taille_tick = Val_arg ; 
            case 'taille_stick'
                taille_stick = Val_arg ; 
            case 'add_prior'
                add_prior = Val_arg ; 
            case 'neg_pos_prior'
                neg_pos_prior = Val_arg ; 
            case 'col_prior'
                col_prior = Val_arg ; 
            case 'shaded_std'
                shaded_std = Val_arg ; 
            case 'show_max'
                show_max = Val_arg ; 
            case 'col_mean'
                col_mean = Val_arg ; 
            case 'col_std'
                col_std = Val_arg ; 
            case 'caption_lgd'
                caption_lgd = Val_arg ;
            case 'logx_curve'
                logx_curve = Val_arg ; 
            case 'xlim_val'
                xlim_val = Val_arg ; 
            case 'create_fig'
                create_fig = Val_arg ; 
            case 'styles'
                styles = Val_arg ; 
        end
    end
end
n_cols = size(col_mean,1) ; 
n_cols_std = size(col_std,1) ; 
n_styles = length(styles) ; 

if create_fig
    fig_h = figure('units','centimeters','outerposition',fig_sz,...
        'Name', fig_name) ; hold on ;
end

if add_prior
    x_tmp = xlim_val ;
    y_tmp = prior_val*ones(1,2) ; 
    plot(x_tmp, y_tmp, 'linewidth', lwdt, 'Color',col_prior) ; hold on ; 
    if neg_pos_prior
        plot(x_tmp, -y_tmp, 'linewidth', lwdt, 'Color',col_prior) ; hold on ; 
    end
end

% ====*==== Indicate pm STD
for idx_curve = 1:n_curves
    idx_col = mod(idx_curve-1,n_cols_std) + 1 ; 
    col_std_tmp = col_std(idx_col,:) ; 

    x_tmp = x_data(idx_curve,:) ; 
    y_tmp = y_data(idx_curve,:) ; 
    y_std_tmp = y_std(idx_curve,:) ;

    if shaded_std
        fill([x_tmp, fliplr(x_tmp)], [y_tmp-y_std_tmp, fliplr(y_tmp+y_std_tmp)], col_std_tmp, ...
            'EdgeColor','none', 'FaceAlpha', 0.2) ; 
        hold on ;    
    else
        plot(x_tmp, y_tmp-y_std_tmp, '-.-','color', col_std_tmp, 'Linewidth', lwdt) ; hold on ; 
        plot(x_tmp, y_tmp+y_std_tmp, '-.-','color', col_std_tmp, 'Linewidth', lwdt) ; 
    end
end

% ====*==== Indicate MEAN

for idx_curve = 1:n_curves
    idx_col = mod(idx_curve-1,n_cols) + 1 ; 
    col_tmp = col_mean(idx_col,:) ; 
    
    x_tmp = x_data(idx_curve,:) ; 
    y_tmp = y_data(idx_curve,:) ; 
    is = mod(idx_curve-1, n_styles) + 1 ; 
    h(idx_curve) = plot(x_tmp, y_tmp, styles{is},'color', col_tmp, 'Linewidth', lwdt_mean) ; hold on ; 
end

if show_max
    sz_pt = 49 ; % 25 
    [y_max, max_idx] = max(y_data, [], 2) ; 
    for idx_curve = 1:n_curves
        idx_col = mod(idx_curve-1,n_cols_std) + 1 ; 
        col_std_tmp = col_std(idx_col,:) ; 
        
        scatter(x_data(idx_curve,max_idx(idx_curve)), y_max(idx_curve), sz_pt, 'o', 'LineWidth', 1.5, ...
            'MarkerEdgeColor', [0,0,0], 'MarkerFaceColor', col_std_tmp) ; hold('on') ;
    end
end
  
if logx_curve
   set(gca,'XScale','log') ;  
end

if n_curves>1 && ~isempty(caption_lgd)
    lgd = legend(h, caption_lgd, 'location', 'best') ;  set(lgd,'FontSize',taille_stick) ; 
end

set(gca,'YGrid', 'off', 'XGrid','off','XMinorGrid','off', ...
    'FontSize', taille_stick, 'XLim', xlim_val) ;
xlabel(x_lab, 'FontSize', taille_tick,'FontName','Arial') ; %'Interpreter','Latex'); 
ylabel(y_lab, 'FontSize', taille_tick,'FontName','Arial') ; %'Interpreter','Latex'); 


if save_fig
    my_save_fig(fig_h, fn_save, eps_fig, fig_fig, pdf_fig) ;
end

end

