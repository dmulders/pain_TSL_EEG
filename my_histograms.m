function my_histograms(data_struct, varargin)
% Plot histograms of input data.
% 
% used in TSL_fit_on_ratings. 
%
% data_struct
%   Structure array with one entry per group of data, which will be
%   associated with a different histogram. Required fields: 
%       'x' for the data set.

n_gr = length(data_struct) ; 

% === Default arguments
eps_fig = 1 ; fig_fig = 0 ; pdf_fig = 0 ; 
fn_save         = 'test' ; 
save_fig        = 1 ; 
fig_name        = '' ; 
x_lab           = 'x' ; 
y_lab           = 'Count' ; 
lwdt            = 2 ; 
taille_tick     = 20 ;  %
taille_stick    = 16 ;  % small ticks
fig_sz          = [1 2 18 12] ; %[1 2 12 12] ;
add_gauss_fit   = 0 ;
axis_equal      = 0 ; 
facealpha       = 0.5 ; 
%[ ~, c_pts] = get_some_nice_colors(1) ; 
% c_pts = c_pts(1:n_gr,:) ;
if n_gr==3
    c_pts       = [102, 204, 0; 0, 128, 255; 255, 0, 127]./255 ; 
elseif n_gr==6
    c_pts       = [102, 204, 0; 0, 128, 255; 255, 0, 127; ...
        0, 0, 204; 0, 153, 76; 153, 51, 255]./255 ; 
else
    c_pts       = parula(n_gr) ; % jet(n_gr) ; %parula ; magma
end


all_ranges = zeros(n_gr,2) ; 
caption_lgd = cell(1, n_gr) ; 

for i_gr = 1:n_gr
    tmp_data = data_struct(i_gr).x ; 
    tmp_n = length(tmp_data) ; 
    
    % Ensure col vectors
    data_struct(i_gr).x = reshape(tmp_data,tmp_n,1) ; 
    all_ranges(i_gr,1) = min(tmp_data) ; all_ranges(i_gr,2) = max(tmp_data) ;
    caption_lgd{i_gr} = ['Gr',num2str(i_gr)] ; 
end
full_range = [min(all_ranges(:,1)), max(all_ranges(:,2))] ; 

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
            case 'add_gauss_fit'
                add_gauss_fit = Val_arg ; 
            case 'taille_tick'
                taille_tick = Val_arg ; 
            case 'taille_stick'
                taille_stick = Val_arg ;     
            case 'axis_equal'
                axis_equal = Val_arg ; 
            case 'facealpha'
                facealpha = Val_arg ; 
        end
    end
end


fig_h = figure('units','centimeters','outerposition',fig_sz,...
    'Name', fig_name) ;

% == * == Histograms
for i_gr = 1:n_gr
    tmp_data = data_struct(i_gr).x ; 
    col_used = c_pts(i_gr,:) ; 
    
    % == * == Plot i_gr th histogram
    h(i_gr) = histogram(tmp_data, 'BinMethod','integers', 'DisplayStyle','stairs') ; hold on;
    % DisplayStyle: 'bar', 'stairs'
    % BinMethod: fd, integers
    set(h(i_gr),'EdgeColor',col_used,'Normalization', 'count') ; % pdf 
    %set(h(i_gr),'FaceColor',col_used,'EdgeColor','w','facealpha', ...
    %    facealpha, 'Normalization', 'count') ; % pdf 
    
end

% == * == Legend
if n_gr>1
    lgd = legend(h, caption_lgd) ; set(lgd,'FontSize',taille_stick) ; 
end

% == * == Gaussian fit
if add_gauss_fit
    x_tmp = full_range(1):0.001:full_range(2) ;
    for i_gr = 1:n_gr
        tmp_data = data_struct(i_gr).x ; 
        col_used = c_pts(i_gr,:) ; 
        [mean_g, std_g] = normfit(tmp_data) ;
        plot(x_tmp, normpdf(x_tmp, mean_g, std_g), '--', 'Color', ...
            col_used, 'Linewidth', lwdt) ; hold on
    end
      
end

% == * == Axes
set(gca,'YGrid', 'on', 'XGrid','on','XMinorGrid','off', ...
    'FontSize', taille_stick) ; 
xlabel(x_lab, 'FontSize', taille_tick,'Interpreter','Latex'); 
ylabel(y_lab, 'FontSize', taille_tick,'Interpreter','Latex'); 

if axis_equal
    xL = get(gca,'XLim') ; yL = get(gca,'YLim') ; 
    xL(1) = min(xL(1),yL(1)) ; xL(2) = max(xL(2),yL(2)) ; 
    yL = xL ; 
    xlim(xL) ; ylim(yL) ; 
end

% == * == Saving
if save_fig
    my_save_fig(fig_h, fn_save, eps_fig, fig_fig, pdf_fig) ;
end



end

