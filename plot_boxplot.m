
function plot_boxplot(data, varargin)
% Plot boxplots. 
% used in TSL_analyze_ratings.
%
% data
%   (N_data, N_boxplots) data.

N_groups = size(data,2) ; 

% === Default arguments
eps_fig = 1 ; fig_fig = 0 ; pdf_fig = 0 ; 
fn_save         = 'test' ; 
save_fig        = 1 ; 
fig_name        = '' ; 
x_lab           = 'x' ; 
y_lab           = 'y' ; 
title_str       = '' ; 
lwdt            = 2 ; 
taille_tick     = 20 ;  %
taille_stick    = 16 ;  % small ticks
fig_sz          = [1 2 12 12] ; %[1 2 12 12] ;
labels          = 1:N_groups ; 


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
            case 'title_str'
                title_str = Val_arg ; 
            case 'fig_sz'
                fig_sz = Val_arg ; 
            case 'taille_tick'
                taille_tick = Val_arg ; 
            case 'taille_stick'
                taille_stick = Val_arg ; 
            case 'labels'
                labels = Val_arg ; 
        end
    end
end


fig_h = figure('units','centimeters','outerposition',fig_sz,...
    'Name', fig_name) ;
h = boxplot(data, 'colors', 'k', 'Labels', labels) ; 
set(h,{'linew'},{2}) ; hold on ; 
%gca.XAxis.TickLabelInterpreter = 'tex' ;
% add the means
plot(1:N_groups, mean(data), 's', 'color', 'b', 'MarkerSize', 10, 'LineWidth', 2) ; 
 


set(gca,'YGrid', 'off', 'XGrid','off','XMinorGrid','off', ...
    'FontSize', taille_stick, 'FontName', 'Arial') ; 
bp = gca ; box off ; 
bp.XAxis.TickLabelInterpreter = 'tex' ;
xlabel(x_lab, 'FontSize', taille_tick,'FontName', 'Arial') ; %'Interpreter','Latex'); 
ylabel(y_lab, 'FontSize', taille_tick,'FontName', 'Arial') ; %'Interpreter','Latex'); 
title(title_str, 'FontSize', taille_tick,'Interpreter','Latex')



if save_fig
    %disp([' === Saving fig ', fn_save])
    my_save_fig(fig_h, fn_save, eps_fig, fig_fig, pdf_fig) ;
end

end

