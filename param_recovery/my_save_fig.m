function my_save_fig(fig_h, fn_save, eps_fig, fig_fig, pdf_fig)
% Save figure fig_h.

if nargin<5, eps_fig = 0 ; end
if nargin<4, fig_fig = 0 ; end

if nargin<3, pdf_fig = 0 ; end

if eps_fig
    file_extension = '.eps' ; 
    format_to_save = 'epsc' ; 
elseif fig_fig
    file_extension = '.fig' ; 
    format_to_save = 'fig' ; 
else
    file_extension = '.png' ; 
    format_to_save = 'png' ; 
end


if pdf_fig
    fig_h.PaperPositionMode = 'auto';
    fig_pos = fig_h.PaperPosition ; 
    fig_h.PaperSize = [fig_pos(3) fig_pos(4)] ; 
    % best but will leave some white space...
    %print('fn_save','-dpdf','-fillpage') % -bestfit
    
    % perfect if no white space
    print(fig_h,fn_save,'-dpdf') % -painters
    
    %print(fig_h,'fn_save','-dpdf', '-bestfit')
else
    disp([' === Saving fig ', fn_save, file_extension])
    saveas(fig_h, [fn_save, file_extension], format_to_save) ;
end


end

