function [] = plot_single_topo(data_topo, chanlocs, fn_topo, tit_fig, ...
    avg_cmap_levels, eps_fig, fig_fig, pdf_fig, topo_style, numC_topo, cmap_topo, add_hc, hc_title)
% Plot data_topo on a topoplot & save. 

taille_axis = 24 ; 
taille_axis_cbar = 15 ; 
if nargin<13
    hc_title = '' ;
end
if nargin<12
   add_hc = 1 ;  
end

fig_h = figure('units','centimeters','outerposition',[1 2 8 8] ,...
    'Name', tit_fig) ; 

if avg_cmap_levels
    my_topoplot(data_topo, chanlocs, 'style', 'fill', 'numcontour', numC_topo) ; 
else
    my_topoplot(data_topo, chanlocs, 'style', topo_style) ; 
end
colormap(cmap_topo) ;
%caxis([min_cmap, max_cmap]) ; 
axis tight ; 
%title(['I1: ', num2str(t_I1), ' ms'], 'FontSize', taille_axis) ;
if add_hc
    hc = colorbar() ;
    title(hc, hc_title, 'FontSize', taille_axis) ;
    hc.Label.FontSize = taille_axis ; 
    hc.FontSize = taille_axis_cbar ; % for the ticklabels 
end
my_save_fig(fig_h, fn_topo, eps_fig, fig_fig, pdf_fig) ;


end

