function test_my_colormaps( )
% Display the colors of some colormaps.


% -- test a hand define colormap for curves
n_colors = 64 ; % 14 ; 
cm_data = viridis(n_colors) ; %my_blues3(n_colors) ; 
cm_data = parula(n_colors) ;
cm_data = redblue(n_colors) ; 

% magma, viridis, inferno, plasma, my_PiYG

fig_tmp = figure() ; set(gca, 'ColorOrder', cm_data, 'NextPlot', 'replacechildren');
plot(1:50, [1:n_colors]'*[1:50], '-') ; % plot(1:50, randn(14,50) , '-')

% show heatmap of a colormap
n_data = 100 ; 
cmap = cm_data ; %parula ; % jet
fig_tmp = figure('units','normalized','outerposition',[0.2 0.5 0.5 0.2]) ; 
imagesc(1:n_data) ; 
colormap(cmap) ; caxis([1, n_data]) ; 
%Caxis rescales the colormap so anything below caxis(1) is the 
% minimal color and anything above caxis(2) is the top of the colormap
set(gca,'ytick',[]) ; %,'yticklabel', descr_first_sites)
set(gca,'xtick',[]) ; %1:nb_second_sites,'xticklabel', descr_second_sites)


idx_color = 1 ; 
fig_last = figure('units','centimeters','outerposition',[5, 15, 5, 5]) ; 
plot(0,0,'.','color', cm_data(idx_color,:), 'MarkerSize', 100)

end



