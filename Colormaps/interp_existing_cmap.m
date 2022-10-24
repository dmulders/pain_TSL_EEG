function [cmap_interp] = interp_existing_cmap(cmap_orig,m)
% From the colormap cmap_orig (with one color in each row)

hsv = rgb2hsv(cmap_orig); 
cmap_interp=interp1(linspace(0,1,size(cmap_orig,1)),hsv,linspace(0,1,m));
cmap_interp=hsv2rgb(cmap_interp);

end
