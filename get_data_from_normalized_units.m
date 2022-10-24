function [xpos_data,ypos_data] = get_data_from_normalized_units(ax_orig,X_norm,Y_norm)
% Convert the xpos_data and ypos_data vectors, expressed in data units, in
% normalized units.

ax_ylim = ylim(ax_orig) ;
ax_xlim = xlim(ax_orig) ;


% Annotation: normalized units so that (0,0) is the lower left corner
% of the FIGURE (not the axis), (1,1) the upper right corner
%get axes position in normalized units relative to figure
oldunits = get(ax_orig, 'Units');set(ax_orig, 'Units', 'normalized');
% position of axis within the figure
axpos = get(ax_orig, 'Position') ;            
set(ax_orig, 'Units', oldunits) ;


% conversion ratio between data and normalized units, based on axis
% position
xdata_per_norm = diff(ax_xlim)/axpos(3) ;
ydata_per_norm = diff(ax_ylim)/axpos(4) ; 


xpos_data = (X_norm - axpos(1)).* xdata_per_norm + ax_xlim(1) ;
ypos_data = (Y_norm - axpos(2)).* ydata_per_norm + ax_ylim(1) ;
end
