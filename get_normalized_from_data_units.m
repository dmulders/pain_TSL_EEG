function [X_norm,Y_norm] = get_normalized_from_data_units(ax_orig,xpos_data,ypos_data,axpos)
% Convert the xpos_data and ypos_data vectors, expressed in data units, in
% normalized units.

ax_ylim = ylim(ax_orig) ;
ax_xlim = xlim(ax_orig) ;

if nargin<4
    % Annotation: normalized units so that (0,0) is the lower left corner
    % of the FIGURE (not the axis), (1,1) the upper right corner
    %get axes position in normalized units relative to figure
    oldunits = get(ax_orig, 'Units');set(ax_orig, 'Units', 'normalized');
    % position of axis within the figure
    axpos = get(ax_orig, 'Position') ;            
    set(ax_orig, 'Units', oldunits) ;
end

% conversion ratio between normalized and data units, based on axis
% position
norm_per_xdata = axpos(3) ./ diff(ax_xlim) ;
norm_per_ydata = axpos(4) ./ diff(ax_ylim) ; 

X_norm = (xpos_data - ax_xlim(1)).* norm_per_xdata + axpos(1) ;
Y_norm = (ypos_data - ax_ylim(1)).* norm_per_ydata + axpos(2) ;

end
