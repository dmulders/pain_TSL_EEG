function [] = reduce_hspace_subplot(h_subplots, percent_removed, increase_width)
% Reduce the horizontal space between subplots by percent_emoved percents.

if nargin<3
    increase_width = 1 ; % otherwise, shift towards center
end

pos = get(h_subplots,'position') ; 
n_subplots = length(pos) ; 

right1 = pos{1}(1) + pos{1}(3) ;
left2 = pos{2}(1) ; space_hor = left2-right1 ;
height1 = pos{1}(4) ; 
y1 = pos{1}(2) ; 

if ~increase_width
    all_xshift = zeros(1,n_subplots) ; 
    if mod(n_subplots,2)==0
        max_shift = n_subplots/2 - 0.5 ; 
        all_xshift = max_shift:-1:-max_shift ; 
    else
        max_shift = (n_subplots-1)/2 ; 
        all_xshift = max_shift:-1:-max_shift ; 
    end
    all_xshift = all_xshift.*space_hor.*percent_removed ; 
end

for i_plot=1:n_subplots
    % ==*== Remove some yticklabels
    set(h_subplots(i_plot), 'xlabel', []) ;
    if i_plot>1
        set(h_subplots(i_plot), 'ylabel', []) ; % 'yticklabel',[]
    end
    
    if increase_width
        pos{i_plot}(3) = pos{i_plot}(3) + space_hor.*percent_removed ;
    else
        pos{i_plot}(1) = pos{i_plot}(1) + all_xshift(i_plot) ;
    end

    % ==*== Shift upwards
    pos{i_plot}(2) = y1 + 0.05*height1 ;
    % ==*== Reduce the height
    %pos{i_plot}(4) = height1 - 1.2*delta_h ; 

    set(h_subplots(i_plot),'position',pos{i_plot});
end


end

