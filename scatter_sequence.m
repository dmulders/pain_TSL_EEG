function [] = scatter_sequence(all_trials, s,varargin)
% Display the input sequence on the current axes. 

% Default plot parameters
mSize           = 16 ;      %
colors          = [0,0,179; 255,67,46]./255 ; 

small_version   = 1 ;       % no labels for the 2 I
cond_name       = '' ;       % ticklabels in small_version
show_ytick_labs = 1 ;       %
taille          = 24 ;      %
taille_tick     = 20 ;      %
taille_stick    = 16 ;      % small ticks
add_xlabel      = 1 ;       % display xlabel

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
            case 'msize'
                mSize = Val_arg ;
            case 'colors'
                colors = Val_arg ;
            case 'show_ytick_labs'
                show_ytick_labs = Val_arg ; 
            case 'taille'
                taille = Val_arg ;
            case 'taille_tick'
                taille_tick = Val_arg ;
            case 'taille_stick'
                taille_stick = Val_arg ;
            case 'add_xlabel'
                add_xlabel = Val_arg ; 
            case 'cond_name'
                cond_name = Val_arg ; 
        end
    end
end

c1 = colors(1,:) ;
c2 = colors(2,:) ; 

idx1 = s==1 ; 
idx2 = s==2 ; 


if small_version
    v1 = 1 ; 
    v2 = 1.6 ; 
    mSize = 18 ; 
    %plot([0,all_trials(end)],[v1,v1],':', 'color', [140,140,140]./255, 'linewidth',1) ; hold on
    %plot([0,all_trials(end)],[v2,v2],':', 'color', [140,140,140]./255, 'linewidth',1) ; hold on
    scatter(all_trials(idx1), s(idx1).*v1, mSize, c1, 'filled', 'linewidth', 1,...
        'MarkerEdgeColor','k', 'LineWidth', 0.1) ; hold on ; 
    scatter(all_trials(idx2), s(idx2)*v2/2, mSize, c2, 'filled', 'linewidth', 1,...
        'MarkerEdgeColor','k', 'LineWidth', 0.1) ; ylim([1,2]) ; hold on ; 
    
    curr_yticks = [v1,v2] ;  
    yticks(1.1*v2) ;
    set(gca,'yticklabels',cond_name,'fontsize',...
        taille_tick,'Ticklength', [0,0.1]) ;
    axis off ;
else
    scatter(all_trials(idx1), s(idx1), mSize, c1, 'filled', 'linewidth', 1,...
        'MarkerEdgeColor','k') ; hold on ; 
    scatter(all_trials(idx2), s(idx2), mSize, c2, 'filled', 'linewidth', 1,...
        'MarkerEdgeColor','k') ; ylim([1,2]) ; hold on ; 
    
    if show_ytick_labs
        curr_yticks = [1,2] ;  
        yticks(curr_yticks) ;
        set(gca,'yticklabels',{'I_1','I_2'},'fontsize',...
            taille_tick,'Ticklength', [0.01,0.05]) ;
    else
        set(gca,'ytick',[])
        set(gca,'FontSize',taille_stick) ; 
        h_lgd = legend('I_1','I_2','location','northeastoutside','orientation',...
            'horizontal') ; %,'fontsize',taille_tick)
        set(h_lgd, 'FontSize', taille_stick)
    end
end
if add_xlabel
    xlabel('Stimuli','FontSize',taille,'Interpreter','Latex') ; 
else
    xlabel(' ','FontSize',taille,'Interpreter','Latex') ; 
    set(gca, 'XTickLabels', '') ; 
end

pos = get(gca,'Position') ; 
h_plot = pos(4) ; 
if add_xlabel
    pos(2) = pos(2) + 0.4*h_plot ; 
end
if show_ytick_labs
    pos(4) = pos(4) - 0.6*h_plot ; 
else
    pos(4) = pos(4) - 0.86*h_plot ; 
end
set(gca,'Position',pos) ;


end

