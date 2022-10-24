function cm_data = my_blues3(m)
% Created b y hand to span blue colors, cfr test_my_colormaps

if nargin<1
    m = 100 ;
end

c1 = [102, 255, 255]./255 ;
c2 = [0, 153, 153]./255 ; 
c3 = [51, 153, 255]./255 ; 
c4 = [102,178,255]./255 ; 
c5 = [0, 76,153]./255 ; 
c6 = [0,0,153]./255 ; 
c7 = [51, 51, 255]./255 ; 

cm_data = [c2; c1; c4;c5; c6; c7];

cm_data = interp_existing_cmap(cm_data,m) ; 

end

