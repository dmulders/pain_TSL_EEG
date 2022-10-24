function cm_data = my_blues2(m)
% Copied from 
%https://github.com/matplotlib/matplotlib/blob/master/lib/matplotlib/_cm.py
%

cm = [...
0.96862745098039216,  0.98431372549019602,  1.0                ;
0.87058823529411766,  0.92156862745098034,  0.96862745098039216;
0.77647058823529413,  0.85882352941176465,  0.93725490196078431;
0.61960784313725492,  0.792156862745098  ,  0.88235294117647056;
0.41960784313725491,  0.68235294117647061,  0.83921568627450982;
0.25882352941176473,  0.5725490196078431 ,  0.77647058823529413;
0.12941176470588237,  0.44313725490196076,  0.70980392156862748;
0.03137254901960784,  0.31764705882352939,  0.61176470588235299;
0.03137254901960784,  0.18823529411764706,  0.41960784313725491];

if nargin < 1
    cm_data = cm;
else
    hsv=rgb2hsv(cm);
    cm_data=interp1(linspace(0,1,size(cm,1)),hsv,linspace(0,1,m));
    cm_data=hsv2rgb(cm_data);  
end

end
