function [mse,rmse] = compute_mse_rmse(x,y)
% Compute mean squared error and root mean squared error between vector x
% and y
% Input
% - x, y : vectors of same length
% Output :
% - mse : mean squared error
% - rmse : root mse

[x,Sx] = make_col_vector(x) ;
[y,Sy] = make_col_vector(y) ;

if max(min(Sx),min(Sy)) > 1
    error('x and y must be vectors of the same size') ;
end

if Sx(1) ~= Sy(1)
    error('In compute_mse_rmse : x and y must have same length')
end

v = x-y ;

mse = (v')*v/Sx(1) ;
rmse = sqrt(mse) ;

end

