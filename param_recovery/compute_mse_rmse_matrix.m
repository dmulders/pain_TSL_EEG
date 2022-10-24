function [mse,rmse,eps] = compute_mse_rmse_matrix(X,Y)
% Compute mean squared error and root mean squared error between each pair
% of columns in X and Y. 
% 
% Input
% - X, Y : matrices of same size
% Output :
% - mse : mean squared error
% - rmse : root mse
% - eps: residuals (differences btw X and Y). 


eps = X-Y ;
N = size(X,1) ; 

mse = sum(eps.^2,1)./N ;
rmse = sqrt(mse) ;

end

