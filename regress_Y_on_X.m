function [beta, mse, Rs, eps] = regress_Y_on_X(Y, X, random_inter, groups)
% Perform multiple regressions using \ and return key metrics.
%   Y = X*beta
% A column of ones is appended to X to add an intercept term. 
%
% Inputs:
%   - Y: (N,M) matrix, where 
%       N = nb of observations and 
%       M = nb of response variables.
%   - X: (N,V) matrix where 
%       V = nb of variables in the model (covariates). 
%
% Outputs:
%   - beta: regression coefficients, without the intercept. (V,M).
%   - mse: mean squared error. (1,M).
%   - Rs: R squared for the regression .(1,M).
%   - eps: residuals (N, M). 
%
if nargin<3
   random_inter = 0 ;  
end

N = size(Y,1) ; 
if random_inter
    % Random intercept (for each group)
    % With full dummy variable coding 
    % https://nl.mathworks.com/help/stats/dummy-indicator-variables.html#btisngt
    rand_groups = unique(groups) ; 
    n_inter = length(rand_groups) ; 
    inter = zeros(N,n_inter) ; 
    for i_gr=1:n_inter
       inter(groups==rand_groups(i_gr),i_gr) = 1 ; 
    end
    X = [inter, X] ; 
else
    % Fixed intercept
    n_inter = 1 ; 
    X = [ones(N,1), X] ; 
end



beta = X\Y ; 
Y_hat = X*beta ; % N_obs x n_data
[mse,~,eps] = compute_mse_rmse_matrix(Y, Y_hat) ; % row vec
% Explained variance = coefficient of determination = 1 - SS_res/SS_tot
% https://en.wikipedia.org/wiki/Coefficient_of_determination
Rs = 1 - mse./ (sum((Y-repmat(mean(Y,1),N,1)).^2,1)./N) ; 

% Remove the intercept term(s)
beta = beta(n_inter+1:end,:) ; 

end

