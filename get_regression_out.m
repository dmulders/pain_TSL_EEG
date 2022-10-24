function [beta, B, mse, Rs, Rs_adj] = get_regression_out(lmem)
% Return useful outputs from a linear regressors lmem (fitlme outcome). 

beta = fixedEffects(lmem) ; 
B = randomEffects(lmem) ;

mse = lmem.MSE ; 
EV = lmem.Rsquared ; 
Rs = EV.Ordinary ; 
Rs_adj = EV.Adjusted ; 
% adjusted to accoutn for the fact that R squared increases as we add predictors
% https://en.wikipedia.org/wiki/Coefficient_of_determination
end

