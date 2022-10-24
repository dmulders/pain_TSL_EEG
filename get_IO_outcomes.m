function [surprise, distUpdate, conf_p1, unpred, p1_mean, pred_error, N_used, p1_sd] = ...
    get_IO_outcomes(in, s)
% Return key outputs inferred by the Bayesian IO with parameters specified
% in in. 
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ===== * ===== Optimal inference (IO)
in.s = s ; 
out_SO = IdealObserver(in) ;
p1_mean = out_SO.p1_mean ; %p1_MAP = out_SO.p1_MAP ; 
p1_sd = out_SO.p1_sd ; conf_p1 = -log2(p1_sd) ;
N_used = out_SO.N_used + 2 ; % convert counts to (alpha+beta) in beta distrib

distUpdate = out_SO.distUpdate ; 
surprise = out_SO.surprise ; 
unpred = out_SO.unpred ; 
pred_error = out_SO.pred_error ; % to test predictive coding

%conf_p1     = make_row_vector(conf_p1) ; 


end
