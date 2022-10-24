function [in] = get_default_IO_in()
% Return default inputs for the IdealObserver code. 

MemParam        = {'Decay', [8]} ; 
AboutFirst      = 'WithoutFirst' ;
learned_param   = 'transition' ; %'transition' ; % 'frequency', 'transition', 'alternation'



in.learned = learned_param ;        % estimate TPs
in.jump = 0 ;                       % without jump
%in.mode = 'HMM' ;% use the HMM algo (cfr Meyniel2016 p17) (not sampling, Meyniel2015)
in.opt.MemParam = MemParam ;        % memory limit
in.opt.AboutFirst = AboutFirst ;    % use p(y1|theta) = 0.5 instead of a fct 
                                    % of theta to have analytical solutions
n = 40 ;                            % resolution of the univariate probability grid
in.opt.pgrid = linspace(0,1,n) ;    % grid to return full distributions
in.opt.ReturnDist = 0 ;             % Return full posterior distributions
in.opt.ReturnDistUpdate = 0 ; 
eps_prior = 0.1 ; %0.01 ; 
in.opt.priorp1g2 = [1 1] +eps_prior;% uniform Beta prior (when learning TP)
                                    % Add eps to 1 to the prior: avoid NaN
                                    % in MAP when learning TP and AF after
                                    % the first observation (cfr MAP in
                                    % mar_ComputeMAPandPrediction and
                                    % bern_ComputeMAPandPrediction). 
in.opt.priorp2g1 = [1, 1]+eps_prior;% uniform Beta prior (when learning TP)
in.opt.priorp1 = [1,1]+eps_prior ;  % uniform Beta prior (when learning IF, AF)
in.verbose = 1;                     % to check that no default values are used.


end

