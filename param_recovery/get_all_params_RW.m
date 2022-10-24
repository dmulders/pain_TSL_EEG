function all_params_RW = get_all_params_RW()

% learning rates, in [0, 1]
all_a = 0.05:0.05:1 ; % 20 values; 3*21: 63
all_a = [0.005:0.005:0.04, 0.05:0.05:(1-0.05)] ; % 27 values
all_a = [0.005:0.01:(1-0.05)] ; % 95 values
all_a = [0.005:0.005:0.04, 0.05:0.01:(1-0.05)] ; % 99 values
n_a = length(all_a) ; 
all_learned_param = {'transition', 'frequency', 'alternation'} ; 
n_learned_params = length(all_learned_param) ;

n_params = n_a*n_learned_params ; 


all_params_RW = repmat(struct('learned_param', 'transition',...
    'a_RW', 0.1, 'str_model', '_dTP_a10'), n_params,1) ;  
% Preallocate the array which will contain all the sets of params


idx_param = 1 ; 


for idx_learned_param = 1:n_learned_params
    learned_param = all_learned_param{idx_learned_param} ;      
    for a_RW = all_a
        all_params_RW(idx_param).learned_param = learned_param ;
        all_params_RW(idx_param).a_RW = a_RW ;
        
        str_model = get_model_param_str('', learned_param,0,0, 0, 0, 0, 0, ...
            0, 0, 0,[],'RW',a_RW) ; 
        all_params_RW(idx_param).str_model = str_model ; 

        idx_param = idx_param + 1 ;
    end
end


end