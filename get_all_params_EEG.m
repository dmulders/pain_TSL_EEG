function all_params = get_all_params_EEG(consider_flip, combine_decay, ...
    relative_weight, without_AF)

if nargin<4 ; without_AF = 0 ; end

all_AboutFirst  = {'WithoutFirst'} ; 
% since our sequences have length 100, decay should not be way larger than
% 100...
%all_decays      =  [1:20] ; 
% Maheu2019: 55 values tested (cfr fig 5E & p18)

%all_decays      =  [1:20, 22:2:200] ;
all_weights = [1] ; 

if combine_decay
    % include 0 to also have the cases of only F or B integration!
    if relative_weight
        all_decays  = [0:20, Inf] ; % 22 values
        all_weights = [0.5, 1, 2, 3] ; % weights for the bwd leaks (compared to fwd leaks)
    else
        all_decays      = [0:20, 22:2:40, 45:5:100, Inf] ; % 44 values
        all_decays      = [0:20, 22:2:40, 45:5:90, Inf] ; % 42 values
        % 42*42*3 = 5292
    end
elseif consider_flip
    all_decays      = [1:20, 22:2:40, 45:5:150, Inf] ; % 53 values
else
    %all_decays      = [1:20, 22:2:40, 45:5:150, Inf] ; % 53 values, 53*3 = 159
    %all_decays      = [1:20, 25:5:150, Inf] ; % 47 values,  47*3 = 141 ;
    all_decays      = [1:20, 25:5:120, Inf] ; % 41 values, 41*3 = 123 ;
    
    %%%%%%% TMP: to test codes!!!
    %all_decays = [1:20] ; 
end

all_windows     = [] ;
if without_AF
    all_learned_param = {'transition', 'frequency'} ; 
else
    all_learned_param = {'transition', 'frequency', 'alternation'} ; 
% 'frequency', 'transition', 'alternation' (cfr [NEWW] in IdealObserver)
end

if consider_flip
    all_flips = [0,1] ; 
else
    all_flips = 0 ; 
end

n_flips = length(all_flips) ; 
n_first = length(all_AboutFirst) ; 
n_decays = length(all_decays) ; 
n_weights = length(all_weights) ; 
n_windows = length(all_windows) ; 
n_learned_params = length(all_learned_param) ; 

if combine_decay
    all_windows     = [] ; %all_windowsF = [] ; 
    all_flips       = 0 ;
    
    all_decaysF     = all_decays ;
    n_params = n_first*n_learned_params*n_weights*(n_decays^2) ; 
else
    all_decaysF = [0] ; 
    n_params = n_first*n_learned_params*(n_flips*(n_windows + n_decays)) ; 
end

n_flips = length(all_flips) ; 

all_params = repmat(struct('AboutFirst','WithoutFirst', 'learned_param', 'transition',...
    'decay',[Inf],'decayF', [Inf], 'window', [Inf], 'windowF', [Inf], ...
    'leaky_int', 0, 'leakyF_int', 0, 'window_int',0, 'windowF_int', 0, ...
    'MemParam', [], 'weight', 1, 'str_model', '_TP_WithoutFirst'), n_params,1) ;  
% Preallocate the array which will contain all the sets of params


idx_param = 1 ; 

for idx_first = 1:n_first    
    AboutFirst = all_AboutFirst{idx_first} ; 
    for idx_learned_param = 1:n_learned_params
        learned_param = all_learned_param{idx_learned_param} ;      
        for weight = all_weights
            for idx_flip=1:n_flips
                curr_flip = all_flips(idx_flip) ; 

                % ===> Leaky integration
                for decay = all_decays
                    for decayF=all_decaysF   
                        all_params(idx_param).AboutFirst = AboutFirst ;
                        all_params(idx_param).learned_param = learned_param ;
                        all_params(idx_param).weight = weight ;

                        MemParam = {} ;  

                        if combine_decay
                            dB = decay ; dF = decayF ; 
                            leaky_bool = 1 ; leakyF_bool = 1 ; 
                        else
                            % can have all_flips = [0,1] ; when ~combine_decay
                            if curr_flip
                                % use decay (NOT decayF) for the flipped decay!!
                                dB = 0 ; dF = decay ; 
                                leaky_bool = 0 ; leakyF_bool = 1 ;                            
                            else
                                dB = decay ; dF = 0 ; 
                                leaky_bool = 1 ; leakyF_bool = 0 ;
                            end                        
                        end

                        all_params(idx_param).decay = dB ;
                        all_params(idx_param).decayF = dF ; 
                        all_params(idx_param).leaky_int = leaky_bool ;
                        all_params(idx_param).leakyF_int = leakyF_bool ;                      
                        MemParam(end+1:end+2) = {'Decay', dB} ;
                        MemParam(end+1:end+2) = {'DecayF', dF} ;
                        if weight~=1 ; MemParam(end+1:end+2) = {'Weight',weight} ; end
                        all_params(idx_param).MemParam = MemParam ;

                        str_model = get_model_param_str(AboutFirst, learned_param, ...
                            dB, dF, Inf, Inf, leaky_bool, leakyF_bool, 0, 0, weight) ; 
                        all_params(idx_param).str_model = str_model ; 

                        idx_param = idx_param + 1 ; 
                    end
                end

                % ===> Windowed integration
                for window = all_windows
                    all_params(idx_param).window = window ;
                    all_params(idx_param).AboutFirst = AboutFirst ;
                    all_params(idx_param).learned_param = learned_param ; 
                    all_params(idx_param).window_int = 1 ;
                    %all_params(idx_param).weight = weight ; % combined
                    %windows not used

                    MemParam = {} ;
                    if curr_flip
                        % (AboutFirst, learned_param, decay, decayF, window, 
                        %  windowF, leaky_int, leakyF_int, window_int, windowF_int)
                        str_model = get_model_param_str(AboutFirst, learned_param, ...
                            Inf, Inf, Inf, window, 0, 0, 0, 1) ; 

                        MemParam(end+1:end+2) = {'LimitedF', window} ; 
                    else
                        str_model = get_model_param_str(AboutFirst, learned_param, ...
                            Inf, Inf, window, Inf, 0, 0, 1, 0) ; 
                        MemParam(end+1:end+2) = {'Limited', window} ; 
                    end  
                    %if weight~=1 ; MemParam(end+1:end+2) = {'Weight',weight} ; end
                    all_params(idx_param).str_model = str_model ; 
                    all_params(idx_param).MemParam = MemParam ;

                    idx_param = idx_param + 1 ; 
                end

            end
        end
    end
end


end
