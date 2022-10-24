function str_model = get_model_param_str(AboutFirst, learned_param, ...
    decay, decayF, window, windowF, leaky_int, leakyF_int, window_int, ...
    windowF_int, weight, prior_data, model_type,a_LR)

% Parameters describing the model specified by the arguments.
% decayF, windowF, leakyF_int, windowF_int: 
%   F stands for Forwards (or flipped), as opposed to Backwards
%       Forwards: largest weights on first observations of the sequence
%       Backwards: largest weights on last observations of the sequence

% used in TSL_analyze_ratings & TSL_fit_on_ratings.

if nargin<13 ; model_type = 'IO' ; end
if nargin<12 ; prior_data = struct() ; end
if nargin<11 ; weight = 1 ; end 
if nargin<10 ; windowF_int = 0 ; end
if nargin<9 ; window_int = 0 ; end
if nargin<8 ; leakyF_int = 0 ; end
if nargin<7 ; leaky_int = 0 ; end
if nargin<6 ; windowF = Inf ; end
if nargin<5 ; window = Inf ; end
if nargin<4 ; decayF = Inf ; end
if nargin<3 ; decay = Inf ; end
if nargin<2
   error('== Need at least 2 arguments...') 
end

if strcmpi(model_type, 'rw')
    str_model = '_d' ; % delta rule
    if strcmpi(learned_param, 'transition')
        str_model = [str_model, 'TP'] ;     
    elseif strcmpi(learned_param, 'frequency')
        str_model = [str_model, 'IF'] ;   
    elseif strcmpi(learned_param, 'alternation')
        str_model = [str_model, 'AF'] ;  
    else
        error(['Unknown learned_param == ', learned_param,' ==']) ; 
    end
    str_model = [str_model, '_a', num2str(round(1000*a_LR))] ; 
    
    return
end

if strcmpi(learned_param, 'transition')
    str_model = ['_TP'] ;     
elseif strcmpi(learned_param, 'frequency')
    str_model = ['_IF'] ;   
elseif strcmpi(learned_param, 'alternation')
    str_model = ['_AF'] ;  
else
    error(['Unknown learned_param == ', learned_param,' ==']) ; 
end
str_model = [str_model, '_', AboutFirst] ; 


if leaky_int %&& decay<Inf
    str_model = [str_model, '_leak',num2str(decay)] ; 
end
if leakyF_int %&& decayF<Inf
    str_model = [str_model, '_leakF',num2str(decayF)] ; 
end
if leaky_int && leakyF_int && weight~=1 %&& (decay>0 && decayF>0)
    % remove last condition for the case of leaks = (0,0) which is latter
    % replaced by (1,1) (cfr CountEventInMemory)
    % indicate weight, when 
    % * 2 leaks exist 
    % * it is diff from 1
    str_model = [str_model, '_A', num2str(round(10*weight))] ; 
    % A: to avoid confounding with the decays w
end

if window_int && window<Inf
    str_model = [str_model,'_win',num2str(window)] ; 
end
if windowF_int && windowF<Inf
    str_model = [str_model,'_winF',num2str(windowF)] ; 
end

if isfield(prior_data,'priorN1')
    if prior_data.priorN1~=1
        str_model = [str_model, '_p1', num2str(prior_data.priorN1)] ; 
    end
end
if isfield(prior_data,'priorN2')
    if prior_data.priorN2~=1
        str_model = [str_model, '_p2', num2str(prior_data.priorN2)] ; 
    end
end
if isfield(prior_data,'factorN1')
    if prior_data.factorN1~=1
        str_model = [str_model, '_f1', num2str(round(100*prior_data.factorN1))] ; 
    end
end
if isfield(prior_data,'factorN2')
    if prior_data.factorN2~=1
        str_model = [str_model, '_f2', num2str(round(100*prior_data.factorN2))] ; 
    end
end

            
end