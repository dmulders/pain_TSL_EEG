function str_out = get_str_param(struct_param)
% Based on input parameters (as given by function get_all_params_to_test),
% return a string describing the model parameters. 
%
% struct_param = struct('AboutFirst','WithoutFirst', 'learned_param', 'transition',...
%    'decay',[Inf],'decayF', [Inf], 'window', [Inf], 'windowF', [Inf], ...
%    'leaky_int', 0, 'leakyF_int', 0, 'window_int',0, 'windowF_int', 0, ...
%    'MemParam', [], 'weight', 1, 'str_model', '_TP_WithoutFirst') ; 
     
str_leak = '' ; 
if struct_param.leaky_int
    str_leak = ['; (bwd) leaky integration with leakB = ', num2str(struct_param.decay)] ; 
end
str_leakF = '' ; 
if struct_param.leakyF_int
    str_leakF = ['; fwd leaky integration with leakF = ', num2str(struct_param.decayF)] ; 
end

str_win = '' ; 
if struct_param.windowF_int
    str_win = ['; (bwd) windowed integration with winB = ', num2str(struct_param.window)] ; 
end
str_winF = '' ; 
if struct_param.windowF_int
    str_winF = ['; fwd windowed integration with winF = ', num2str(struct_param.windowF)] ; 
end
str_weight = '' ; 
if struct_param.weight~=1
    str_weight = ['; weighting bwd vs. fwd = ', num2str(struct_param.weight)] ; 
end

str_out = ['Learned parameter: ', struct_param.learned_param, str_win, str_winF, ...
    str_leak, str_leakF, str_weight] ;

end