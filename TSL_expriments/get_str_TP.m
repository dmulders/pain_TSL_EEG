function str_TP = get_str_TP(p_vec)
% Str describing the TP in inputs.
%
% Input
%   p_vec   length-2 vector of TPs
%

str_TP = ['_TP', num2str(round(100*p_vec(1))), '_', num2str(round(100*p_vec(2)))] ; 


end

