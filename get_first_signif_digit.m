function [f, power10] = get_first_signif_digit(x)
% Return the first significant digit of the input number x. 

power10 = floor(log10(abs(x))) ; 
f = fix(abs(x) .* 10 .^ (-power10)) ;

end