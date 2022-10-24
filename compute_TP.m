function [p_1from2,p_2from1] = compute_TP(p_alt,p_1)
% Compute p(I_1|I_2) and p(I_2|I_1) for the input values of p(alt) and 
% p(I_1)given in arguments.
% p_alt and p_1 must have the same size. 

p_1from2 = p_alt./(2.*(1-p_1)) ; 
p_2from1 = p_alt./(2.*p_1) ; 

end

