function p_alt = compute_p_alt(p_1from2,p_2from1)
% Compute p(alt) for the input values of p(I_1|I_2) and p(I_2|I_1) given in
% arguments.
% p_1from2 and p_2from1 must have the same size. 

p_alt = (2.*p_1from2.*p_2from1)./(p_1from2 + p_2from1) ; 

end