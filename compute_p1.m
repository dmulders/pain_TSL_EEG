function p1 = compute_p1(p_1from2,p_2from1)
% Compute p(I1) for the input values of p(I_1|I_2) and p(I_2|I_1) given in
% arguments.
% p_1from2 and p_2from1 must have the same size. 

p1 = p_1from2./(p_1from2 + p_2from1) ; 

end


