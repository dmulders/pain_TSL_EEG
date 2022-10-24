function [H1_from2,H2_from1] = compute_entropy(p_1from2,p_2from1)
% Compute the entropies of the input TPs p(I_1|I_2) and p(I_2|I_1) given in
% arguments.
% p_1from2 and p_2from1 must have the same size. 

H1_from2 = get_entropy(p_1from2) ; 
H2_from1 = get_entropy(p_2from1) ; 

end

function H = get_entropy(p)
% Returns the entropy of a binary VA with proba p and 1-p. 

H = -p.*log2(p)-(1-p).*log2(1-p) ; 
end

