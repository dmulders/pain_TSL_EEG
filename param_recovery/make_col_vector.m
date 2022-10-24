function [x,Sx] = make_col_vector(x)
% Transpose x if size(x,1) == 1
% Input : 
% - x : vector or matrix
% Output : 
% - x : x transposed if size(x,1) is 1
% - Sx : size(x) after x has been transposed if needed

Sx = size(x) ;
if Sx(1) == 1
    x = x' ;
    Sx = size(x) ;
end

end
