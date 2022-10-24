function [x,Sx] = make_row_vector(x)
% Transpose x if size(x,2) == 1
% Input : 
% - x : vector or matrix
% Output : 
% - x : x transposed if size(x,2) is 1
% - Sx : size(x) after x has been transposed if needed

Sx = size(x) ;
if Sx(2) == 1
    x = x' ;
    Sx = size(x) ;
end

end