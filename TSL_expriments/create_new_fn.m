function final_fn = create_new_fn(curr_fn,extension)
% Create non-existing filename with a minimum idx happened at the end of
% the provided str. 
% 
% Input
%   curr_fn     filename, without the idx and without extension!
%   extension   extension, default is '.mat'
%
% Output
%   final_fn    filename with version idx and with the '.mat' extension.
% 

if nargin<2
    extension = '.mat' ; 
end

idx_test = 1 ; 
final_fn = [curr_fn,'_',num2str(idx_test),extension] ; 
% avoid erasing already existing file:
already_exist = exist(final_fn, 'file') ;
while already_exist
    idx_test = idx_test + 1 ; 
    final_fn = [curr_fn,'_',num2str(idx_test),extension] ; 
    already_exist = exist(final_fn, 'file') ;
end

end

