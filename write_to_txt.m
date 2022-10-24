function [ ] = write_to_txt(filename, str_to_write,permission,nb_lines)
% Add str_to_write and then nb_lines empty lines in the file filename.
if nargin<3
    permission = 'a' ; 
end
if nargin<4
    nb_lines = 0 ; 
end

file_txt = fopen(filename,permission);
fprintf(file_txt, '%s', str_to_write) ;
fclose(file_txt);
write_empty_lines(filename, nb_lines,'a') ; 
end

