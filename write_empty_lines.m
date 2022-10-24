function [ ] = write_empty_lines(filename, nb_lines,permission)
% Add nb_lines empty lines in the file filename.
if nargin<3
    permission = 'a' ; 
end

file_txt = fopen(filename,permission);
for idx_line=1:nb_lines
    fprintf(file_txt, '\n') ;
end

fclose(file_txt);
end

