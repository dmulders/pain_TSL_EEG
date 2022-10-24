function [] = write_latex_table_latencies(row_names, all_data, filename, ...
    col_names, col_header, caption_str, empty_nan, h_mask, h_mask_corr,...
    permission,table_label,varargin)
% write in a txt file the data from all_data so that it forms a table when
% pasted in latex.
%
% In: - row_names: string array containing the name of each row to display.
%     - all_data: matrix, [nb_rows, nb_col]
%     - col_names: string array containing the name of each col to display.
%                  Default: entries are 1:nb_col.
%     - empty_nan: boolean, if true: don't print 'NaN', leave blank space
%     - h_mask: binary matrix indicating which entries in all_data
%               will be indicated in red in the latex table. Default: none
%               of the value is written in red.
%     - col_header: string giving the "label" for the columns 
%                   (for instance 'Parameters').
%
add_std = 0 ;   
all_data_std = NaN ; 
write_integers = 1 ; % otherwise: floats
if write_integers
    fmt = '%d' ; % format of entry
    fmt_long = '%d' ; % %4d
else
    fmt = '%6.3f' ;
    fmt_long = '%4.3f' ; 
end

% Parse varargin
nargs_left = length(varargin) ; 
if nargs_left>0 
    if ~(round(nargs_left/2) == nargs_left/2)
        error('-- There should be an even number of varargin arguments')
    end
    for i = 1:2:nargs_left
        Name_arg = varargin{i};
        Val_arg = varargin{i+1};
        if ~ischar(Name_arg)
            error('Flag arguments must be strings')
        end
        Name_arg = lower(Name_arg);
        switch Name_arg
            case 'add_std'
                add_std = Val_arg ;
            case 'all_data_std'
                all_data_std = Val_arg ; 
        end
    end
end
if sum(size(all_data_std)~=size(all_data)) % one dim is neq
    add_std = 0 ; 
end

if nargin<11; table_label = 'latencies' ; end
if nargin<10; permission = 'w' ; end % could be 'a'
if nargin<9; h_mask_corr = zeros(size(all_data)) ; end
if nargin<8; h_mask = zeros(size(all_data)) ; end

if nargin<7; empty_nan = 0 ; end

if nargin<6
    caption_str = ['Latencies (sec) between the temperature and',...
        ' the rating peaks. '] ; 
end

if nargin<5; col_header = [] ; end

full_str_caption = ['\caption{',caption_str, '} \label{tab:',table_label,'}'] ; 

[nb_rows, nb_col] = size(all_data) ; 
if nargin<4
    col_names = cell(nb_col,1) ;  
    for idx_col=1:nb_col
        col_names{idx_col} = num2str(idx_col) ; 
    end
end

str_col_tab = ['l', repmat(['l'], 1, nb_col)] ; %['l|', repmat(['|l'], 1, nb_col)] ; % l or c
str_header = [col_names{1}, ' '] ; 
for idx_col=2:nb_col
    str_header = [str_header, ' & ', col_names{idx_col}] ; 
end  
file_txt = fopen(filename,permission);
fprintf(file_txt, '%s\n', ['\begin{longtable}[h!]{',str_col_tab,'}']) ;
%fprintf(file_txt, '%s\n', '\hline') ;
if ~isempty(col_header)
    fprintf(file_txt, '%s\n', ['& \multicolumn{',num2str(nb_col),...
        '}{|c}{', col_header,'} \\']) ;
end
fprintf(file_txt, '%s\n', ['& ', str_header, '\\']) ;
fprintf(file_txt, '%s\n', '\hline \hline') ;

for idx_row=1:nb_rows
    curr_row_data = all_data(idx_row,:) ;
    if add_std
        curr_row_std = all_data_std(idx_row,:) ;
    end
    fprintf(file_txt, '%s', [row_names{idx_row}]) ;
    for idx_col=1:nb_col
        curr_data = curr_row_data(idx_col) ;
        if add_std
            curr_std = curr_row_std(idx_col) ;
        end
        if empty_nan && isnan(curr_data)
            fprintf(file_txt, '%s%6', ' & ') ;
        else
            fprintf(file_txt, '%s', ' & ') ;
            
            if h_mask(idx_row,idx_col)
                % blue bullet:
                %fprintf(file_txt, '%s', ' \textcolor{blue}{$\bullet$}') ;
                %fprintf(file_txt, '%6.3f', curr_data) ;
                % or indicate the number in bold!
                if add_std
                    fprintf(file_txt, ['%s',fmt_long,'(',fmt_long,')%s'], ' \textbf{',curr_data,curr_std,'}') ;
                else
                    fprintf(file_txt, ['%s',fmt_long,'%s'], ' \textbf{',curr_data,'}') ;
                end
            else
                if add_std
                    fprintf(file_txt, [fmt,'(',fmt,')'], curr_data,curr_std) ;
                else
                    fprintf(file_txt, fmt, curr_data) ;
                end
            end  
            if h_mask_corr(idx_row,idx_col)
                %fprintf(file_txt, '%s%6.3f%s', ' & \textcolor{red}{', curr_data, '}') ;
                fprintf(file_txt, '%s', ' \textcolor{red}{$\star$} ') ;
            end
            
        end
    end
    fprintf(file_txt, '%s\n', '\\') ;% \hline
end

fprintf(file_txt, '%s\n', full_str_caption) ;
fprintf(file_txt, '%s\n', '\end{longtable}') ;

fclose(file_txt); 
% %6.5f : 6 is the length of the printed string and 5 is the number of
% decimals considered
% %d: integer

end
