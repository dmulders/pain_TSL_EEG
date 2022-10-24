function short_names = get_short_names(long_names)
% Abbreviations for IF, AF, TP.

short_names = long_names ; 
for i=1:length(long_names)
    curr_n = long_names{i} ; 
    
    if strcmpi(curr_n, 'transition')
        short_names{i} = 'TP' ; 
    elseif strcmpi(curr_n, 'frequency')
        short_names{i} = 'IF' ;
    elseif strcmpi(curr_n, 'alternation')
        short_names{i} = 'AF' ; 
    end

end

end
