function [ output_nb ] = num2str_padd(input_nb, nb_digits)
% Add 0 at the beginning of the created str from input_nb in order to have
% at leat nb_digits digits.
%
% Eg: num2str_padd(5,2) --> '05'

output_nb = num2str(input_nb) ; 
nb_init_digits = length(output_nb) ; 

if nb_init_digits<nb_digits
    output_nb = [repmat('0', 1,nb_digits-nb_init_digits), output_nb] ; 
end

end

