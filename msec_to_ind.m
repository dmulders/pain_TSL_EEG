function ind_values = msec_to_ind(time_vec, msec_values)
% Convert msec times to indices within the input time vector.

nval = numel(msec_values) ; 
ind_values = zeros(size(msec_values)) ; 

for i_val = 1:nval
    [~,ind_values(i_val)] = min(abs(time_vec-msec_values(i_val)));  
end

end

