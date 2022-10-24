function [windowed, windowedF, MemSpan, MemSpanF, leaky_int, leakyF_int, ...
MemDecay, MemDecayF, Weight, Gaussian] = get_memParams(param)
% Extract the memory params from the input argument.
%
% From CountEventInMemory.m. 
% Used in MarkovNoJump.m, CountEventInMemory.m
% Dounia Mulders -- 2019. 

windowed = 0; windowedF = 0 ; MemSpan = Inf ; MemSpanF = Inf ; 
leaky_int = 0 ; MemDecay = Inf; leakyF_int = 0 ; MemDecayF = Inf ; 
Gaussian = 0 ; Weight = 1 ; 

if exist('param', 'var')
    nb_param = length(param) ; 
    if iscell(param) && any(nb_param == [2 4 6 8 10 12])
        for k = 1:2:nb_param % take only odd indices
            if strcmp(param{k}, 'Limited')
                windowed = 1;
                MemSpan = param{k+1};
            end    
            if strcmp(param{k}, 'LimitedF')
                windowedF = 1;
                MemSpanF = param{k+1};
            end  
            
            if strcmp(param{k}, 'Decay')
                leaky_int = 1 ; 
                MemDecay = param{k+1};
            end
            if strcmp(param{k}, 'DecayF')
                leakyF_int = 1 ; 
                MemDecayF = param{k+1};
            end
            if strcmp(param{k}, 'Weight')
                Weight = param{k+1} ; 
            end
            
            if strcmp(param{k}, 'Gaussian')
               Gaussian = param{k+1} ;  
            end
        end
    elseif ~isempty(param)
        error('check parameters')
    end
end

if ~(leaky_int && leakyF_int)
    Weight = 1 ; 
end

end

