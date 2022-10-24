function [L,pLL,pLH]=seq_source(trial)
% Define the chunk lengths (geometrically distributed). 
%
% Input
%   trial   nb of trials in the sequence. 
%
% Output
%   L   (iterations x 1) cell with the chunks length in each entry. 
%   pLL     (iterations x 1) vector of transition proba L-L.
%   pLH     (iterations x 1) vector of transition proba L-H.

% chunk = sub-sequence with fixed TPs. 
% The chunk lengths follow a geometric distribution (ie it is the number of
% indep Bernouilli experiments required before a first success) of proba
% p_change (Meyniel2015 p16). 


%%% Parameters %%%%%%
iterations = 500;   % nb of times to define the chunks. 
p_change=0.014;     % 0.016 gives a mean of 61 trials without change; 
                    % 0.013 gives a mean of 76 trials without change
%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) GENERATE CHUNKS (L)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iter = 1:iterations
    % draw from geometric distribution
    
    Ndraw = 1000;
    R = rand(Ndraw,1)>=p_change;
    
    % calculate chunk length and truncate chunks > 80 trials & < 5 trials
    
    change_idx = find(R==0);
    %     CL = change_idx(1);
    for i = 2:length(change_idx)
        CL(i) =  change_idx(i)-change_idx(i-1);
    end
    % CL: chunk lengths
    trunc=find(CL>200);
    CL(trunc)=[];
    trunc=find(CL<5);
    CL(trunc)=[];
    
    if ~isempty(CL)        
        % select chuncks for N trials per session
        
        CL_sum = cumsum(CL);
        stop_idx = find(CL_sum <= trial);
        last_C=trial-CL_sum(length(stop_idx));
        tmp=CL(stop_idx);
        tmp(length(stop_idx)+1)=last_C;
        L{iter}=tmp;
        
        if length(L{iter}) <= 2
            L{iter} = [];
        end
        
        %         N_chuncks(iter) = length(L{iter});
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) GENERATE TRANSITION PROBABILITIES (p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p_base=[0.15:0.05:0.85];

for iter = 1:iterations
    
    p_base=p_base(randperm(length(p_base)));
    pLL(iter)=p_base(1);
    p_base=p_base(randperm(length(p_base)));
    pLH(iter)=p_base(1);
    
end

minchange=0.2;
% discard those p that change < 0.2 from one to the next
jj=1;
for j = 2:iterations
    
    jj=jj+1;
    if abs(pLL(jj)-pLL(jj-1)) < minchange
        pLL(jj)=[];
        jj=jj-1;
    else
    end
    
end

jj=1;
for j = 2:iterations
    
    jj=jj+1;
    if abs(pLH(jj)-pLH(jj-1)) < minchange
        pLH(jj)=[];
        jj=jj-1;
    else
    end
    
end

end
