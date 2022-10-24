function [NBgA, NAgB, NAgA, NBgB, s1, NA, NB] = CountEventInMemory(s, param, prior_data)
% Compute transitions in sequence s.
%
% Usage: [NBgA, NAgB, NAgA, NBgB, s1, NA, NB] = CountEventInMemory(s, param)
%   Input:
%           s: sequence of event (1s and 2s)
%           param: cell with paired arguments (name, value)
%               * the default is to use all events without decay
%               * {'Limited', 10}: the memory window size is limited to 10 events
%               * {'Decay', 2} an exponential decay exp(-n/2) is applied 
%                 to event n in the past. Note: the latest even is n=1
%               * [NEWW] {'DecayF', 1}: use the largest weights for the first
%                 observations.
%               * [NEWW] {'LimitedF', 5}: window for the first
%               observations.
%               * can be combined, e.g. {'Limited', 10, 'Decay', 2}
%
% Output:
%           Nigj: transition count from j to i
%           Ni: event count for event type i
%           s1: 1st event of the sequence considered.
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 
 
[windowed, windowedF, MemSpan, MemSpanF, leaky_int, leakyF_int, ...
    MemDecay, MemDecayF, Weight, ~] = get_memParams(param) ; 

factorN1 = prior_data.factorN1 ; 
factorN2 = prior_data.factorN2 ; 

bias_weight = 0 ; 
% 1: initial weighting; first weight = exp(-1/MemDecay) => different for
%     each leak...
% 0: always starts the observation weighting with a weight of 1
% ATTENTION: should be adapted accordingly in ComputeLikelihood1stEvent in
% MarkovNoJump.m

% Conditions to window the sequence forwards and backwards
cond_winB = windowed == 1 && length(s) >= MemSpan ; 
cond_winF = windowedF == 1 && length(s) >= MemSpanF ;
cond_two_win = cond_winB && cond_winF && length(s) >= (MemSpan+MemSpanF) ; 
% Compute the windowed sequence to consider
if cond_two_win
    % get only recent events AND first events
    trn      = [diff(s(1:MemSpanF)), diff(s(end-MemSpan+1:end))] ; % 1: A -> B; -1: B -> A;
    subs4trn = [s(1:MemSpanF-1), s(end-MemSpan+1:end-1)] ;
    subs     = [s(1:MemSpanF), s(end-MemSpan+1:end)] ;
    
    % 1st event
    s1 = s(1);    
elseif cond_winB && ~cond_winF
    % get only recent events
    trn      = diff(s(end-MemSpan+1:end)); % 1: A -> B; -1: B -> A;
    subs4trn = s(end-MemSpan+1:end-1);
    subs     = s(end-MemSpan+1:end);
    
    % 1st event
    s1 = s(end-MemSpan+1);
elseif cond_winF && ~cond_winB
    % get only first events
    trn      = diff(s(1:MemSpanF)) ; % 1: A -> B; -1: B -> A;
    subs4trn = s(1:MemSpanF-1) ;
    subs     = s(1:MemSpanF) ;
    
    % 1st event
    s1 = s(1);
else
    % get all event if they are fewer events than the memory span
    trn      = diff(s(1:end)); % 1: A -> B; -1: B -> A;
    subs4trn = s(1:end-1);
    subs     = s;
    
    % 1st event
    s1 = s(1);
end

if MemDecay==0 && MemDecayF==0
    % Avoid no integration
    %warning('Setting MemDecay and MemDecayF to 1 to avoid no integration') ;
    MemDecay = 1 ; MemDecayF = 1 ; 
end
    
if MemDecay < Inf || MemDecayF < Inf 
    % if both MemDecay are Inf, it is a perfect integration
    
    if nargout <=5 % as in MarkovNoJump.m
        % the caller requests count about transitions
        % for speed, compute only that
        
        % Compute the decay factor (values: from small to large)
        TrnDecay = zeros(1,length(subs4trn)) ;         
        if leaky_int && MemDecay>0
            TrnDecay = TrnDecay + Weight.*exp(-(1/MemDecay)*(length(subs4trn)+bias_weight - ...
                [1:length(subs4trn)]) );
        end        
        if leakyF_int && MemDecayF>0
            TrnDecay = TrnDecay + fliplr(exp(-(1/MemDecayF)*(length(subs4trn)+bias_weight - ...
                [1:length(subs4trn)])) );
        end
        
        % compute observed transition (with decay & memory span)
        NBgA = factorN2*sum( (trn(subs4trn==1)==1)  .* TrnDecay(subs4trn==1));
        NAgB = factorN1*sum( (trn(subs4trn==2)==-1) .* TrnDecay(subs4trn==2));
        NAgA = factorN1*sum( (trn(subs4trn==1)==0)  .* TrnDecay(subs4trn==1));
        NBgB = factorN2*sum( (trn(subs4trn==2)==0)  .* TrnDecay(subs4trn==2));
        
    elseif nargout == 2 % as in BernoulliNoJump.m
        % the caller requests count about stimuli
        % for speed, compute only that
        
        % Compute the decay factor
        SDecay = zeros(1,length(subs)) ; 
        if leaky_int && MemDecay>0
            SDecay = SDecay + Weight.*exp(-(1/MemDecay)*(length(subs) + bias_weight - ...
                [1:length(subs)]) );
        end
        if leakyF_int && MemDecayF>0
            SDecay = SDecay + fliplr(exp(-(1/MemDecayF)*(length(subs) + bias_weight - ...
                [1:length(subs)]) )) ; 
        end
        
        % compute observed event count (with decay & memory span)
        NA = factorN1*sum(SDecay(subs==1));
        NB = factorN2*sum(SDecay(subs==2));
        
    else
        % the call request everything
        
        % Compute the decay factor
        SDecay = zeros(1,length(subs)) ; 
        TrnDecay = zeros(1,length(subs4trn)) ; 
        if leaky_int && MemDecay>0
            SDecay   = SDecay + Weight.*exp(-(1/MemDecay)*(length(subs) +bias_weight - ...
                [1:length(subs)    ]) );
            TrnDecay = TrnDecay + Weight.*exp(-(1/MemDecay)*(length(subs4trn)+bias_weight - ...
                [1:length(subs4trn)]) );
        end
        if leakyF_int && MemDecayF>0
            SDecay = SDecay + fliplr(exp(-(1/MemDecayF)*(length(subs) +bias_weight - ...
                [1:length(subs)]) )) ;
            TrnDecay = TrnDecay + fliplr(exp(-(1/MemDecayF)*(length(subs4trn)+bias_weight - ...
                [1:length(subs4trn)])) );
        end
        
        % compute observed transition (with decay & memory span)
        NBgA = factorN2*sum( (trn(subs4trn==1)==1)  .* TrnDecay(subs4trn==1));
        NAgB = factorN1*sum( (trn(subs4trn==2)==-1) .* TrnDecay(subs4trn==2));
        NAgA = factorN1*sum( (trn(subs4trn==1)==0)  .* TrnDecay(subs4trn==1));
        NBgB = factorN2*sum( (trn(subs4trn==2)==0)  .* TrnDecay(subs4trn==2));
        
        % compute observed event count (with decay & memory span)
        NA = factorN1*sum(SDecay(subs==1));
        NB = factorN2*sum(SDecay(subs==2));
    end
    
else
    % although the previous lines are correct whatever MemDecay, for speed
    % avoid computing decay when there is not (= always 1)
    % compute observed transitions (within memory span)
    NBgA = factorN2*sum( trn(subs4trn==1)==1 );
    NAgB = factorN1*sum( trn(subs4trn==2)==-1);
    NAgA = factorN1*sum( trn(subs4trn==1)==0 );
    NBgB = factorN2*sum( trn(subs4trn==2)==0 );
    
    % compute only if asked
    if nargout == 7 || nargout == 2
        NA = factorN1*sum(subs==1);
        NB = factorN2*sum(subs==2);
    end
end

end
