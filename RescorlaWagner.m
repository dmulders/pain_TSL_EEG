function [p1_mean, p1_update, confidence] = RescorlaWagner(s, learned_param, a_LR) 
% Learns the statistics of the input sequence s according to a
% Rescorla-Wagner model (delta rule) with learning rate a_LR. 
%
% Inputs:
%   s     binary sequence of stimuli. Entries are in {1, 2}. 
%   learned_param   in {'transition', 'frequency', 'alternation'}. 
%   a_LR    learning rate, in [0, 1]. 

n_obs = length(s) ; 
p1_mean = 0.5*ones(1, n_obs) ;
p1_update = zeros(1, n_obs) ; 

switch learned_param
    case 'frequency'
        % ======= Learning item frequencies (IF).
        % p1_mean is equal to state value. 
        rewards = s ; 
        rewards(s==2) = 0 ; 
        
        % --*-- 1st step
        V = 0.5 ; % prior for the state value
        update = rewards(1) - V ; 
        p1_update(1) = abs(update) ; 
        p1_mean(1) = V + a_LR*update ; 
        
        % --*-- subsequent steps
        for j=2:n_obs
            update = rewards(j) - p1_mean(j-1) ;          
            p1_update(j) = abs(update) ; 
            p1_mean(j) = p1_mean(j-1) + a_LR*update ; 
        end
        
    case 'alternation'
        % ======= Learning alternation frequencies (AF).
        all_V = 0.5*ones(1,n_obs) ; % P(repetition)
        
        rewards = 1 - [1, abs(diff(s))] ; 
        % when repetition: 1
        % when alternation: 0
                
        % After 1 observation: no alt/rep observed!
        % --*-- Iterate
        for j=2:n_obs
            update = rewards(j) - all_V(j-1) ;          
            p1_update(j) = abs(update) ; 
            all_V(j) = all_V(j-1) + a_LR*update ; 
        end
        
        % --*-- Deduce p1_mean    
        for k = 1:n_obs
            if s(k) == 1
                p1_mean(k) = all_V(k) ;            
            else
                p1_mean(k) = 1-all_V(k) ;
            end
        end
        
    case 'transition'
        % ======= Learning transition probabilities (TP).
        all_V_g1 = 0.5*ones(1,n_obs) ; % state value when previous stimulus was 1
        all_V_g2 = 0.5*ones(1,n_obs) ; % state value when previous stimulus was 2
        
        rewards = s ; 
        rewards(s==2) = 0 ;
                
        % After 1 observation: no transition observed! V stays at 0.5
        % --*-- Iterate
        for j=2:n_obs
            last_s = s(j-1) ; 
            if last_s==1
                % only update all_V_g1
                update = rewards(j) - all_V_g1(j-1) ;          
                p1_update(j) = abs(update) ; 
                all_V_g1(j) = all_V_g1(j-1) + a_LR*update ; 
                % keep last estimate
                all_V_g2(j) = all_V_g2(j-1) ; 
            elseif last_s==2
                % only update all_V_g2
                update = rewards(j) - all_V_g2(j-1) ;          
                p1_update(j) = abs(update) ; 
                all_V_g2(j) = all_V_g2(j-1) + a_LR*update ; 
                % keep last estimate
                all_V_g1(j) = all_V_g1(j-1) ; 
            end
        end
        
        % --*-- Deduce p1_mean    
        for k = 1:n_obs
            if s(k) == 1
                p1_mean(k) = all_V_g1(k) ;            
            else
                p1_mean(k) = all_V_g2(k) ;
            end
        end
        
end

confidence = abs(p1_mean-0.5) ; 
% probability for the choice to be correct = decision confidence

end

