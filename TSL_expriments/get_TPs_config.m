function [all_pts, str_fn,labels] = get_TPs_config(p_small, plus_sessions)
% Gives the TPs for some configurations.
% To indicate sessions on the TP matrix
% 
% Inputs
%   p_small         smallest TP considered
%   plus_sessions   in {0,1,2} with 
%                       0: sessions on a cross, 
%                       1: sessions on a plus,
%                       2: bias selectively IF and AF
%
% Outputs
%   labels          labels for each TP. 

p_large = 1-p_small ; 
% str_fn:  to add to filename
switch plus_sessions
    case 1
        str_fn = '_plus' ; 
        all_pts = [0.5,0.5;...
            0.5, p_large;...
            0.5, p_small;...
            p_small,0.5;...
            p_large,0.5] ; % defining the sessions
    case 0
        str_fn = '_cross' ; 
        all_pts = [0.5,0.5;...
            p_small, p_large;...
            p_large, p_small;...
            p_small, p_small;...
            p_large,p_large] ;
    case 2
        str_fn = '_levels' ; % on the level curves of IF and AF
        % --- Bias alterantion frequency:
        [p_1from2_s,p_2from1_s] = compute_TP(0.5,p_small) ; 
        [p_1from2_l,p_2from1_l] = compute_TP(0.5,p_large) ; 
        
        all_pts = [0.5,0.5;...
            % level of p(alt)=0.5
            p_1from2_s, p_2from1_s;...
            p_1from2_l, p_2from1_l;...            
            % level of p(I_1) = p(I_2) = 0.5
            p_small, p_small;...
            p_large,p_large] ;
end
str_fn = [str_fn,num2str(round(10*p_small))] ; 

n_pts = size(all_pts,1) ; 
labels = cell(1,n_pts) ;
for idx_pt = 1:n_pts
    curr_p = round(all_pts(idx_pt,:),2) ; 
    labels{idx_pt} = ['(',num2str(curr_p(1)),', ',num2str(curr_p(2)),')'] ; 
end

end

