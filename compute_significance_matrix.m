function [p_values, h_mask, p_values_new_n, h_mask_new_n, conf_intervals, ...
    test_statistic, cohend] = ...
    compute_significance_matrix(all_data, alpha_level, new_n, ...
    paired_tests, mean_H0, always_t_test, two_sided, cell_pairwise_comp, ...
    angles_in_rad)

% Compute all the paired-significances between the COLUMNS of the input
% matrix.
%
% Inputs:
%   paired_tests:   compute all pairwise signficance. Otherwise: consider
%                   only significance wrto a constant (default=1).
%   mean_H0:    if ~paired_tests, gives the mean of the distribution under
%               H_0.
%   cell_pairwise_comp:     cell with one entry for each column of all_data
%                           s.t. the comparisons btw columns i and 
%                           cell_pairwise_comp{i} are performed. 
%   angles_in_rad:  indicates whether the compared quantities are angles
%                   (expressed in radians). If yes and that the tests are
%                   paired, the differences between angles are computed in
%                   [-pi,pi].
%
% Outputs:
%   p_values: matrix of p-values for each column or for each pair of
%               columns, depending on paired_tests.
%   h_mask: 1 when H0 is rejected by the statistical test.
%   conf_intervals: confidence intervals around the mean sample values,
%                   only valid when t-test are performed (and not Wilcoxon 
%                   signed rank). 
%
% Used in TCS2_analyze_ratings and TCS2_analyze_SSEP.
if nargin<9
   angles_in_rad = 0 ;  
end
if nargin<7
   two_sided = 1 ;  
end
% ! only for single sample t-tests
if two_sided
    tail_value = 'both' ; % test whether the mean is different than "M"
else
    tail_value = 'right' ; % test whether the mean is greater than "M"
end
if nargin<6
    always_t_test = 0 ; % never Wilcoxon
end
if nargin<5
   mean_H0 = 0 ;  
end
if nargin<4
    paired_tests = 1 ; 
end

if nargin<3
    new_n = size(all_data,1) + 2 ; 
end

nb_col = size(all_data,2) ; 

if paired_tests
    p_values = NaN*ones(nb_col) ; 
    h_mask = zeros(nb_col) ;  
    %t_stat = NaN*ones(nb_params) ;
    h_mask_new_n = zeros(nb_col) ;
    p_values_new_n = NaN*ones(nb_col) ;
    
    %conf_intervals = NaN*ones(nb_col,nb_col,2) ; 
    conf_intervals = repmat((nanmean(all_data,1))', 1, nb_col,2) - ...
       repmat(nanmean(all_data,1), nb_col,1, 2) ; 
    
    test_statistic = NaN*ones(nb_col) ; 
    cohend = NaN(nb_col) ; 

    for i_param=2:nb_col
    %for i_param=1:nb_params-1    
        row_data = all_data(:,i_param) ;         
        if nargin<8
            curr_j_params = 1:i_param-1 ; 
            % to explore lower triangular part
        else
            curr_j_params = cell_pairwise_comp{i_param} ; 
            curr_j_params = curr_j_params(curr_j_params<=i_param-1) ; 
        end
        for j_param=curr_j_params%1:i_param-1        
        %for j_param=1:i_param-1% to fill lower triangular part
        %for j_param=i_param+1:nb_params (to fill upper triangular part)
            col_data = all_data(:, j_param) ;
            curr_ind = ~isnan(row_data) & ~isnan(col_data) ; 
            % if there is at least data from 2 subjects to compare
            curr_n_subj = sum(curr_ind) ; 
            if curr_n_subj>1
                diff_data = row_data(curr_ind)-col_data(curr_ind) ; 
                meandiff = mean(diff_data) ; 
                stddiff = std(diff_data) ; 
                % https://www.datanovia.com/en/lessons/t-test-effect-size-using-cohens-d-measure/
                cohend(i_param, j_param) = meandiff/stddiff ; 

                if kstest((diff_data-meandiff)./(stddiff)) && ...
                        ~always_t_test
                    % 1: Kolmogorov-Smirnov test rejects the null hypothesis
                    % that the data are (standard) normally distributed
                    disp(['Do a Wilcoxon signed-rank test for parameter sets ', ...
                        num2str(i_param), ' and ', num2str(j_param)])

                    if angles_in_rad
                        [p, h, stats] = signrank(angle_in_interval(diff_data,[-pi,pi]), 0,...
                             'alpha', alpha_level, 'tail', tail_value) ; 
                    else
                        [p, h, stats] = signrank(row_data(curr_ind), ...
                            col_data(curr_ind), 'alpha', alpha_level,'tail', tail_value) ; 
                    end
                    % Wilcoxon signed-rank test: alternative to paired t-test
                    % when population cannot be assumed to be normally 
                    % distributed.
                    curr_test_stat = stats.signedrank ; 
                else
                    % 1: Kolmogorov-Smirnov test does NOT reject the null 
                    % hypothesis that the data follow N(0,1)
                    if angles_in_rad
                        [h,p,ci,stats] = ttest(angle_in_interval(diff_data,[-pi,pi]), 0, ...
                            'Alpha',alpha_level, 'tail', tail_value) ; % [h,p,ci,stats]
                    else
                        [h,p,ci,stats] = ttest(row_data(curr_ind),...
                            col_data(curr_ind),'Alpha',alpha_level, 'tail', tail_value) ; % [h,p,ci,stats]
                        % p = 2*tcdf(abs(stats.tstat),stats.df,'upper')
                    end
                    %t_stat(i_param, j_param) = stats.tstat ;
                    
                    if nargout>2
                        % theoretical p-val obtained if new_n subjects 
                        tstat_new = stats.tstat*sqrt(new_n/curr_n_subj) ; 
                        p_values_new_n(i_param, j_param) = 2*tcdf(abs(tstat_new), ...
                            new_n-1, 'upper') ; 
                        h_mask_new_n(i_param, j_param) = ...
                            p_values_new_n(i_param, j_param)<alpha_level ; 
                    end
                    
                    conf_intervals(i_param,j_param, :) = ci ; 
                    
                    curr_test_stat = stats.tstat ; 
                end
                p_values(i_param, j_param) = p ; 
                h_mask(i_param, j_param) = h ; 
                test_statistic(i_param, j_param) = curr_test_stat ; 
            end
        end
    end
    
else % test against 0
    
    % For each column (=param) of the input matrix, compute significance
    % of the observations wrto 0.
    p_values = NaN*ones(nb_col, 1) ; 
    h_mask = zeros(nb_col, 1) ; 
    
    h_mask_new_n = zeros(nb_col, 1) ;
    p_values_new_n = NaN*ones(nb_col, 1) ;
    
    meand = nanmean(all_data,1) ; stdd = nanstd(all_data, 0, 1) ; 
    cohend = (meand-mean_H0)./stdd ; cohend = cohend' ;  
    conf_intervals = (meand)'*ones(1,2) ; 
    test_statistic = NaN*ones(nb_col, 1) ; 
    
    if always_t_test
        % do all tests at once!
        [h_mask,p_values,conf_intervals,stats] = ttest(all_data,...
            mean_H0,'Alpha',alpha_level, 'tail', tail_value) ;
        test_statistic = (stats.tstat)' ; 
        conf_intervals = conf_intervals' ; 
        h_mask = h_mask' ; 
        p_values = p_values' ; 
    else
        for i_param=1:nb_col
            row_data = all_data(:,i_param) ; 
            curr_ind = ~isnan(row_data) ; 
            % if there is at least data from 2 subjects to compare
            curr_n_subj = sum(curr_ind) ; 
            if curr_n_subj>1
                curr_data = row_data(curr_ind) ;

                if std(curr_data)<10^(-8) || ...
                        kstest((curr_data-mean(curr_data))./(std(curr_data))) &&...
                        ~always_t_test
                    % 1: Kolmogorov-Smirnov test rejects the null hypothesis
                    % that the data are (standard) normally distributed
                    disp(['Do a Wilcoxon signed-rank test for parameter set ', ...
                        num2str(i_param)])

                    if angles_in_rad
                        [p, h, stats] = signrank(angle_in_interval(curr_data,[-pi,pi]), mean_H0,...
                            'alpha', alpha_level, 'tail', tail_value) ;
                    else
                        [p, h, stats] = signrank(curr_data, ...
                            mean_H0, 'alpha', alpha_level,'tail', tail_value) ;
                    end

                    % Wilcoxon signed-rank test: alternative to paired t-test
                    % when population cannot be assumed to be normally 
                    % distributed.
                    curr_test_stat = stats.signedrank ; 
                else
                    % 1: Kolmogorov-Smirnov test does NOT reject the null 
                    % hypothesis that the data follow N(0,1)
                    if angles_in_rad
                        [h,p,ci,stats] = ttest(angle_in_interval(curr_data,[-pi,pi]), mean_H0, ...
                            'Alpha',alpha_level, 'tail', tail_value) ; % [h,p,ci,stats]
                    else
                        [h,p,ci,stats] = ttest(curr_data,...
                            mean_H0,'Alpha',alpha_level, 'tail', tail_value) ;
                    end
                    % p = 2*tcdf(abs(stats.tstat),stats.df,'upper')

                    %t_stat(i_param) = stats.tstat ;

                    if nargout>2
                        % theoretical p-val obtained if new_n subjects 
                        tstat_new = stats.tstat*sqrt(new_n/curr_n_subj) ; 
                        p_values_new_n(i_param) = 2*tcdf(abs(tstat_new), ...
                            new_n-1, 'upper') ; 
                        h_mask_new_n(i_param) = ...
                            p_values_new_n(i_param)<alpha_level ; 
                    end
                    conf_intervals(i_param,:) = ci ; 
                    curr_test_stat = stats.tstat ; 
                end
                p_values(i_param) = p ; 
                h_mask(i_param) = h ; 
                test_statistic(i_param) = curr_test_stat ; 
            end
        end
    
    end
end

end

