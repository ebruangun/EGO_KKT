%This function performs the pre-processing step that requires SK parameters estimation  
function [output_all, replications] = pre_process(X, covariance_matrix, replications, output_all, alpha_E, ...
    alpha_L, gamma, t)
    
stop = 0;
n = size(X, 1);
B_predict = 1;
B = ones(n - 1, 1);
z = norminv(1 - (alpha_E / (2 * n * t)));
test = zeros(n, t);

for i = 1 : t
    [Sigma_e_bar((i - 1) * n + 1 : i * n), mean_output((i - 1) * n + 1 : i * n)] = ...
        intrinsic_variances(output_all((i - 1) * n + 1 : i * n, :), n, replications);
end

while stop == 0
    %Apply LOO-CV
    for i = 1 : t
        mean_output_sub = mean_output((i - 1) * n + 1 : i * n)';
        Sigma_e_bar_sub = Sigma_e_bar((i - 1) * n + 1 : i * n)'; 
        for j = 1 : n
            X_j = X;
            X_j(j, :) = [];
            mean_output_j = mean_output_sub;
            mean_output_j(j) = [];
            Sigma_e_bar_j = Sigma_e_bar_sub;
            Sigma_e_bar_j(j) = [];
            x_predict = X(j,:);
            %Insert your own SK parameters estimation and prediction here 
            %E.g., we use SKfit_new and SKpredict_new: See Acknowledgment  
            dmodel_j = SKfit_new(X_j, mean_output_j, B, Sigma_e_bar_j, 2, 3);
            [y_hat_j, mspe_j] = SKpredict_new(dmodel_j, x_predict, B_predict);
            %compute test statistics: See equation (16)
            test(j, i) = (mean_output_sub(j) - y_hat_j) / sqrt(Sigma_e_bar_sub(j) + mspe_j);
        end
    end
    test_statistic= max(abs(test), [], 'all'); %See equation (16)
    if (test_statistic > z)%The LOO-CV test is rejected
        count=0;
        %Apply Law's rule sequentially
        for i = 1 : n
            mean_output_sub = zeros(1, t);
            Sigma_e_bar_sub = zeros(1, t);
            for j = 1 : t
                mean_output_sub(j) = mean_output((j - 1) * n + i); 
                Sigma_e_bar_sub(j) = Sigma_e_bar((j - 1) * n + i);
            end
            critical_value = tinv(1 - (alpha_L / 2), replications(i) - 1);
            half_intervals = critical_value * sqrt(Sigma_e_bar_sub); 
            if ~(half_intervals / abs(mean_output_sub) < (gamma / (gamma+1)))
                replications(i) = replications(i) + 1;
                %Insert your simulation here
                %E.g., simulation for Gramacy's problem
                output_all_new = simulate_Gramacy(X(i,:), covariance_matrix(i,:), 1, t);
                for j = 1 : t
                    output_all((j - 1) * n + i, replications(i)) = output_all_new(j);
                    mean_output((j - 1) * n + i) = mean(nonzeros(output_all((j - 1) * n + i, :)));
                    Sigma_e_bar((j - 1) * n + i) = ((nonzeros(output_all((j - 1) * n + i, :)) - ...
                        mean_output((j - 1) * n + i))' * (nonzeros(output_all((j - 1) * n + i, :)) - ...
                        mean_output((j - 1) * n + i ))) / ((replications(i) - 1) * replications(i));  
                end
            else
                count=count+1;
            end
        end
        if count == n
            %This is the case where one of the GP model is rejected but no
            %additional replication is required by Law's rule, so the KKT
            %algorithm starts with "unvalidated" models
            stop =1;
        end
    else
        stop=1;
        %This is the case when the LOO-CV test fails to be rejected, so the
        %KKT algortihm starts with "validated" models
    end
end

    
