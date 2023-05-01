%This function provides the estimated covariance matrix: see equation (9) 
function [Sigma_e_bar, mean_output] = intrinsic_variances(output, n, replications)
   
mean_output = zeros(n, 1);
Sigma_e_bar = zeros(n, 1);

for i = 1 : n
    mean_output(i) = mean(nonzeros(output(i, :)));
    Sigma_e_bar(i) = ((nonzeros(output(i, :)) - mean_output(i))' * (nonzeros(output(i, :)) - mean_output(i))) / ...
        ((replications(i) - 1) * replications(i));
end


        