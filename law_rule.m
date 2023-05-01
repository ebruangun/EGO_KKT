%This function applies Law's rule to determine the number of replications
function [all_outputs_new, replicate] = law_rule(infill, variance_new, all_outputs_new, replicate,...
    alpha_L, gamma, t)

stop = 0;
Sigma_e_bar_new = zeros(1, t); mean_output_new = zeros(1, t);
n = size(infill, 1);

while (stop == 0)
    critical_value = tinv(1 - (alpha_L / 2), replicate - 1);
    for i = 1 : t
        [Sigma_e_bar_new(i), mean_output_new(i)] = intrinsic_variances(all_outputs_new(i, :), n, replicate);
    end
    half_intervals = critical_value * sqrt(Sigma_e_bar_new);
    if ~(half_intervals / abs(mmean_output_new) < (gamma / (gamma+1)))
        replicate = replicate + 1;
        %Insert your simulation here
        %E.g., simulation for Gramacy's problem
        output_new = simulate_Gramacy(infill, variance_new, 1, t);
        all_outputs_new(:, replicate) = output_new;
    else
        stop=1;
    end
end