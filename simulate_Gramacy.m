%This function provides simulation replications for Gramacy's problem
function W_all = simulate_Gramacy(X, covariance_matrix, replicates, t)

n = size(X, 1);
e = zeros(n, replicates); 

for i = 1 : t
    for j = 1 : n
        e(j, :) = normrnd(0, sqrt(covariance_matrix(j, i)), 1, replicates);
    end
    if (i == 1)
        W_all((i - 1) * n + 1 : i * n, :) = X(:, 1) + X(:, 2) + e;
    elseif (i == 2)
        W_all((i - 1) * n + 1 : i * n, :) = 1.5 - X(:, 1) - 2 * X(:, 2) - 0.5 * sin(2 * pi * ...
            (X(:, 1).^2 - 2 * X(:, 2))) + e;
    else
        W_all((i - 1) * n + 1 : i * n, :) = X(:, 1).^2 + X(:, 2).^2 -3/2 + e;
    end
end



