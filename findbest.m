%This function determines the best feasible solution so far
%If no point is feasible, then use Monte Carlo to find a feasible solution;
%See Hartmann-6 problem
function kriging_best = findbest(X, all_models, t, alpha_infe, rhs, low, up)

n = size(X, 1);
y_hat = zeros(n, t);  
mspe = zeros(n, t);  
B = 1;
%Insert your SK prediction here.
%E.g., we use SKpredict_new
for i = 1 : n
    for j = 1 : t
        [y_hat(i, j), mspe(i, j)] = SKpredict_new(all_models(j, :), X(i, :), B);
    end
end

%Find the "feasible" points
feasible = check_feasible(y_hat(:, 2 : t), mspe(:, 2 : t), n, alpha_infe, rhs, t);
y_hat_0_min = 10^9;

if (sum(feasible) >= 1)
    for i = 1 : n
        if (feasible(i))
            if (y_hat(i, 1) < y_hat_0_min)
                y_hat_0_min = y_hat(i, 1);
                index = i;
            end
        end
    end
    kriging_best = [X(index, :) y_hat(index, :) sqrt(mspe(index, :))];
else 
    feasible_new = 0;
    y_hat_new = zeros(1, t);  
    mspe_new = zeros(1, t);  
    while (feasible_new == 0)
        new_point = unifrnd(low, up);
        %Insert your SK prediction here.
        %E.g., we use SKpredict_new 
        for j = 2 : t
            [y_hat_new(1, j), mspe_new(1, j)] = SKpredict_new(all_models(j, :), new_point, B);
        end
        feasible_new = check_feasible(y_hat_new(1, 2 : t), mspe_new(1, 2 : t), 1, alpha_infe, rhs, t);
        if (feasible_new)
            %Insert your SK prediction here.
            [y_hat_new(1, 1), mspe_new(1, 1)] = SKpredict_new(all_models(1, :), new_point, B);
            kriging_best = [new_point y_hat_new(1, :) sqrt(mspe_new(1, :))];
        end
    end
end
