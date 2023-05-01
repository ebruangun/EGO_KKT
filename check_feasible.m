%This function determines which points in X_design are "feasible" 
function feasible = check_feasible(y_hat_cons, mspe_cons, n, alpha_infe, rhs, t)

feasible = zeros(n, 1);
for i = 1 : n
    if (y_hat_cons(i, :) + norminv(1 - (alpha_infe/(t - 1))) * sqrt(mspe_cons(i, :)) <= rhs(:))
        feasible(i) = 1;
    end
end