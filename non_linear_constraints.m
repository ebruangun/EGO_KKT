%This function provides the nonlinear constraints with their gradients in
%equation (19)
function  [c, ceq, grad_c, grad_ceq] = non_linear_constraints(x, all_constraint_models, alpha_A, t, k, X_d, rhs)

n = size(X_d, 1);
c = zeros(t-1, 1); 
grad_c = zeros(k, t-1);
ceq = [];
z = norminv(1-(alpha_A / (2 * (t-1))));
B_predict=1;
%Insert your own SK predictor here
for i = 1 : t - 1
    [y_hat, mspe] = SKpredict_new(all_constraint_models(i, :), x, B_predict);
    c(i) = y_hat - z * sqrt(mspe) - rhs(i);
end
 
if nargout>3 %Check if fmincon really asks for the gradients
    for i = 1 : t - 1
        grad_c(:, i) = gradient_estimation(all_constraint_models(i, :), x, n, k, X_d);
    end
end
grad_ceq=[];