%This function provides the acquisition function MEI * cos
function [obj, cosine, modified_exp_imp] = determinecosvalue(x, all_models, X_d, alpha_A, low, up, k,...
    rhs, n, t, modified_exp_imp)
    
%Check if input constraints are binding
lower_grad = zeros(k, 1);
upper_grad = zeros(k, 1);
grad_binding_input = [];
for i = 1 : k
    if (x(i) == low(i))
        lower_grad(i) = -1;
        grad_binding_input = [grad_binding_input lower_grad];
    end
    if (x(i) == up(i))
        upper_grad(i) = 1;
        grad_binding_input = [grad_binding_input upper_grad];
    end
end

%Insert your SK predictor here instead of SKpredict_new
y_hat = zeros(t - 1, 1); mspe = zeros(t - 1, 1);
B_predict=1;
for i = 1 : t - 1
    [y_hat(i), mspe(i)] = SKpredict_new(all_models(i + 1, :), x, B_predict);
end
%Find binding output constraints
binding_output = check_binding_output_constraints(y_hat, mspe, alpha_A, rhs, t);

if (sum(binding_output) == 0)%No binding output constraint
    obj = -1; 
    cosine = NaN;
    modified_exp_imp = NaN;
else 
    %Determine goal gradient 
    goal_grad = gradient_estimation(all_models(1, :), x, n, k, X_d);
    %Determine the gradients of the binding constraints
    Big_delta = zeros(k, sum(binding_output));
    column = 1;
    for i = 1 : t - 1
        if (binding_output(i))
            cons_grad = gradient_estimation(all_models(i + 1, :), x, n, k, X_d);
            Big_delta(:, column) = cons_grad;
            column = column + 1;
        end
    end
    %Add binding input constraints if there are any
    if (isempty(grad_binding_input) == 0)
        Big_delta = [Big_delta grad_binding_input];
    end
    %Determine the acquisition function MEI * cos
    cosine = determine_cos(goal_grad, Big_delta, k);
    obj = -modified_exp_imp * cosine;
end



