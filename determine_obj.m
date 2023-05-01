%This function provides the acquisition function that is cosine times MEI
function objective = determine_obj(x, all_models, X_d, alpha_A, t, low, up, k, rhs, fmin)

n=size(X_d, 1);

%Compute MEI 
[modified_exp_imp, ~] = compute_modified_exp_imp(all_models(1,:), x, fmin);

%Provide acquisition function objective = MEI * cos
[objective, ~, ~] = determinecosvalue(x, all_models, X_d, alpha_A, low, up, k, rhs, n, t, modified_exp_imp);









