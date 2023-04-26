%This function provides the cosine in equation (21)
function cosine = determine_cos(goal_grad, Big_delta, k)

A = Big_delta;
if size(A,1) == size(A,2) %square matrix
    % make it not square: otherwise mldivide does not give the ls
    % estimates
    A=[A zeros(size(A,1),1)];
end
%Determine LS estimation of Lagrange multipliers in equation (4)
nu = A \ (-goal_grad); 

if ((nu >= 0)) %If all Lagrange multipliers are non-negative, compute the cosine
    goal_grad_tilde = zeros(k, 1);
    for i = 1 : size(Big_delta, 2)
        goal_grad_tilde = goal_grad_tilde + nu(i) * Big_delta(:, i);
    end
    %Determine cosine
    dot_product = dot(goal_grad_tilde, -goal_grad); 
    cosine = dot_product / (norm(goal_grad_tilde) * norm(-goal_grad));
else
    cosine = -1; %KKT conditions are not satisfied at current x 
end
