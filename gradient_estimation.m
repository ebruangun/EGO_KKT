%This function estimates the gradients when the correlation function is
%specified by the anisotropic Gaussian in equation (11)
function grad_estimates = gradient_estimation(model, x, n, k, X_d)
    
partial_derivatives = zeros(n, k);
tau_squared = model.tausquared;
residuals = model.residuals;
thetas = model.theta;
Sigmahatinv = model.sigmahatinv;
constant_vector = Sigmahatinv * residuals;
minX = min(X_d);  
maxX = max(X_d);
X_d = (X_d - repmat(minX, n, 1)) ./ repmat(maxX - minX, n, 1);
x = [(x(1) - minX)/(maxX - minX) (x(2) - minX)/(maxX - minX)]; 
%See equation (12)
for j = 1 : n
    partial_derivatives(j, 1) = tau_squared * exp(-thetas(1) * (x(1) - X_d(j,1))^2 -...
        thetas(2) * (x(2) - X_d(j,2))^2) * (-2 * thetas(1) * (x(1) - X_d(j,1)));
    partial_derivatives(j, 2) = tau_squared * exp(-thetas(1) * (x(1) - X_d(j,1))^2 -...
        thetas(2) * (x(2) - X_d(j,2))^2) * (-2 * thetas(2) * (x(2) - X_d(j,2)));
end
grad_estimates = partial_derivatives' * constant_vector;


