%This function determines whether there is a binding output constraint at x
function binding_output = check_binding_output_constraints(y_hat, mspe, alpha_A, rhs, t)

z = norminv(1-(alpha_A / (2 * (t - 1)))); 
binding_output = zeros(t - 1, 1);

for i = 1 : t - 1
    if (abs(y_hat(i) - rhs(i)) / sqrt(mspe(i)) <= z)
        binding_output(i) = 1;
    end
end

        
