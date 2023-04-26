%Step 1 of the KKT algorithm: Enter program parameters
k = input("Enter the dimension of input vector x, 1 <= k <= 20: ");
t = input("Enter the number of simulated outputs (i.e., objective plus constraint(s)): ");
alpha_A = input("Enter the alpha_A value in equation (18): ");
alpha_A_min = input("Enter the alpha_A_min: ");
alpha_E = input("Enter the alpha_E value in equation (15): ");
alpha_infe = input("Enter the alpha_infe value in equation (16): ");
alpha_L = input("Enter the alpha_L value in equation (13): ");
gamma = input("Enter the gamma value in equation (13): ");
epsilon = input("Enter the epsilon value below equation (17): ");
m_0 = input("Enter the initial number of replications: ");

%Enter the lower and upper bounds on x: See l <= x <= u in equation (2)
lower=zeros(k, 1);
upper=zeros(k, 1);
for j = 1 : k
    lower(j) = input("Enter the lower bound on x_j, j = 1, ..., k: ");
    upper(j) = input("Enter the upper bound on x_j, j = 1, ..., k: ");
end

%Enter the right_hand_side(s) of output constraints in equation (1): c_h'
right_hand_sides=zeros(t-1, 1);
for h = 1 : t - 1
    right_hand_sides(h) = input("Enter the right-hand-side of constraint h', h' = 1, ..., t-1: ");
end

%If your simulation is for a synthetic problem (like Gramacy's example,
%Hartmann-6 example, or Angün's example), then insert variances and
%correlations of your intrinsic noises
multivariate_versus_univariate = input("If the intrinsic noises are cross-correlated, enter 1; otherwise", ... 
    "enter 0: ");
intrinsic_noise_variances = zeros(t, 1);
for h = 1 : t
    intrinsic_noise_variances(h) = input("Enter the intrinsic noise variance, h = 1, ..., t: ");
end

%Determine the covariance matrix for the synthetic example
Sigma_e = diag(intrinsic_noise_variances);
if (multivariate_versus_univariate)
    intrinsic_noise_cross_correlations = input("Enter the intrinsic noise cross-correlations: ");
    for i = 1 : t - 1
        for j = i + 1 : t
            Sigma_e(i, j) = sqrt(intrinsic_noise_variances(i) * intrinsic_noise_variances(j)) * intrinsic_noise_cross_correlations;
        end
    end
end
    
%Enter the number of parallel processors to be used 
parallel_process = input("Enter the number of parallel processors to be used: ");

%Determine the number of start points for fmincon
fmincon_start = 10 * k;

%Step 2 of the KKT algorithm: determine the number of input combinations n according to Tao et al. (2020)
if (k <= 6)
    n_initial = (k + 1) * (k + 2) / 2;
else
    n_initial = 5 * k;
end

%Collect best inputs and outputs
best_so_far = [];

%Options for Matlab's fmincon
options=optimoptions('fmincon','SpecifyConstraintGradient', true);

%Use of parallel threads
%Use parallel.pool.Constant to hold a RandStream on each worker
sc = parallel.pool.Constant(RandStream('mrg32k3a'));

%Delete parfor if only one thread is used
parfor restart = 1 : parallel_process
    
    % Get the stream value
    stream = sc.Value;
    % Set the Substream
    set(stream,'Substream',restart);
    % Make this stream the default (for normrnd etc.), and store
    % the old value for later.
    oldGlobalStream = RandStream.setGlobalStream(stream);

    %Set numbers of input combinations, iterations to their initial values
    n = n_initial;
    iteration = 0;
    stopflag = 0;
    
    %Determine the initial LHS design that changes per restart 
    X_initial = lhs_design(n, k, lower, upper); 

    %Obtain initial simulation runs: Insert your own simulation here
    %E.g.: Gramacy's problem with t = 3 outputs
    W_all = simulate_Gramacy(X_initial, Sigma_e, m_0, t);
    
    %Apply pre-process: LOOCV test
    replication_all_inputs = repmat(m_0, n, 1);
    [W_all, replication_all_inputs] = pre_process(X_initial, Sigma_e, replication_all_inputs, W_all, alpha_E, ...
        alpha_L, gamma, t);
    total_budget_used = sum(replication_all_inputs);
    X_design = X_initial;
    
    while stopflag==0 
        
        %Insert your own SK parameters estimation here 
        %E.g., we use SKfit_new: See Acknowledgment 
        B = ones(n, 1);
        all_models = [];
        for h = 1 : t
            [Sigma_e_bar, mean_W] = intrinsic_variances(W_all((h - 1) * n + 1 : h * n, :), n, ...
                replication_all_inputs);
            model = SKfit_new(X_design, mean_W, B, Sigma_e_bar, 2, 3);
            all_models = [all_models; model];
        end
        

        %Step 6 of the KKT algorithm: Among current points, determine the best feasible
        kriging_best = findbest(X_design, all_models, t, alpha_infe, rhs, lower, upper);
        fmin = kriging_best(k+1); 

        %Collect best inputs and their outputs and mspe
        best_so_far = [best_so_far; kriging_best];
        
        flagalpha_A = 0;
        %Steps 7 to 9 of the KKT algorithm: Search for a new infill point
        while (flagalpha_A == 0) 
            %E.g., we use Matlab's fmincon that requires the nonlinear
            %constraints and objective function to be defined as follows
            problem_nonlcons = @(x)non_linear_constraints(x, all_models(2 : t, :), alpha_A, t, k, X_design, rhs);
            problem_obj = @(x)determine_obj(x, dmodel0,dmodel1,dmodel2,X_design,alpha,q,bounds,k,rhs,fmin);
            benchmark = 0;
            %Determine starting points for fmincon
            X_start = lhs_design(fmincon_start, k, lower, upper);
            %Start with an empty infill point
            infill =[];
            
            for prob = 1 : size(X_start, 1) 
               x0 = X_start(prob, :);
               %Insert your numerical optimizer here
               %E.g., we use Matlab's fmincon
               [x, fval, exitflag] = fmincon(problem_obj, x0, [], [], [], [], lower, upper, ...
                   problem_nonlcons, options);
               if (exitflag ~= -2)
                   [modified_exp_imp, ~] = compute_modified_exp_imp(all_models(1, :), x, fmin);
                   if and(modified_exp_imp > epsilon * abs(fmin), fval < benchmark)
                       infill = x; 
                       benchmark = fval;
                   end
               end
            end
            if (isempty(infill) == 0)
                flagalpha_A = 1;
                %If the infill is an old point
                if (size(unique([X_design; infill], 'rows'), 1) < size([X_design; infill], 1))
                    %Halve alpha_A
                    alpha_A = alpha_A / 2;
                else
                    %Step 9 of the KKT algorithm
                    X_design = [X_design; infill];
                    n = n + 1;
                    replicates_new = min(replication_all_inputs);
                    %Insert your own simulation program here
                    %E.g., Gramacy's problem
                    W_all_new = simulate_Gramacy(X, infill, variance_infill, replicates_new, t);
                    %Law's rule to determine the number of replications 
                    [W_all_new, replicates_new] = law_rule(infill, variance_infill, W_all_new,...
                        replicates_new, alpha_L, gamma, t);
                    replication_all_inputs = [replication_all_inputs; replicates_new];
                    column = size(W_all, 2);
                    if (column > replicates_new)
                        for i = 1 : t
                            W_all((i - 1) * n + n, :) = [W_all_new(i, :) zeros(1, column - replicates_new)];
                        end
                    else
                        W_all(:, column + 1 : replicates_new) = zeros(size(W_all, 1), replicates_new - column); 
                        for i = 1 : t
                            W_all((i - 1) * n + n, :) = W_all_new(i, :);
                        end
                    end
                    total_budget_used = total_budget_used + replicates_new;
                end
            else 
                if alpha_A >= alpha_A_min 
                  alpha_A = alpha_A / 2;
                else
                    infill = [];
                    flagalpha_A = 1; 
                    stopflag = 1;
                end
            end
        end
    end
end
delete(gcp);
    



