%This function provides MEI in equation (18) that uses deterministic mspe   
function [modified_exp_imp, mspe_goal] = compute_modified_exp_imp(dmodel_goal, x, fmin)
    
%Use your own SK predictor and OK predictor for mspe
B_predict = 1;
[y_hat, mspe_goal, mspe_goal_det] = SKpredict_new(dmodel_goal, x, B_predict);
sigma = sqrt(mspe_goal_det);

%Compute MEI
diff = fmin - y_hat;
if (sigma < 10^(-8))
    if (diff >= 0)
        modified_exp_imp = diff;
    else
        modified_exp_imp = 0;
    end
else
    modified_exp_imp = diff * normcdf(diff / sigma, 0, 1) + sigma * normpdf(diff / sigma, 0, 1);
end


