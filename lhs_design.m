%This function provides the initial LHS design with midpoints
function X = lhs_design(n, k, low, up)

X_start = lhsdesign(n, k, 'criterion', 'maximin', 'smooth','off');
X = low' + X_start .* (up - low)';


