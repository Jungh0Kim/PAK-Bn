% Design of Experiments by LHS sampling

function [x_doe] = DoE_byLHS(n_dim,N_iniDOE,LHS_bound)

Coded_value = lhsdesign(N_iniDOE,n_dim);
bounds = repmat([-LHS_bound LHS_bound],n_dim,1);
x_doe = zeros(size(Coded_value));
for i = 1:size(Coded_value,2)  % Convert coded values to real-world units
    zmax = max(Coded_value(:,i));
    zmin = min(Coded_value(:,i));
    x_doe(:,i) = interp1([zmin zmax],bounds(i,:),Coded_value(:,i));   
end

end % function end