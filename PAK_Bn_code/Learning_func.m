% Minimize the learning functon of PAK-Bn

function [LF_min,i_min] = Learning_func(mu, sigma,penalty_factor)

% Learning function
U_val = abs(mu)./sigma;
LF_val = U_val.* penalty_factor;
[LF_min, i_min] = min(LF_val);      

end % function end