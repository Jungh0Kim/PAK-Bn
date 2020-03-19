% n-ball domain

function [fx] = nB_domain(u, r)

% u is m x n matrix (n is dimension)
u_dim = size(u,2);
u_square  = zeros(size(u,1),1);
 for i=1:u_dim
     u_square = u_square + u(:,i).^2;
 end
fx = u_square - r.^2;

end % function end