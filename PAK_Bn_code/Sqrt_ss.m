% Square Root of Sum of Squares

function [fx] = Sqrt_ss(x)

% u function should be m x n matrix (n is dimension)
x_dim = size(x,2);
x_square_sum  = zeros(size(x,1),1);
 for i=1:x_dim
     x_square_sum = x_square_sum + x(:,i).^2;
 end
fx = sqrt(x_square_sum);

end % function end