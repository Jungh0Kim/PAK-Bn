% Limit-state function

function [gx] = G_function(u1, u2)

% 4 component function
k = 7;
g1x = 3+(u1-u2).^2/10-(u1+u2)/sqrt(2);
g2x = 3+(u1-u2).^2/10+(u1+u2)/sqrt(2);
g3x = u1-u2+k/sqrt(2);
g4x = u2-u1+k/sqrt(2);

gx = min(g1x, min(g2x, min(g3x, g4x)));

end % function end