function [output] = F(b)
global r C m b_s Sq

% Geometry factor
b_r = [0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0]'; % b/r
f = [0.595, 0.558, 0.554, 0.564, 0.608, 0.695, 0.858]'; % Geometry factor

cur_b_r = b/r;  % Current b/r

output = interp1(b_r,f,cur_b_r); % Current geometry factor

end

