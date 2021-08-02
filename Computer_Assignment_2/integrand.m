function [CdKn] = integrand(b)
global r C m b_s Sq

% Stress intensity factor
DK = F(b).*Sq.*sqrt(pi.*b/b_s); % DK for growth in meters
CdKn = 1./(C.*DK.^m);
end
