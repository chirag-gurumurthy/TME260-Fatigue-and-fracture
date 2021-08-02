% Given design data:
% Load: rotating bending Material:

sigma_u = 800; % [MPa] or 116 ksi
% Surface: machined
% Geometry:
d1 = 38; % [mm]
d2 = 40; % [mm]
rho = 6; % [mm]

% Derivation
% accourding to Juvinall

% m=m(e)m(t)m(d)m(s); data from Table 10.1 and Figure 10.10
m_e = 0.5; % steel, sigma_u <= 1400
m_t = 1; % bending
m_d = 0.9; % 10<d<50 [mm]
m_s = 0.73; % sigma_u = 116 ksi and machined surface
m = m_e * m_t * m_d * m_s;
m_prim = 0.9; % bending
% stress koncentration factor:
% Figure A.12 c)
parameter1=rho/d1;
parameter2=d2/d1;
% usage of parameter1 and parameter2 in A.12 c) gives
k_t=1.55;
% derivation of k_f:
% Peterson gives:
log_alpha = 2.654*1e-7*sigma_u^2-1.309*1e-3*sigma_u+0.01103;
alpha = 10^log_alpha;
q = 1/(1+(alpha/rho));
k_f = 1+q*(k_t-1);

% Postprocessing
% Craetion of a Wholer (SN) Curve
loglog([1 1e3 1e6 1e7],[sigma_u m_prim*sigma_u m*sigma_u m*sigma_u])
hold on
loglog([1 1e3 1e6 1e7],[sigma_u m_prim*sigma_u/k_f m*sigma_u/k_f m*sigma_u/k_f],'--')
legend('SN according to Juvinall','reduced for notch')
xlabel('log(N)')
ylabel('log(S_{a})')
title('Wholer (SN) Curve')