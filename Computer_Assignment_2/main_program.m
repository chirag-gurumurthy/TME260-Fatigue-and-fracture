clear all
close all
clc
global r C m b_s Sq

% INPUT DATA
% Material data

C = 5e-12; % for da/dN in m/cycle, DK in MPa(m)^(1/2)
m = 3; % Crack growth exponent
E = 210e3; % youngs modulus [MPa]
Kic = 60; % Fracture toughness [MPa(m)^(1/2)]

% Geometry
r = 0.08; % Axle radius [m]
b_s = 0.4; % Ratio b/s

% Crack geometry
b_ini = 0.001; % Initial crack depth [m]

% A.1 TASK 1a
% Constant amplitude

M = [60].*1e3; % Constant bending moment

% Variable amplitude (in operations this could be a measured series of values)
prompt1 = sprintf(' If task 1a,b or 2 WRITE 1, if task 1c WRITE 2: ');
prompt2 = sprintf(' If solving task 1a,b or 2 WRITE 0. Do you wish to \n sole 1c using consistant'); 

task = input(prompt1);
disp('----------------------------------------------------')
version = input(prompt2);

if version==1
    %Fist version
    M = [40, 60, 60, 75].*10^3;
elseif version==2
%Second version
    M1 = linspace(40,40,1e5);
    M2 = linspace(60,60,2e5);
    M3 = linspace(75,75,1e5);
    M = [M1,M2,M3].*10^3;
else
    % Constant bending moment
    M = [60].*1e3;
end

% Nominal stress
S = 4*M/(pi*r^3);

% Critical crack depth
Ki=0; % initiate loop
b1a=b_ini; % initial guess for critical crack depth
while abs(Ki/Kic)<1 % failure when Ki>=Kic
    b1a=b1a+b_ini; % update crack depth
    [F1] = F(b1a); % compute new design factor
    Ki = S*F1*sqrt(pi*b1a/b_s)*1e-6; % compute stress intensity factor
end

% Postprocessing task 1a
disp('----------------------------------------------------')
disp('Answer to 1a) is:')
disp('The critical crack depth is in millimeters equal to:')
disp(b1a*1000)
disp('----------------------------------------------------')

% A.2 TASK 1 b,c
% Equivalent bending stress

S = 4.*M/(pi*r^3)*1e-6;
Smax = max(S);
Smin = min(S);

Sq = 0;
for i = 1:size(S,2)
    Sq = Sq + S(i)^m;
end
Sq = (Sq./size(S,2)).^(1/m);

% CRITICAL CRACK SIZE
est_bc = [1e-6 r]; % Estimated interval of critical crack depth
bc = fzero(@(b) Kic-Smax*F(b)*sqrt(pi*b/0.4), est_bc);

% INTEGRATE CRACK GROWTH LIFE
% Evaluate crack growth through integration
bvec=[b_ini:0.001:bc]; % Vector with final crack sizes
noCycles = zeros(size(bvec)); % no of cycles to the (final) crack sizes
for i = 1:size(bvec,2)
    Q = quad(@integrand,b_ini,bvec(i)); % integrate eq 11.34
    noCycles(i) = Q;
    K(i) = F(bvec(i)).*Sq.*sqrt(pi.*bvec(i)/b_s);
end
disp('==== Integration ====')
disp('Number of cycles to failure according to integration:')
Life_int = noCycles(size(bvec,2));
disp(Life_int)
disp('Critical crack depth [m]:')
disp(bc)
disp('Equivalent stress intensity [MPa(m)^(1/2)]:')
disp(K(size(bvec,2)))
disp('----------------------------------------------------')

% Save plot variables
b_int = bvec; % Crack depths in mm
K_int = K; % Stress intensity factors at all studeid crack depths
cycles_int = noCycles; % Crack propagation time to each time instant

% CYCLE-BY-CYCLE EVALUATION
b_evol = [];
b_evol(1) = b_ini;
DK = [];
DK(1) = 0;
i = 1;
b_cur = b_ini;
while b_cur < bc
    DK(i+1) = S((i+1)-(size(S,2).*(floor(i./size(S,2))))).*F(b_cur)...
    *sqrt(pi.*b_cur/b_s); % Stress intensity range (= max)
    dadN(i) = C*DK(i+1)^m; % Crack growth increment [m]
    b_cur = b_cur + dadN(i); % Current crack length [m]
    i = i+1;
    b_evol(i) = b_cur;

    if i>1e6
        disp('non convergent')
        break
    end
end
disp('==== cycle by cycle ====')

disp('Number of cycles to falure according to cycle counting:')
Life_cyc = i-1;
disp(Life_cyc)
disp('Fracture at crack depth [m]:')
disp(b_evol(i-1))
disp('At stress intensity range [MPa(m)^(1/2)]:')
disp(DK(i-1))
disp('----------------------------------------------------')
% Save plot variables
b_cyc = b_evol; % Crack depths in mm
K_cyc = DK; % Stress intensity factors at all studeid crack depths
cycles_cyc = [1:i];

% Plot results

fig=0;

if task == 1
    fig = fig+1;
    figure(fig)
    plot(b_int.*1000,K_int,'r--','LineWidth',2)
    hold on
    plot(b_cyc.*1000,K_cyc,'k.','LineWidth',2)
    hold off
    set(gca,'FontSize',12)
    grid
    title('SIF vs crack depth')
    xlabel('Crack depth [mm]','FontSize',18)
    ylabel('Stress intensity factor [MPa\surdm]','FontSize',18)
legend('Integration','Cycle-by-cycle')
end

fig = fig+1;
figure(fig)
plot(cycles_int,b_int.*1000,'r--','LineWidth',2)
hold on
plot(cycles_cyc,b_cyc.*1000,'k-','LineWidth',2)
hold off
set(gca,'FontSize',12)
grid
title('Crack depth vs cycle')
xlabel('Number of cycles','FontSize',18)
ylabel('Crack depth [mm]','FontSize',18)
legend('Integration','Cycle-by-cycle')
if task == 1
    fig = fig+1;
    figure(fig)
    plot(cycles_int,K_int,'r--','LineWidth',2)
    hold on
    plot(cycles_cyc,K_cyc,'k.','LineWidth',2)
    hold off
    set(gca,'FontSize',12)
    grid
    title('SIF vs cycle')
    xlabel('Number of cycles','FontSize',18)
    ylabel('Stress intensity [MPa\surdm]','FontSize',18)
    legend('Integration','Cycle-by-cycle')
end  

% A.3 TASK 2: DIGITAL TWIN
b_test = [];
N_test = floor(Life_int/6):floor(Life_int/6):Life_int-100; % Inspection cycles
for i = 1:size(N_test,2)
    b_test(i) = interp1(cycles_int,b_int,N_test(i)); % "Ideal" crack depths
end
% Add noise to crack depth data
for i = 1: size(b_test,2)
    pd = makedist('Normal','mu',b_test(i),'sigma',b_test(i)/10);
    b_test(i) = random(pd); % Distorted crack depth data
end
% EVALUATE INSPECTION INTERVALS
for i = 1:size(b_test,2)
    Q = quad(@integrand,b_test(i),bc);
    remainingCycles(i) = Q;
    inspectionCycles(i) = remainingCycles(i)./3;
    totalLife(i) = remainingCycles(i)+N_test(i);
end

% Plot results
fig = fig+1;
figure(fig)
plot(N_test,totalLife,'-','LineWidth',2)
set(gca,'FontSize',12)
hold on
plot(N_test,remainingCycles,'--','LineWidth',2)
plot(N_test,inspectionCycles,':','LineWidth',2)
grid
hold off
legend('Predicted total life', 'Predicted remaining life',...
'Predicted inspection interval')
title('Updated life prediction')
xlabel('Number of cycles')
ylabel('Life estimations')
   
% ESTIMATE MATERIAL PARAMETERS
% Establish m and C by curve fitting
% Evaluate and take logarithm of dadN and DK
for i = 2:size(b_test,2)
    % Average crack growth rate in the inspection interval NOTE dimensions!
    dadN_test(i-1) = (b_test(i)-b_test(i-1))/(N_test(i)-N_test(i-1));
    % Average stress intensity factor range in the inspection interval
    % NOTE dimensions and use that the loading is alternating!
    DK_test(i-1) = (F(b_test(i)).*S.*sqrt(pi.*b_test(i)/b_s)+...
    F(b_test(i-1)).*S.*sqrt(pi.*b_test(i-1)/b_s))/2*1e6;
end
log_dadN_test = log10(dadN_test);
log_DK_test = log10(DK_test);
p = polyfit(log_DK_test, log_dadN_test, 1);
m_fit = p(1);
C_fit = exp(p(2));
disp('m-value predicted by curve fitting (da/dN in mm/cycle and DK in MPa):')
disp(m_fit)
disp('C-value predicted by curve fitting (da/dN in mm/cycle and DK in MPa):')
disp(C_fit)

% Plot results
fig = fig+1;
figure(fig)
plot(log_DK_test,log_dadN_test,'x','LineWidth',2)
% set(gca,'FontSize',18)
hold on
y = polyval(p,log_DK_test);
plot(log_DK_test,y,'-','LineWidth',2)
grid
legend('Measured', 'Curve fit')
title('da/dN vs DK fit')
xlabel('log \Delta K')
ylabel('log da/dN')