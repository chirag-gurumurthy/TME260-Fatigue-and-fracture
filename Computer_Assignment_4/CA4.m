% A.1 StressMaker.m

graphicx color
% Creates a stress matrix and calls for Dang Van evaluation
clear all
close all
clc
% stress - [s11 s22 s33 s12 s13 s23; time 1
% s11 s22 s33 s12 s13 s23; time 2
% etc ]

% Crossland fatigue limit
SeC=318*1e6;

i_vector=[0:0.01:4.*pi];
stress=zeros(size(i_vector,2),6);
% Insert stress components below
Ktx = 2;
Ktxy = 1.5;
Sx = 97.8*1e6*Ktx;
Sy = 0;
Txy = 48.9*1e6*Ktxy;
% Insert residual stress
Sx_res = 100*1e6;
% Insert phase angle
n=(0:7);
phi = n.*pi./8;
for ii=1:size(n,2)
    j=0;
    for i = i_vector
        j=j+1;
        stress(j,:)=[Sx.*sin(i)+Sx_res Sy.*sin(i) 0 Txy.*sin(i+phi(ii)) 0 0];
    end
% Insert material parameter
    cC=2/3;
    SeqC(ii)=Crossland(stress,cC);
    
end
SF=SeC/max(SeqC);

% A.2 Crossland.m

function MaxC = Crossland(stress, cC)

% Evaluate Crossland stress
% stress - [s11 s22 s33 s12 s13 s23; time 1
%           s11 s22 s33 s12 s13 s23; time 2
%           etc ]
% cC - material parameter in the crossland criterion
% Evaluate hydrostatic stress

sh = (stress(:,1)+stress(:,2)+stress(:,3))./3;
sh_max = max(sh);

% Deviatoric stress tensor
sdev = stress;
for i = 1:3
    sdev(:,i) = sdev(:,i)-sh;
end

% Find
% min { max(t) [J2(sdev-smid)]}
% Find midvalue without an extra condition that smid is deviatoric
% smid = fminsearch(@(smid) J2(smid,sdev),(max(sdev)+min(sdev))/2);
% Find midvalue while assuring that smid is deviatoric
% Use s5tosdev to translate from a vector with 5 values to deviatoric
% tensor with 6 values
% Use sdevtos5 to go from a deviatoric stress tensor with 6 values to a 5
% value vector
smid5s = fminsearch(@(s5) J2(s5tosdev(s5),sdev),sdevtos5((max(sdev)+min(sdev))/2));
smid=s5tosdev(smid5s);

% Evaluate Sa = Sdev-Smid
for i=1:size(sdev,1)
    Sa(i,:) = sdev(i,:)-smid;
end

% Make stress tensor (3x3) from the state of str5ess at the current instant
% in time
for i=1:size(Sa,1)
    s = [Sa(i,1) Sa(i,4) Sa(i,5); Sa(i,4) Sa(i,2) Sa(i,6); Sa(i,5) Sa(i,6) Sa(i,3)];
    % Evaluate SUM(SijSji)
    t=0;
    for k=1:3
        for l=1:3
            t=t+(s(k,l)*s(k,l));
        end
    end
    svM(i) = sqrt(3/2*t);
end

% Evaluate the Crossland stress
seqC = svM+cC.*sh_max;

% Create a time coordinate to plot the Crossland stress against
x=1:size(seqC,2);

% Plot stress components
figure()
cax = newplot();
set(cax,'FontName','Times','FontSize',14);
subplot(2,1,1)
plot(x,stress(:,1)','-',x,stress(:,4)','-.',x,stress(:,2)','--','LineWidth',2)
grid
legend('S1','S4','S2')
title('Stress components')
ylabel('Stress')
xlabel('time increment')

% Plot Crossland stress
% figure()
% cax = newplot()
% set(cax,'FontName','Times','FontSize',14)
subplot(2,1,2)
plot(x,seqC,'-',x,svM,'-.',x,sh','--','LineWidth',2)
grid
legend('Crossland equivalent stress','von Mises amplitude','Hydrostatic stress')
title('Crossland equivalent stress')
ylabel('Stress')
xlabel('time increment')

% Print Maximum Crossland equivalent stress
disp('========= run =========')
MaxC = max(seqC)
MaxvonMisesAmp = max(svM)
disp('===== end of run ======')


% A.3 sdevtos5.m

function s5=sdevtos5(s)
% Maps time history of a 6D stress tensor to a 5D stress tensor
% by removing Sz

s5=[s(:,1) s(:,2) s(:,4:6)];
    

% A.4 s5tosdev.m

function s=s5tosdev(s5)
% Maps time history of a 5D stress tensor to a DEVIATORIC 6D stress tensor
% This is assured by setting Sz = Sy since both are zero before Sh is
% subtracted

s=[s5(:,1) s5(:,2) -s5(:,1)-s5(:,2) s5(:,3:5)];


% A.5 J2.m

function f=J2(smid,sdev)
% max J2 (or von Mises norm of sa)
% smid is the evaluated mid value of the deviatoric stress tensor
% (dimension: 1 x 6)
% sdev contains the deviatoric stress tensor for all time increments
% (dimension: time increments x 6)

% Create empty matrix with deviatoric stress "amplitudes"
sa = zeros(size(sdev));
% Create vector sa with deviatoric stress "amplitude"
for i=1:size(sdev,1)
    sa(i,:) = sdev(i,:)-smid;
end

% make stress tensor
% Make stress tensor (3x3) from the state of str5ess at the current instant
% in time
for i=1:size(sa,1)
    s = [sa(i,1) sa(i,4) sa(i,5); sa(i,4) sa(i,2) sa(i,6); sa(i,5) sa(i,6) sa(i,3)];
    % Evaluate SUM(SijSji)
    t=0;
    for k=1:3
        for l=1:3
            t=t+(s(k,l)*s(k,l));
        end
    end
    % Normalise
    sJ2(i) = sqrt(3/2*t);
end
% Evaluate the maximum J2 magnitude
f=max(sJ2);

    
    
    
    
    