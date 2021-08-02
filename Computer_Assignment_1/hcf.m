function [Nf]=hcf(mid, range, criterion, fatlim)

% Evaluate HCF life
% mid        - stress mid values 
% range      - stress ranges
% criterion  - HCF-criterion
% Nf         - Fatigue life--output data

% S_ult=800; %ultimate tensile strength
% WOHLER CURVE
% Give two combinations of stress magnitudes and corresponding fatigue life
% to derive a Wohler curve
Nf_mtrl = [10^3 10^6];
Sa_mtrl = [453.7719 167.8956];
% Take the logarithm of (at least) the fatigue life
logNf_mtrl = log10(Nf_mtrl);

% Curve fit Wohler curve or alternatively derive Basquin coefficeints
p = polyfit(logNf_mtrl,Sa_mtrl,1);

% STRESS MEASURES
Sa = range/2;
Smax = mid+Sa;

if criterion == 'SWT' % Smith-Watson-Topper
    if Smax > 0
        Saeq = sqrt(Smax.*Sa)
    else
        Saeq = 0
    end
    D = 0;
    logNf_vec=[]; % vector with fatigue lives for plotting
    if fatlim == 'y'
        % The damage evaluation accounts for a fatigue limit!
        for i=1:size(Saeq,2)  % Loop through the stress cycles
            if Saeq(i)> Sa_mtrl(2)
                logNf = (Saeq(i)-p(2))/p(1);
                D =D+1/(10^logNf);
                logNf_vec = [logNf_vec,logNf];
            else
                logNf =Inf;
                D =D+1/(10^logNf);
                logNf_vec = [logNf_vec,logNf];
            end
        end
    else
        % The damage evaluation does not account for a fatigue limit!
        for i=1:size(Saeq,2)  % Loop through the stress cycles
            logNf = (Saeq(i)-p(2))/p(1);
            D = D+1/(10^logNf);
            logNf_vec = [logNf_vec,logNf];
        end
    end
else
    disp('Curently, the criterion')
    disp(criterion)
    disp('is not implemented')
end
Nf = 1/D;

% Plot
figure(3)
plot(logNf_mtrl,Sa_mtrl,'-k');
hold on
for i=1:size(Saeq,2)
    Nf_plot_data = [0,logNf_vec(i),logNf_vec(i)];
    Saeq_plot_data = [Saeq(i), Saeq(i), 0];
    plot(Nf_plot_data,Saeq_plot_data,'-r')
end
hold off
