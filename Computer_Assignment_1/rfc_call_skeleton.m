% RFC-call
% Read input data in the form of max-min-stresses
clear all
close all
clf

% RAINFLOW COUNT
inputdata_86
[mid,range]=rfc(stress_history); % Call the function that does the RFC 

% Check that the largest range is captured
max_range = range(end);
if max_range ~= max(range)
    disp('Missed largest range ');
end

amplitude =range/2;
figure(1)
errorbar(mid,amplitude);
figure(2)
h = normplot(amplitude);


% FATIGUE LIFE EVALUATION
fatlim = 'y' % Account for the fatigue limit? (y / n)
criterion = 'SWT' % Mid stress criterion
% Parameters: mid stress, stress range, criterion, account for fatigue limit
[Nf]=hcf(mid, range, criterion, fatlim);
disp('=====')
disp('Fatigue life evaluation using and SN-curve and the SWT mid stress criterion.')
if fatlim == 'y'
    disp('The fatigue limit is accounted for.')
else
    disp('The fatigue limit is not accounted for.')
end
disp('The component can sustain ')
disp(num2str(Nf))
disp('repetitions of the current load sequence, which contains')
disp(num2str(length(amplitude)))
disp('stress cycles. Consequently, the component can in total sustain')
disp(num2str(length(amplitude).*Nf))
disp('stress cycles.')
disp('')
disp('The equivalent damage of the load sequence is ')
D =1/Nf;
disp(num2str(D))
disp('[-]')