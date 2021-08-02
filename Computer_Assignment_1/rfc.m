 function [mid,range]=rfc(stress_history)

% Rainflow count by the ASTM method
% presumes a time series with ONLY extreme values
% [1 0 1] works
% [1 0 0 1] does not work
% [1 0 -1 1] does not work
% Should be fixed
%
% history - vector with time history in the form of extreme values--indata
% mid     - vector with mid values
% range   - vector with range

% Check if start and end value of indata are the same (true cycle)

true_cycle=history(size(history,2))==history(1); % 1 if true
if true_cycle % If first value equals last value
    [sig_max,sig_max_pos] = max(history); % Find position of largest element
    % re-organize series starting with largest max
    if sig_max_pos ~= 1 % If largest value is the first, rearrangement is not needed
        history(end) = []; % Remove last value to later avoid dubplication
        his_start = history(1:sig_max_pos-1); % Extract all elements infront of largest history(1:sig_max_pos-1) = []; % Remove all elements infront of largest element
        history = [history his_start sig_max]; % Add extracted values at the end
    end
else
    disp('This is not a true cycle.')
end


% Start identification of cycles
k=1; % counter for identified stress cycles
i=2; % counter to use in the cycle counting
safety = 0; % Avoid infinite loops due to malformatted input values
max_safety = size(history,2)+1e3;
while size(history,2)>3 % Repeat thorugh history until three items remain
    safety = safety + 1;
    if abs(history(i-1)-history(i)) <= abs(history(i)-history(i+1)) % Cycle
        mid(k) = 1/2*(history(i-1)+history(i)); % Mid value of cycle
        range(k) = abs((history(i-1)-history(i))); % Range of cycle
        history(i) = []; % Remove affected cycle values
        history(i-1) = []; % Remove affected cycle values
        k = k+1;
        i = 1; % Start from beginning
    end
    i = i+1;
    if safety >= max_safety % You have looped through more than 1000 cycles
        disp ('Long loop - check your indata')
        break
    end
end
% Final cycle (only three values remain)
mid(k)= (S(2)+S(1))/2; 
range(k)= abs(S(2)-S(1)); 