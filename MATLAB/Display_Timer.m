% This function displays the elapsed time and estimation of the remaining
% time for the rest of the simulation.
% 
% The function receives the following inputs:
% mode:             A string which specifies the mode of the simulation.
%                   Can be either "Fisher" or "Comparison".
% Elapsed_Time:     The elapsed time so far (in seconds).
% runs:             Number of total simualations runs.
% i:                Index which is associated with the current run.
% X:                This parameter has double meaning. If mode = 
%                   "Comparison", X is a cell array, while each cell is
%                   associated with different run. If mode = "Fisher", X is
%                   a logical variable that indicates whether the variaton
%                   done on the parameter is positive (false) or negative
%                   (true).
% 
function Display_Timer(mode,Elapsed_Time,runs,i,X)
% Calculate the percentage of calculations performed so far
if strcmp(mode,'Fisher')
    ind = 2*i+1*X;
elseif strcmp(mode,'Comparison')
    I = zeros(1,runs);
    for j=1:runs
        I(j) = X{j}('I');
    end
    ind = sum(I(1:i))/mean(I);
end
percents = 100*ind/runs;
% Calculate the remaining time (in seconds)
Remaining_Time = Elapsed_Time*(100/percents-1);
% Make the string for the elapsed time
Elapsed_Hours = floor(Elapsed_Time/3600);
Elapsed_Minutes = floor((Elapsed_Time-3600*Elapsed_Hours)/60);
Elapsed_Seconds = floor(Elapsed_Time-3600*Elapsed_Hours-60*Elapsed_Minutes);
if Elapsed_Hours == 0
    Elapsed_Hours = '';
else
    Elapsed_Hours = [num2str(Elapsed_Hours),' hours, '];
end
if Elapsed_Minutes == 0
    Elapsed_Minutes = '';
else
    Elapsed_Minutes = [num2str(Elapsed_Minutes),' minutes, '];
end
if Elapsed_Seconds == 0
    Elapsed_Seconds = '';
else
    Elapsed_Seconds = [num2str(Elapsed_Seconds),' seconds'];
end
Elapsed_Time_str = [Elapsed_Hours,Elapsed_Minutes,Elapsed_Seconds];
% Make the string for the remaining time
Remaining_Hours = floor(Remaining_Time/3600);
Remaining_Minutes = floor((Remaining_Time-3600*Remaining_Hours)/60);
Remaining_Seconds = floor(Remaining_Time-3600*Remaining_Hours-60*Remaining_Minutes);
if Remaining_Hours == 0
    Remaining_Hours = '';
else
    Remaining_Hours = [num2str(Remaining_Hours),' hours, '];
end
if Remaining_Minutes == 0
    Remaining_Minutes = '';
else
    Remaining_Minutes = [num2str(Remaining_Minutes),' minutes, '];
end
if Remaining_Seconds == 0
    Remaining_Seconds = '';
else
    Remaining_Seconds = [num2str(Remaining_Seconds),' seconds'];
end
Remaining_Time_str = [Remaining_Hours,Remaining_Minutes,Remaining_Seconds];
% Print the elapsed and remaining time
fprintf('\n');
disp([num2str(round(percents)),'% of calculations have been completed.']);
disp(['Elapsed time: ',Elapsed_Time_str,'.']);
if percents ~=100
    disp(['Estimation of remaining time: ',Remaining_Time_str,'.']);
end