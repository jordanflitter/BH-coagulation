% This function defines the appropriate names (or labels) of the parameters
% to be presented at the stair plot.
% 
% The function receives the following inputs:
% free_params: an array of strings which specify the free parameters.
% 
% The function returns the following output:
% names: an array of strings of the appropriate parameters names to be
%        presented at the stair plot.
% 
function names = NameOrganizer(free_params)
Npar = length(free_params);
names = cell(1,Npar);
for i = 1:Npar
    if (strcmp(free_params(i),'alpha') || strcmp(free_params(i),'beta') || strcmp(free_params(i),'lambda') || strcmp(free_params(i),'gamma'))
        names(i) = {['\',char(free_params(i))]};
    elseif strcmp(free_params(i),'Mmin')
        names(i) = {'M_{min}'};
    elseif strcmp(free_params(i),'Mmax')
        names(i) = {'M_{max}'};
    elseif strcmp(free_params(i),'mu_pp')
        names(i) = {'\mu_{pp}'};
    elseif strcmp(free_params(i),'sigma_pp')
        names(i) = {'\sigma_{pp}'};
    elseif strcmp(free_params(i),'mu_NS')
        names(i) = {'\mu_{NS}'};
    elseif strcmp(free_params(i),'sigma_NS')
        names(i) = {'\sigma_{NS}'};
    elseif strcmp(free_params(i),'R11')
        names(i) = {'R_{11}'};
    elseif strcmp(free_params(i),'beta2')
        names(i) = {'\beta_2'};
    elseif strcmp(free_params(i),'f_loss')
        names(i) = {'f_{loss}'};
    elseif strcmp(free_params(i),'M_esc')
        names(i) = {'M_{esc}'};
    elseif strcmp(free_params(i),'v_esc')
        names(i) = {'v_{esc}'};
    elseif strcmp(free_params(i),'N_SBH')
        names(i) = {'N_{SBH}'};
    elseif strcmp(free_params(i),'N_NS')
        names(i) = {'N_{NS}'};
    elseif strcmp(free_params(i),'t_NS')
        names(i) = {'t_{NS}'};
    else
        names(i) = free_params(i);
    end
end