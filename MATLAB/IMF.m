% This function generates the Initial Mass Function (IMF).
% 
% The function receives the following inputs:
% M:        a vector containing linearly spaced mass values (in units of 
%           solar masses).
% params:   a dictionary containing numerical parameters. See documantation
%           at Main.m for more details.
% str_params: a dictionary containing string parameters. See documantation
%             at Main.m for more details.
% 
% The function returns the following outputs:
% f_0:    a vector containing samples of the IMF at specific mass values 
%         which are given by the vector M.
% mean:   the IMF's mean (in units of M0).
% 
function [f_0,mean] = IMF(M,params,str_params)
% Extract parameters from the dictionaries
i_max = params('i_max');
IMF_TYPE = str_params('IMF_TYPE');
alpha = params('alpha');
Mmin = params('Mmin');
Mmax = params('Mmax');
mu_pp = params('mu_pp');
sigma_pp = params('sigma_pp');
lambda = params('lambda');
mu_NS = params('mu_NS');
sigma_NS = params('sigma_NS');
N_SBH = params('N_SBH');
N_NS = params('N_NS');
t_NS = params('t_NS');
% Generate the IMF based on the value of IMF_TYPE
switch IMF_TYPE
    case 'Seed'
        f_0 = (M==M(1));
    case 'Model'
        f_0 = power(M,-alpha).*(M>=Mmin & M<=Mmax);
end
% Normalize the IMF
f_0 = f_0/sum(f_0);
% Add the PPSNe contribution
if strcmp(IMF_TYPE,'Model')
    f0_PPSN = exp(-(M-mu_pp).^2/(2*sigma_pp^2)).*(M>=Mmin);
    f0_PPSN = f0_PPSN/sum(f0_PPSN);
    f_0 = (1-lambda)*f_0 + lambda*f0_PPSN;
    % Add the NSs contribution and normalize the IMF to N_SBH+N_NS
    if t_NS == 0
        f0_NS = exp(-(M-mu_NS).^2/(2*sigma_NS^2)).*(M>=0);
        f0_NS = f0_NS/sum(f0_NS);
        f_0 = (N_SBH*f_0 + N_NS*f0_NS);
    else
        f_0 = N_SBH*f_0;
    end
elseif strcmp(IMF_TYPE,'Seed')
    f_0 = N_SBH*f_0;
end
% Calculate the IMF's mean
mean = sum(f_0.*(1:i_max));