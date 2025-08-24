% This function calculates the events density 2D function with respect to
% M1 and M2. Figures and workspace are saved according to the user's desire
% - See documantation at Main.m for more details.
% 
% The function receives the following inputs:
% params:   a dictionary containing numerical parameters. See documantation
%           at Main.m for more details.
% str_params: a dictionary containing string parameters. See documantation
%             at Main.m for more details.
% run:      the index of the current scenario that is being simulated.
% 
% The function returns the following outputs:
% f_det_M1M2z:  a 3D array of the events density 3D function with respect 
%               to M1, M2 and z. The pages correspond to the vector z_vec 
%               while the rows and the columns of each page correspond 
%               to the vector M_vec.
% M_vec:        a vector of mass values (in units of solar masses).
% z_vec:        a vector of redshift values.
% 
function [f_det_M1M2z,M_vec,z_vec] = BHMF(params,str_params,run)
% Extract the parameters from the dictionaries
M0 = params('M0');
i_max = params('i_max');
I = params('I');
Nt = params('Nt');
Nbins = params('Nbins');
observation_flag = params('observation_flag');
plot_coagulation = params('plot_coagulation');
plot_M1_dist = params('plot_M1_dist');
plot_M2_dist = params('plot_M2_dist');
plot_q_dist = params('plot_q_dist');
plot_z_dist = params('plot_z_dist');
save_workspace = params('save_workspace');
save_figures = params('save_figures');
SYM_TYPE = str_params('SYM_TYPE');
IMF_TYPE = str_params('IMF_TYPE');
LIGO_RUN = str_params('LIGO_RUN');
alpha = params('alpha');
Mmin = params('Mmin');
Mmax = params('Mmax');
mu_pp = params('mu_pp');
sigma_pp = params('sigma_pp');
lambda = params('lambda');
mu_NS = params('mu_NS');
sigma_NS = params('sigma_NS');
X = params('X');
R11 = params('R11');
R_0 = params('R_0');
beta = params('beta');
beta2 = params('beta2');
gamma = params('gamma');
t_c = params('t_c');
f_loss = params('f_loss');
M_esc = params('M_esc');
v_esc = params('v_esc');
N_SBH = params('N_SBH');
N_NS = params('N_NS');
t_NS = params('t_NS');
% Generate the IMF
M = M0*(1:i_max);
if ~strcmp(SYM_TYPE,'Test')
    [~,mean] = IMF(M,params,str_params);
    i_max = round(mean);
    params('i_max') = i_max;
    M = M0*(1:i_max);
end
f_0 = IMF(M,params,str_params);
% Set the size of the rate matrices.
if strcmp(LIGO_RUN,'O3a')
    load('P_values_aLIGO_O3a');
elseif strcmp(LIGO_RUN,'O5')
    load('P_values_aLIGO_O5');
end
ind_max = min(round(M_vec(end)/M0),i_max);
% Calculate the static channel density rate matrix R_stat(M1,M2)
[i_mat,j_mat] = meshgrid(1:ind_max);
R2 = power(min(i_mat./j_mat,j_mat./i_mat),beta2);
if t_NS~=0
    f0_NS = exp(-(M-mu_NS).^2/(2*sigma_NS^2)).*(M>=0);
    f0_NS = N_NS*f0_NS/sum(f0_NS);
    f_stat = f_0 + f0_NS;
else
    f_stat = f_0;
end
f_stat = f_stat/sum(f_stat);
R_stat = R_0*(f_stat(1:ind_max)'*f_stat(1:ind_max)).*R2; 
% Evolve the cluster's BHMF by solving the coagulation equation. 
% Calculate Gamma_i,j(t)
if X~=0
    [f_src,Gamma,Gamma_tot,moments,P_500,F1,F2,F3,F4] = Coagulation(M,f_0,params,str_params);
    Gamma = Gamma(1:ind_max,1:ind_max,:);
else
    Gamma = zeros(ind_max,ind_max,Nt);
end
% Calculate the detected merger events distributions
if observation_flag
%   Calculate the events denstity function with respect to M1, M2 and z.
%   f_det(M1,M2,z) is the combined contributions of both the static and
%   dynamic channels to the density function.
%   f_det_stat(M1,M2,z) is the contribution from only the static channel
%   to the density function.
    [f_det_M1M2z,f_det_M1M2z_stat,M_vec,z_vec] = Bias(M(1:ind_max),Gamma,R_stat,params,str_params);
%   Bin the distributions and plot them
    if plot_M1_dist
        [N_obs_M1,Mbins] = Binning_1D(f_det_M1M2z,M_vec,z_vec,'M1',Nbins,LIGO_RUN);
        N_obs_M1_stat = Binning_1D(f_det_M1M2z_stat,M_vec,z_vec,'M1',Nbins,LIGO_RUN,Mbins);
        F5 = Nobs_plotter(Mbins,N_obs_M1_stat,N_obs_M1,'M1',LIGO_RUN);
    end
    if plot_M2_dist
        [N_obs_M2,Mbins] = Binning_1D(f_det_M1M2z,M_vec,z_vec,'M2',Nbins,LIGO_RUN);
        N_obs_M2_stat = Binning_1D(f_det_M1M2z_stat,M_vec,z_vec,'M2',Nbins,LIGO_RUN,Mbins);
        F6 = Nobs_plotter(Mbins,N_obs_M2_stat,N_obs_M2,'M2',LIGO_RUN);
    end
    if plot_q_dist
        [N_obs_q,qbins] = Binning_1D(f_det_M1M2z,M_vec,z_vec,'q',Nbins,LIGO_RUN);
        N_obs_q_stat = Binning_1D(f_det_M1M2z_stat,M_vec,z_vec,'q',Nbins,LIGO_RUN);
        F7 = Nobs_plotter(qbins,N_obs_q_stat,N_obs_q,'q',LIGO_RUN);
    end
    if plot_z_dist
        [N_obs_z,zbins] = Binning_1D(f_det_M1M2z,M_vec,z_vec,'z',Nbins,LIGO_RUN);
        N_obs_z_stat = Binning_1D(f_det_M1M2z_stat,M_vec,z_vec,'z',Nbins,LIGO_RUN);
        F8 = Nobs_plotter(zbins,N_obs_z_stat,N_obs_z,'z',LIGO_RUN);
    end
else
    f_det_M1M2z = 0;
end
% Save figures
if save_figures
    if ~exist(['Run ',num2str(run)], 'dir')
        mkdir(['Run ',num2str(run)]);
    end
    if plot_coagulation && X ~= 0
        savefig(F1,[pwd,'/Run ',num2str(run),'/Run ',num2str(run),' - Coagulation']);
        savefig(F2,[pwd,'/Run ',num2str(run),'/Run ',num2str(run),' - N_tot']);
        savefig(F3,[pwd,'/Run ',num2str(run),'/Run ',num2str(run),' - M_tot']);
        savefig(F4,[pwd,'/Run ',num2str(run),'/Run ',num2str(run),' - Total Merger Rate']);
    end
    if observation_flag && plot_M1_dist
        savefig(F5,[pwd,'/Run ',num2str(run),'/Run ',num2str(run),' - M1 Distribution']);
    end
    if observation_flag && plot_M2_dist
        savefig(F6,[pwd,'/Run ',num2str(run),'/Run ',num2str(run),' - M2 Distribution']);
    end
    if observation_flag && plot_q_dist
        savefig(F7,[pwd,'/Run ',num2str(run),'/Run ',num2str(run),' - q Distribution']);
    end
    if observation_flag && plot_z_dist 
        savefig(F8,[pwd,'/Run ',num2str(run),'/Run ',num2str(run),' - z Distribution']);
    end
end
% Save workspace
if save_workspace
    clear F1 F2 F3 F4 F5 F6 F7 F8;
    if ~exist(['Run ',num2str(run)], 'dir')
        mkdir(['Run ',num2str(run)]);
    end
    save([pwd,'/Run ',num2str(run),'/Run ',num2str(run),' - Workspace']);
end