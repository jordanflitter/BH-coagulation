% This function calculates the Fisher matrix and plots the confidence
% bounds ellipses of the free parameters.

% The function receives the following inputs:
% f_det_M1M2z_0: a 3D array of the events density 3D function with respect 
%                to M1, M2 and z (which corresponds to the fiducial
%                parameters). The pages correspond to the vector z_vec
%                while the rows and the columns of this matrix correspond
%                to the vector M_vec.
% M_vec:         a vector of mass values (in units of solar masses).
% z_vec:         a vector of redshift values.
% free_params:   an array of strings which specify the free parameters.
% params:        a dictionary containing numerical parameters. See 
%                documantation at Main.m for more details.
% str_params:    a dictionary containing string parameters. See 
%                documantation at Main.m for more details.
% Elapsed_Time:  The program's elapsed time so far (in seconds).
% 
function Fisher_Analysis(f_det_M1M2z_0,M_vec,z_vec,free_params,params,str_params,Elapsed_Time)
% Extract various quantities from the inputs
M0 = params('M0');
delta = params('delta');
Npar = length(free_params);
Mmin = params('Mmin');
Nbins = params('Nbins');
LIGO_RUN = str_params('LIGO_RUN');
save_workspace = params('save_workspace');
save_figures = params('save_figures');
% Make sure to not override the already saved figures and workspace (which
% correspond to the fiducial parameters)
params('plot_coagulation') = 0;
params('plot_M1_dist') = 0;
params('plot_M2_dist') = 0;
params('plot_q_dist') = 0;
params('plot_z_dist') = 0;
params('save_workspace') = 0;
% Display the elapsed and remaining time
Display_Timer('Fisher',Elapsed_Time,2*Npar+1,0,true);
% Set experiment
if strcmp(LIGO_RUN,'O3a')
    T = 177.3/365; %year
elseif strcmp(LIGO_RUN,'O5')
    T = 0.8*6; %years
end
% Calculate f_det(M1,M2) for the fiducial values
f_det_M1M2_0 = trapz(z_vec,f_det_M1M2z_0,3);
% Bin the 2D distribution
i_gap = find(sum(f_det_M1M2_0)>0,1,'First');
ibins_edges = round(logspace(log10(M_vec(i_gap)),log10(M_vec(end)),Nbins+1)/M0);
[N_obs_0,N_ind] = Binning_2D(f_det_M1M2_0,ibins_edges);
% Allocate the data matrix and insert the results for the fiducial values 
% in the appropriate row
N_obs_matrix = zeros(2*Npar+1,length(N_obs_0));
N_obs_matrix(Npar+1,:) = N_obs_0;
% At each iteration, only one of the parameters is multiplied by a factor
% of 1+-delta
for i = 1:Npar
%   Calculate the events density function for the positively modified
%   parameter
    disp(['Now simulating for (1+delta)',char(free_params(i)),'...']);
    params_p = containers.Map(params.keys,params.values);
    if params(char(free_params(i))) ~= 0
        params_p(char(free_params(i))) = params(char(free_params(i)))*(1+delta);
    else
        params_p(char(free_params(i))) = delta;
    end
    tic;
    f_det_M1M2z_p = BHMF(params_p,str_params,1);
    f_det_M1M2_p = trapz(z_vec,f_det_M1M2z_p,3);
    N_obs_p = Binning_2D(f_det_M1M2_p,ibins_edges,N_ind);
    Elapsed_Time = Elapsed_Time + toc;
    Display_Timer('Fisher',Elapsed_Time,2*Npar+1,i,false);
%   Calculate the events density function for the negatively modified
%   parameter
    disp(['Now simulating for (1-delta)',char(free_params(i)),'...']);
    params_m = containers.Map(params.keys,params.values);
    if params(char(free_params(i))) ~= 0
        params_m(char(free_params(i))) = params(char(free_params(i)))*(1-delta);
    else
        params_m(char(free_params(i))) = -delta;
    end
    tic;
    f_det_M1M2z_m = BHMF(params_m,str_params,1);
    f_det_M1M2_m = trapz(z_vec,f_det_M1M2z_m,3);
    N_obs_m = Binning_2D(f_det_M1M2_m,ibins_edges,N_ind);
    Elapsed_Time = Elapsed_Time + toc;
    Display_Timer('Fisher',Elapsed_Time,2*Npar+1,i,true);
%   Insert the results in the appropriate rows of the data matrix
    N_obs_matrix(i,:) = N_obs_p; N_obs_matrix(Npar+1+i,:) = N_obs_m;
end
% Calculate the Fisher Matrix and the covariance matrix
disp('Now calculating the Fisher Matrix and making the stair plot...');
fiducial_params = zeros(1,Npar);
for i = 1:Npar
    fiducial_params(i) = params(char(free_params(i)));
end
[F,scales] = Fisher_Matrix(fiducial_params,N_obs_matrix,delta);
cov_matrix = inv(F);
% Rename the free parameters
names = NameOrganizer(free_params);
% Save the results
if save_workspace
    load([pwd,'/Run 1/Run 1 - Workspace']);
    save([pwd,'/Run 1/Run 1 - Workspace']);
end
% Make the stair plot
out = MakeStairPlot(cov_matrix,fiducial_params./scales,names,['aLIGO: ',num2str(T),'yrs'],scales);
annotation('textbox',[0.4,0.5,0.5,0.5],'String',['Number of total events: ',num2str(floor(sum(N_obs_0)))],...
    'FitBoxToText','on','EdgeColor','none','FontName','Times New Roman','FontSize',18);
fprintf('\n');
% Save the stair plot
if save_figures
    savefig(out(1),[pwd,'/Run 1/Run 1 - Stair Plot']);
end