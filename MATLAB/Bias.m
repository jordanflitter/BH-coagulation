% This function calculates the merger events density function
% f_det(M1,M2,z).
% 
% The function receives the following inputs:
% M:          a vector containing linearly spaced mass values
%            (in units of solar masses).
% Gamma:      a 3D array of the merger rates (in units of year^-1).
%             In each page, the rows and columns corrspond to different 
%             mass values (which are given by the vector M), while 
%             different pages correspond to different time samples (at 
%             each t_c/Nt).
% R_stat:     a matrix of the merger rates in the static channel (in 
%             units of Gpc^-3*year^-1). The rows and columns corrspond to 
%             different mass values (which are given by the vector M).
% params:   a dictionary containing numerical parameters. See documantation
%           at Main.m for more details.
% str_params: a dictionary containing string parameters. See documantation
%             at Main.m for more details.
% 
% The function returns the following outputs:
% f_det_M1M2z:  a 3D array of the events density 3D function with respect 
%               to M1, M2 and z which include contributions from both 
%               static and dynamic channels. The pages correspond to the 
%               vector z_vec while the rows and the columns of this matrix 
%               correspond to the vector M_vec.
% f_det_M1M2z_stat:  a 3D array of the events density 3D function with 
%               respect to M1, M2 and z which include the contribution from 
%               only the static channels. The pages correspond to the 
%               vector z_vec while the rows and the columns of this matrix 
%               correspond to the vector M_vec.
% M_vec:        a vector of mass values (in units of solar masses).
% z_vec:        a vector of redshift values.
% 
function [f_det_M1M2z,f_det_M1M2z_stat,M_vec,z_vec] = Bias(M,Gamma,R_stat,params,str_params)
% Extract parameters from the dictionaries
M0 = params('M0');
Mmin = params('Mmin');
X = params('X');
LIGO_RUN = str_params('LIGO_RUN');
% Set experiment
if strcmp(LIGO_RUN,'O3a')
    load('P_values_aLIGO_O3a');
    T = 177.3/365; %year
elseif strcmp(LIGO_RUN,'O5')
    load('P_values_aLIGO_O5');
    T = 0.8*6; %years
end
% Make a grid for P(M1,M2,z)
[M1_mat,M2_mat,z_mat] = ndgrid(M_vec,M_vec,z_vec);
P_grid = griddedInterpolant(M1_mat,M2_mat,z_mat,P_values,'linear','none');
% Redefine M_vec (in order to sample better below the lower mass cutoff)
i_min = find(abs(M_vec-Mmin)==min(abs(M_vec-Mmin)),1,'last');
M_vec = [M_vec(1):M0:M_vec(i_min),M_vec(i_min+1):1:M_vec(end)];
% Interpolate P at M_vec and z_vec values
[M1_mat,M2_mat,z_mat] = ndgrid(M_vec,M_vec,z_vec);
P = P_grid(M1_mat,M2_mat,z_mat); P(isnan(P))=0;
% Calculate A_dyn(z)=4*pi*c*chi(z)^2*nGC(z)/((1+z)*H(z)), 
% A_stat(z)=4*pi*c*chi(z)^2*(1+z)^3/((1+z)*H(z)) and z(t) where t is
% accosicated with the time samples of Gamma_i,j(t)
[A_dyn,A_stat,z_t] = Cosmo(M_vec,z_vec,params);
% Interpolate Gamma at M_vec and z_vec values
[M1_mat,M2_mat,z_mat] = ndgrid(M,M,z_t(end:-1:1));
Gamma_grid = griddedInterpolant(M1_mat,M2_mat,z_mat,Gamma(:,:,end:-1:1)/(M0^2));
[M1_mat,M2_mat,z_mat] = ndgrid(M_vec,M_vec,z_vec);
Gamma_z = Gamma_grid(M1_mat,M2_mat,z_mat);
% Interpolate R_stat at M_vec and z_vec values
[M1_mat,M2_mat,z_mat] = ndgrid(M,M,z_vec);
R_stat_grid = griddedInterpolant(M1_mat,M2_mat,z_mat,repmat(R_stat,1,1,length(z_vec))/(M0^2));
[M1_mat,M2_mat,z_mat] = ndgrid(M_vec,M_vec,z_vec);
R_stat_z = R_stat_grid(M1_mat,M2_mat,z_mat);
% Calculate f_det(M1,M2,z)
f_det_M1M2z_dyn = T*A_dyn.*Gamma_z.*P;
f_det_M1M2z_stat = T*A_stat.*R_stat_z.*P;
f_det_M1M2z = X*f_det_M1M2z_dyn + (1-X)*f_det_M1M2z_stat;