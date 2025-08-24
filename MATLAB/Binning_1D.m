% This function returns the detected merger events 1D distribution with 
% respect to a chosen variable (M1, M2, q or z).
% 
% The function receives the following inputs:
% f_det_M1M2z:  a 3D array of the events density 3D function with respect 
%               to M1, M2 and z. The pages correspond to the vector z_vec 
%               while the rows and the columns of this matrix correspond 
%               to the vector M_vec.
% M_vec:        a vector of mass values (in units of solar masses).
% z_vec:        a vector of redshift values.
% x_str:        a string which specifies the corresponding distribution. 
%               The options are 'M1', 'M2', 'q' or 'z'.
% Nbins:        Number of bins to be shown in the events distributions.
% LIGO_RUN:     a string which specifies the corresponding experiment.
% Mbins:        an optional input. It's a vector of the mass bins' centers 
%               (in units of solar masses).
% 
% The function returns the following outputs:
% N_obs_x:    a vector containing number of detected events at each bin.
% xbins:      a vector of the bins' centers (in units of solar masses in 
%             case x_str is either 'M1' or 'M2').
% 
function [N_obs_x,xbins] = Binning_1D(f_det_M1M2z,M_vec,z_vec,x_str,Nbins,LIGO_RUN,Mbins)
% Calculate relative errors.
if strcmp(LIGO_RUN,'O3a')
    load('O3a_data');
    sigma_M1 = median(0.5*(M1_data(:,2)+M1_data(:,3))./M1_data(:,1));
    sigma_M2 = median(0.5*(M2_data(:,2)+M2_data(:,3))./M2_data(:,1));
    sigma_q = sqrt((sigma_M1)^2+(sigma_M2)^2);
    sigma_z = median(0.5*(z_data(:,2)+z_data(:,3))./z_data(:,1));
elseif strcmp(LIGO_RUN,'O5')
    sigma_M1 = 0.05; sigma_M2 = 0.05; sigma_z = 0.05;
    sigma_q = sigma_M1*sqrt(2);
end
% Calculate f_det(M1,M2)
f_det_M1M2 = trapz(z_vec,f_det_M1M2z,3);
i_gap = find(sum(f_det_M1M2)>0,1,'First');
%% Calculate either N_obs(M1) or N_obs(M2)
if strcmp(x_str,'M1') || strcmp(x_str,'M2')
%   Calculate f_det(M)
    f_det_M = zeros(1,length(M_vec));
    for i = i_gap:length(M_vec)
        if strcmp(x_str,'M1') && i>i_gap
            f_det_M(i) = trapz(M_vec(i_gap:i),f_det_M1M2(i,i_gap:i));
        elseif strcmp(x_str,'M2') && i<length(M_vec)
            f_det_M(i) = trapz(M_vec(i:end),f_det_M1M2(i:end,i));
        end
    end
%   Calculate f_obs(M)
    if strcmp(x_str,'M1')
        sigma_M = sigma_M1;
    elseif strcmp(x_str,'M2')
        sigma_M = sigma_M1;
    end
    f_det_Mdet = @(Mdet) interp1(M_vec(i_gap:end),f_det_M(i_gap:end),Mdet,'pchip',0);
    f_G_x = @(x) sqrt(2/(pi*sigma_M^2))*(1+erf(1/sqrt(2*sigma_M^2)))^-1*exp(-(x-1).^2/(2*sigma_M^2)).*(x>0);
    f_obs_M = @(Mobs) integral(@(Mdet) f_det_Mdet(Mdet).*f_G_x(Mobs./Mdet)./Mdet,M_vec(i_gap),M_vec(end));
%   Calculate N_obs(M)
    if nargin<7
        Mbins_edges = logspace(log10(M_vec(i_gap)/2),log10(M_vec(end)),Nbins+1);
        Mbins = sqrt(Mbins_edges(1:end-1).*Mbins_edges(2:end));
    else
        Mbins_edges = sqrt(Mbins(1:end-1).*Mbins(2:end));
        Mbins_edges = [Mbins_edges(1)^2/Mbins_edges(2),Mbins_edges,Mbins_edges(end)^2/Mbins_edges(end-1)];
    end
    N_obs_M = zeros(1,Nbins);
    for i=1:Nbins
        N_obs_M(i) = integral(f_obs_M,Mbins_edges(i),Mbins_edges(i+1),'ArrayValued',true);
        if i == 1
            N_obs_M(1) = N_obs_M(1) + integral(f_obs_M,0,Mbins_edges(1),'ArrayValued',true);
        end
        if i == Nbins
            N_obs_M(end) = N_obs_M(end) + integral(f_obs_M,Mbins_edges(end),Inf,'ArrayValued',true);
        end
    end
%   Prepare outputs
    N_obs_x = N_obs_M; xbins = Mbins;
end
%% Calculate N_obs(q)
if strcmp(x_str,'q')
%   Calculate f_det(q,Mtot)
    [M1_mat,M2_mat] = ndgrid(M_vec,M_vec);
    f_det_M1M2_grid = griddedInterpolant(M1_mat(i_gap:end,i_gap:end),M2_mat(i_gap:end,i_gap:end),f_det_M1M2(i_gap:end,i_gap:end),'linear','none');
    q_vec = linspace(M_vec(i_gap)/M_vec(end),1,100);
    Mtot_vec = 2*M_vec;
    [q_mat,Mtot_mat] = ndgrid(q_vec,Mtot_vec);
    M1_mat = Mtot_mat./(q_mat+1);
    M2_mat = q_mat.*Mtot_mat./(q_mat+1);
    f_det_interp = f_det_M1M2_grid(M1_mat,M2_mat); f_det_interp(isnan(f_det_interp))=0;
    f_det_qMtot = f_det_interp.*Mtot_mat./(q_mat+1).^2;
%   Calculate f_det(q)
    f_det_q = trapz(Mtot_vec,f_det_qMtot,2)';
%   Calculate f_obs(q)
    f_det_qdet = @(qdet) interp1(q_vec,f_det_q,qdet,'pchip',0);
    f_G_x = @(x) exp(-(x-1).^2/(2*sigma_q^2)).*(x>0);
    f_obs_q = @(qobs) integral(@(qdet) f_det_qdet(qdet).*f_G_x(qobs./qdet)./qdet,q_vec(1),q_vec(end));
    A = integral(f_det_qdet,q_vec(1),q_vec(end))/integral(f_obs_q,q_vec(1),q_vec(end),'ArrayValued',true);
    f_obs_q = @(qobs) f_obs_q(qobs)*A;
%   Calculate N_obs(q)
    qbins_edges = linspace(0,1,Nbins+1);
    qbins = 0.5*(qbins_edges(1:end-1)+qbins_edges(2:end));
    N_obs_q = zeros(1,Nbins);
    for i=1:Nbins
        N_obs_q(i) = integral(f_obs_q,qbins_edges(i),qbins_edges(i+1),'ArrayValued',true);
    end
%   Prepare outputs
    N_obs_x = N_obs_q; xbins = qbins;
end
%% Calculate N_obs(z)
if strcmp(x_str,'z')
%   Calculate f_det(z)
    f_det_z = reshape(0.5*trapz(M_vec(i_gap:end),trapz(M_vec(i_gap:end),f_det_M1M2z(i_gap:end,i_gap:end,:),2),1),1,length(z_vec));
%   Calculate f_obs(z)
    f_det_zdet = @(zdet) interp1(z_vec,f_det_z,zdet,'pchip',0);
    f_G_x = @(x) sqrt(2/(pi*sigma_z^2))*(1+erf(1/sqrt(2*sigma_z^2)))^-1*exp(-(x-1).^2/(2*sigma_z^2)).*(x>0);
    f_obs_z = @(zobs) integral(@(zdet) f_det_zdet(zdet).*f_G_x(zobs./zdet)./zdet,z_vec(1),z_vec(end));
%   Calculate N_obs(z)
    zbins_edges = logspace(log10(z_vec(2)),log10(z_vec(end)),Nbins+1);
    zbins = sqrt(zbins_edges(1:end-1).*zbins_edges(2:end));
    N_obs_z = zeros(1,Nbins);
    for i=1:Nbins
        N_obs_z(i) = integral(f_obs_z,zbins_edges(i),zbins_edges(i+1),'ArrayValued',true);
        if i == 1
            N_obs_z(1) = N_obs_z(1) + integral(f_obs_z,0,zbins_edges(1),'ArrayValued',true);
        end
        if i == Nbins
            N_obs_z(end) = N_obs_z(end) + integral(f_obs_z,zbins_edges(end),Inf,'ArrayValued',true);
        end
    end
%   Prepare outputs
    N_obs_x = N_obs_z; xbins = zbins;
end