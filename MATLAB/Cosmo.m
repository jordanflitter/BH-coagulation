% This functions calculates several cosmological quantities as described 
% below.
% 
% The function receives the following inputs:
% M_vec:    a vector of masses (in units of solar masses).
% z_vec:    a vector of redshift values.
% params:   a dictionary containing numerical parameters. See documantation
%           at Main.m for more details.
% 
% The function returns the following outputs:
% A_dyn:    a 3D array of A_dyn(z)=4*pi*c*chi(z)^2*nGC(z)/((1+z)*H(z)).
%           In each page, the rows and columns corrspond to different mass 
%           values (which are specified by the vector M_vec), while 
%           different pages correspond to different redshift samples (which 
%           are specified by the vector z_vec).
% A_stat:   a 3D array of A_stat(z)=4*pi*c*chi(z)^2*(1+z)^3/((1+z)*H(z)) 
%           (in units of Gpc^3). In each page, the rows and columns 
%           corrspond to different mass values (which are specified by the 
%           vector M_vec), while different pages correspond to different 
%           redshift samples (which are specified by the vector z_vec).
% z_t:      a vector of redshift values which correspond to time samples
%           which areassociated with Gamma_i,j(t).
% 
function [A_dyn,A_stat,z_t] = Cosmo(M_vec,z_vec,params)
% Extract parameters from the dictionary
Nt = params('Nt');
t_c = params('t_c');
% Set cosmological parameters
Omega_r0 = 0;
Omega_k0 = 0;
Omega_Lambda0 = 0.685;
Omega_m0 = 0.315;
h = 0.674;
c = 3e+5; %Km/sec
H0 = 100*h; %Km*Mpc^-1*sec^-1
nGC0 = 3; %Mpc^-3
% Compute H(z)
H = H0*sqrt(Omega_r0*(1+z_vec).^4 + Omega_m0*(1+z_vec).^3 + Omega_k0*(1+z_vec).^2 + Omega_Lambda0);
% Compute nGC(z)
nGC = nGC0*(H/H0).^3;
% Compute chi(z)
chi = zeros(1,length(z_vec));
for i = 1:length(z_vec)
    if i>1
        z_tag = z_vec(1:i);
        chi(i) = c*trapz(z_tag,1./H(1:i)); %Mpc
    end
end
% Compute A_dyn(z) and A_stat(z)
[~,~,A_dyn] = ndgrid(M_vec,M_vec,4*pi*c*chi.^2.*nGC./((1+z_vec).*H));
[~,~,A_stat] = ndgrid(M_vec,M_vec,4*pi*c*chi.^2.*((1+z_vec).^3)./((1+z_vec).*H)*1e-9); %Gpc^3
% Compute t(z)
t = zeros(1,length(z_vec));
for i = 1:(length(z_vec))
    if i < length(z_vec)
        z_tag = z_vec(i:end);
        t(i) = 978.56*trapz(z_tag,1./((1+z_tag).*H(i:end))); %Gyears
    end
    t0 = t(1); %Gyears
end
% Compute z(t), where t are the time samples associated with Gamma_i,j(t)
z_grid = griddedInterpolant(t(end:-1:1),z_vec(end:-1:1));
z_t = z_grid((0:(1/(Nt-1)):1)*t_c + (t0-t_c));