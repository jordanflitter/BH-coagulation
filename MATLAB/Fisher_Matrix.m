% This function calculates the Fisher matrix.
% 
% The function receives the following inputs:
% fiducial_params:  a vector containing the free parameters' fiducial 
%                   values.
% N_obs_matrix:     a matrix of the number of observed events. Each row
%                   corresponds to a different variation in the parameters
%                   while each column corresponds to a different mass bin.
% delta:            The ratio between the variation in the parameters and
%                   their fiducial values.
% 
% The function returns the following outputs:
% F:        The Fisher matrix.
% scales:   a vector which specifies the rescaling which was done on each
%           parameter.
% 
function [F,scales] = Fisher_Matrix(fiducial_params,N_obs_matrix,delta)
% Extract various quantities from the inputs
Npar = length(fiducial_params);
Nbins = size(N_obs_matrix,2);
% N0 is the model's prediction for the fiducial values
N0 = N_obs_matrix(Npar+1,:);
% Set initial conditions for the loop
F = 0; scales = ones(1,Npar); tries = 0; norm_ind = [];
% In each iteraion, the Fisher matrix is recalculated. If the matrix is
% close to be singular (or misbehaves), one of the paremeters is rescaled
while ((abs(det(F))<0.01 || abs(det(F))>100) && tries<Npar)
%   Calculate the parameters' variation and the derivatives with respect to
%   the parameters
    dpar = (fiducial_params)*delta+(fiducial_params==0)*delta;
    dNdpar = zeros(Npar,Nbins);
    for ii=1:Npar
        N_p = N_obs_matrix(ii,:);
        N_m = N_obs_matrix(Npar+1+ii,:);
        dNdpar(ii,:) = (N_p-N_m)/(2*dpar(ii));
    end
    % Calculate Fisher matrix
    F = zeros(Npar,Npar);
    for ii=1:Npar
        for jj=1:Npar
            d1 = dNdpar(ii,:);
            d2 = dNdpar(jj,:);
            F(ii,jj) = sum(d1.*d2./N0);
        end
    end
%   If necessary, rescale parameters to order of unity so the Fisher matrix
%   is well-behaved
    if (abs(det(F))<0.01 || abs(det(F))>100)
        max_norm = max(abs(log10(fiducial_params(fiducial_params~=0))));
        ind = find(abs(log10(fiducial_params))==max_norm);
        if ~ismember(ind,norm_ind)
            norms = 10.^floor(log10(abs(fiducial_params(ind))));
            fiducial_params(ind) = fiducial_params(ind)/norms;
            scales(ind)=norms;
            norm_ind = [norm_ind,ind];
        end
        tries = tries + 1;
    else
        break;
    end
end