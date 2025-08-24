% This function returns the detected merger events 2D distribution with 
% respect to (M1,M2).
% The function receives the following inputs:
% 
% f_det_M1M2:   a matrix of the events density 2D function with respect 
%               to M1, M2. The rows and the columns of this matrix 
%               correspond to the vector M_vec.
% ibins_edges:  a vector the mass bins edges (in units of M0).
% N_ind:        an optional input. It's a vector of the indices of bins 
%               where the number of events is greater than 10 (with respect
%               to the fiducial values) and M1>=M2.
% 
% The function returns the following outputs:
% N_obs_vec:    a vector containing the number of detected events at each 
%               bin (only if this number is greater than 10 and M1>=M2).
% N_ind:        an optional output. It's a vector of the indices of bins 
%               where the number of events is greater than 10 (with respect
%               to the fiducial values) and M1>=M2.
% 
function [N_obs_vec,N_ind] = Binning_2D(f_det_M1M2,ibins_edges,N_ind)
% Bin the 2D density function f_det(M1,M2)
Nbins = length(ibins_edges)-1;
N_obs_M1M2 = zeros(Nbins,Nbins);
for i = 1:Nbins
    for j = 1:Nbins
        i_vec = ibins_edges(i):ibins_edges(i+1);
        j_vec = ibins_edges(j):ibins_edges(j+1);
        N_obs_M1M2(i,j) = trapz(i_vec,trapz(j_vec,f_det_M1M2(i_vec,j_vec),2),1);
    end
end
% Include only bins where M1>=M2. Divide by 2 for bins where M1=M2
N_obs_M1M2 = tril(N_obs_M1M2.*(1-0.5*eye(Nbins)))+triu(nan*ones(Nbins),1);
% Reshape the matrix into a vector
N_obs_vec = reshape(N_obs_M1M2,1,Nbins^2);
% Return the indices of the bins where the number of events is greater than
% 10 (with respect to the fiducial values case) and M1>=M2
if nargin < 3
    N_obs_vec(N_obs_vec<10) = nan;
    N_ind = ~isnan(N_obs_vec);
end
N_obs_vec = N_obs_vec(N_ind);