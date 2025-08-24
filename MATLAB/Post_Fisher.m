% This script allows to manipulate the (already calculated) Fisher and
% covariance matrices via fixing, priors and marginalization. The stair
% plot of the manipulated covariance matrix is drawn.
close all; clc;
load([pwd,'/Run 1/Run 1 - Workspace']);
% Set the fixed and marginalized parameters for the Fisher analysis
fixed_params = {'X','beta2','gamma','Mmin','Mmax','lambda','mu_pp','sigma_pp','N_SBH'};
marginalized_params = {};
% Apply priors for the chosen parameters
sigma = Inf*ones(1,length(free_params));
sigma(find(strcmp(free_params,'Mmin'))) = Inf;
F = F + diag(1./(sigma.^2));
% Find the fixed parameters indices within the free parameters list
if ~isempty(fixed_params)
    fixed_ind = size(1,length(fixed_params));
    for i = 1:length(fixed_params)
        for j=1:length(free_params)
            if strcmp(fixed_params(i),free_params(j))
                fixed_ind(i) = j;
            end
        end
    end
%   Omit the fixed parameters indices from the Fisher matrix
    F(fixed_ind,:) = [];
    F(:,fixed_ind) = [];
%   Omit the fixed parameters indices from the following variables
    fiducial_params(fixed_ind) = [];
    scales(fixed_ind) = [];
else
    fixed_ind = [];
end
% Calculate the covariance matrix
cov_matrix = inv(F);
% Keep only the remaining free parameters
free_ind = setdiff(1:length(free_params),fixed_ind);
free_params2 = cell(1,length(free_ind));
for i = 1:length(free_ind)
    free_params2(i) = free_params(free_ind(i));
end
% Find the marginalized parameters indices within the free parameters list
if ~isempty(marginalized_params)
    marginalized_ind = size(1,length(marginalized_params));
    for i = 1:length(marginalized_params)
        for j=1:length(free_params2)
            if strcmp(marginalized_params(i),free_params2(j))
                marginalized_ind(i) = j;
            end
        end
    end
%   Omit the marginalized parameters indices from the covariance matrix
    cov_matrix(marginalized_ind,:) = [];
    cov_matrix(:,marginalized_ind) = [];
%   Omit the marginalized parameters indices from the following variables
    fiducial_params(marginalized_ind) = [];
    scales(marginalized_ind) = [];
else
    marginalized_ind = [];
end
% Keep only the remaining free parameters
free_ind = setdiff(1:length(free_params2),marginalized_ind);
free_params3 = cell(1,length(free_ind));
for i = 1:length(free_ind)
    free_params3(i) = free_params2(free_ind(i));
end
% Rename the free parameters
names = NameOrganizer(free_params3);
% Make the stair plot
out = MakeStairPlot(cov_matrix,fiducial_params./scales,names,['aLIGO: ',num2str(T),'yrs'],scales);
annotation('textbox',[0.4,0.5,0.5,0.5],'String',['Number of total events: ',num2str(floor(sum(N_obs_0)))],...
    'FitBoxToText','on','EdgeColor','none','FontName','Times New Roman','FontSize',18);