% This is the main script for simulating the merger events distributions
% which are detected by aLIGO. 2 contributions are considered:
%   1. Static channel: a channel where the Black Hole Mass Function (BHMF) 
%      does not change with time - it is simply the Initial Mass Function 
%      (IMF).
%   2. Dynamic channel: a channel where the BHMF evolves with time. It is 
%      assumed that mergers from this channel come from Globular Clusters
%      (GCs). Only one cluster is simulated, which is considered to be the 
%      typical observable GC by aLIGO. The GC BHMF is evolved via the
%      coagulation equation.
%
% The script can be used to:
%   1. Plot the evolution of the BHMF of the dynamic channel.
%   2. Plot the number of BHs in the dynamic channel vs. time.
%   3. Plot the total BHs mass in the dynamic channel vs. time.
%   4. Plot the total merger rate (of the dynamic channel) vs. time.
%   5. Plot the the detected merger events 1D distributions with respect to 
%      various variables: heavier mass (M1), smaller mass (M2), mass ratio
%      (q), redshift (z).
%   6. Compare between different scenarios. The comparison is done by
%      plotting the BHMFs evolution, mass concentration (eta) vs. time,
%      probability to form Intermediate Mass Black Hole (IMBH) vs. time,
%      total merger rate vs. time, and the detected events distributions.
%   7. Plot the parameters' confidence ellipses from Fisher analysis.
% 
% The script is capable of simulating several different scenarios and save
% all the figures and workspace if specified.
% 
% There are several parameters which are needed to be determined. Some of
% them are model independent and set the simulation accuracy and runtime,
% while some of them are model dependent. Below are listed all the 
% parameters (model dependent and independent) and flags which are used
% in this script.
% 
% M0:     The black hole mass resolution of the simulation (in units of
%         solar masses).
% i_max:  The cutoff for the infinite sum at the coagulation equation. This
%         parameter sets the amount of mass values which participate at 
%         the coagulation process. This value becomes relvant only if
%         SYME_TYPE="Test", otherwise this value will be calculated
%         according to the cluster's total mass.
% I:      Total number of iterations in the coagulation process.
% Nt:     Number of time samples in which the rate matrix shall be saved. 
%         Maximal value is I.
% Nbins:  Number of bins to be shown in the events distributions.
% delta:  The ratio between the variation in the parameters and their 
%         fiducial values (for the Fisher analysis).
% 
% observation_flag: The detected merger events distributions shall be 
%                   calculated only if this flag is true.
% plot_coagulation: The evolution of the BHMF (of the dynamic channel), the 
%                   BHMF's moments, and the total merger rate shall be 
%                   plotted (for the fiducial parameters) only if this 
%                   flag is true.
% plot_M1_dist:     The heavier mass distribution (at detector frame) shall 
%                   be plotted (for the fiducial parameters) only if this 
%                   flag is true.
% plot_M2_dist:     The lighter mass distribution (at detector frame) shall 
%                   be plotted (for the fiducial parameters) only if this 
%                   flag is true.
% plot_q_dist:      The mass ratio distribution (at detector frame) shall 
%                   be plotted (for the fiducial parameters) only if this 
%                   flag is true.
% plot_z_dist:      The redshift distribution (at detector frame) shall 
%                   be plotted (for the fiducial parameters) only if this 
%                   flag is true.
% save_workspace:   The workspace shall be saved only if this flag is true.
% save_figures:     figures shall be saved only if this flag is true.
% num_of_scenarios: The number of scenarios to be simulated. If the list of
%                   free_params is not empty, than this parameter is set to
%                   1.
% Comparison_flag:  Comparison between the different simulated scenarios
%                   shall be performed only if this flag is true, the list
%                   of free_params is empty, and num_of_scenarios is
%                   greater than 1.
% free_params:      This is a list of strings which specify the free
%                   parameters to be varied in the Fisher analysis. If this
%                   list is empty then no Fisher analysis shall be
%                   preformed.
% 
% SYM_TYPE is a string that controls the runtime of the program. If this
%          string is "Test", then the simulation runtime will be reduced
%          significantly according to the value of i_max.
% IMF_TYPE is a string which determines the type of IMF to be used in the
%          simulation. There are 6 options:
%   'Seed':     Corresponds to a "seed" distribution where all the initial
%               BHs start with the same mass M0.
%   'Model':    Corresponds to Salpeter power law which has both lower and 
%               upper mass cutoffs. This power law is then weighted with a
%               Gaussian to account for the contribution of Pulsational
%               Pair Instability Supernovae (PPSNe) to the IMF. Neutron
%               Stars (NSs) might also be considered and their contribution
%               to the IMF also comes in the form of a Gaussian. The
%               contribution of NSs to the BHMF might enter the simulation 
%               at t=0 or later (according to the value ot t_NS).
% LIGO_RUN is a string which specifies the corresponding experiment. There 
%          are 2 options:
%  O3a: corresponds to O3 noise curve and total observation time of 177.3 days.
%  O5:  corresponds to O5 noise curve and total observation time of 6 years.
% 
% alpha:    The Salpeter power law index for the SBHs IMF.
% Mmin:     The lower mass cutoff for the SBHs IMF (in units of solar 
%           masses).
% Mmin:     The upper mass cutoff for the SBHs IMF (in units of solar 
%           masses).
% mu_pp:    The PPSNe Gaussian mean (in units of solar masses).
% sigma_pp: The PPSNe Gaussian statandard deveiation (in units of solar
%           masses).
% lambda:   The proportional contribution of the PPSNe to the IMF.
% mu_NS:    The NSs Gaussian mean (in units of solar masses).
% sigma_NS: The NSs Gaussian statandard deveiation (in units of solar
%           masses).
% 
% X:       The proportional contribution of the dynamic channel to the 
%          detected events distributions. If X=0 then the coagulation 
%          equation shall not be solved (as only the static channel is 
%          relevant for this case).
% R_0:     The merger rate density for the static channel (in units of 
%          Gpc^-3*year^-1).
% R11:     The normalized merger rate for 2 BHs of mass M0 in the dynamic 
%          channel (in units of year^-1).
% beta:    The power index of the quotient term in the rate kernel of the 
%          dynamic channel.
% beta2:   The power index of the quotient term in the rate kernel of the 
%          static channel.
% gamma:   The power index of the sum term in the rate kernel of the 
%          dynamic/static channel.
% t_c:     Overall coagulation time (in units of Gyears).
% f_loss:  The fracation of mass loss in equal masses mergers due to energy 
%          carried away by Gravitational Waves (GWs).
% M_esc:   a parameter which is associated with the ejections rate. The
%          smaller this parameter is the more ejections (in units of solar
%          masses).
% N_BH:    The initial number of BHs in the cluster.
% N_NS:    The initial number of NSs in the cluster.
% t_NS:    The time when NSs enter the simulation (in units of Gyears).
% 
close all; clc; clear;
%% Intialization
% Make a dictionary for numerical parameters 
% Set the values of the simulation parameters
params1 = containers.Map('KeyType','char','ValueType','double');
params1('M0') = 1; %solar masses
params1('i_max') = 1000;
params1('I') = 300;
params1('Nt') = 30;
params1('Nbins') = 25;
params1('delta') = 0.1;
% Set the values of the flag parameters
params1('observation_flag') = 0;
params1('plot_coagulation') = 1;
params1('plot_M1_dist') = 0;
params1('plot_M2_dist') = 0;
params1('plot_q_dist') = 0;
params1('plot_z_dist') = 0;
params1('save_workspace') = 1;
params1('save_figures') = 1;
% Set the number of scenarios to be simulated and decide whether to make
% comparison between the simulated scenarios
num_of_scenarios = 1;
Comparison_flag = 1;
% Define the free parameters for the Fisher analysis
free_params = {};
% Make a dictionary for the string parameters and set their values
str_params1 = containers.Map('KeyType','char','ValueType','char');
str_params1('SYM_TYPE') = 'Test';
str_params1('IMF_TYPE') = 'Model';
str_params1('LIGO_RUN') = 'O3a';
% Set the values of the IMF parameters
params1('alpha') = 2.35;
params1('Mmin') = 5; %solar masses
params1('Mmax') = 50; %solar masses
params1('mu_pp') = 35; %solar masses
params1('sigma_pp') = 3; %solar masses
params1('lambda') = 0.08;
params1('mu_NS') = 1.33; %solar masses
params1('sigma_NS') = 0.09; %solar masses
% Set the values of the other parameters
params1('X') = 1;
params1('R_0') = 129; %Gpc^-3*year^-1
params1('R11') = 6.1e-13; %year^-1
params1('beta') = 2;
params1('beta2') = 2;
params1('gamma') = 0;
params1('t_c') = 13; %Gyear
params1('f_loss') = 0;
params1('M_esc') = Inf; %solar masses
params1('v_esc') = 50; %Km/sec
params1('N_SBH') = 300;
params1('N_NS') = 0;
params1('t_NS') = 0; %Gyear
% Auxiliary code lines
if Comparison_flag
   params1('save_workspace') = 1;
end
Npar = length(free_params);
if Npar~=0
    params1('observation_flag') = 1;
    num_of_scenarios = 1;
    Comparison_flag = 0;
end
%% Preparation
% COPY & PASTE THE 3 FOLLOWING LINES AND CHANGE THEM ACCORDING TO THE
% SCENARIOS TO BE SIMULATED!
params2 = containers.Map(params1.keys,params1.values);
str_params2 = containers.Map(str_params1.keys,str_params1.values);
params2('M_esc') = 5;
params2('R11') = 1.15e-12;
params3 = containers.Map(params1.keys,params1.values);
str_params3 = containers.Map(str_params1.keys,str_params1.values);
params3('v_esc') = 50;
params3('R11') = 1.3e-12;
params4 = containers.Map(params1.keys,params1.values);
str_params4 = containers.Map(str_params1.keys,str_params1.values);
params4('X') = 0.5;
params5 = containers.Map(params1.keys,params1.values);
str_params5 = containers.Map(str_params1.keys,str_params1.values);
params5('v_esc') = 50;
params6 = containers.Map(params1.keys,params1.values);
str_params6 = containers.Map(str_params1.keys,str_params1.values);
params6('gamma') = 2;
params6('beta') = 0;
params6('R11') = 2.75e-16;
params6('M_esc') = 5;
params7 = containers.Map(params1.keys,params1.values);
str_params7 = containers.Map(str_params1.keys,str_params1.values);
params7('gamma') = 2;
params7('beta') = 0;
params7('R11') = 2.75e-16;
params7('M_esc') = 5;
params7('f_loss') = 0.05;
% Prepare super dictionaries
super_params = {params1,params2,params3,params4,params5,params6,params7};
super_str_params = {str_params1,str_params2,str_params3,str_params4,str_params5,str_params6,str_params7};
super_parmas = super_params(1:num_of_scenarios); super_str_parmas = super_str_params(1:num_of_scenarios);
%% Calculation
for i = 1:num_of_scenarios
    if Npar == 0
        disp(['Now simulating scenario ',num2str(i),'...']);
    else
        disp('Now simulating for the fiducial parameters...');
    end
    tic;
    [f_det_M1M2z,M_vec,z_vec] = BHMF(super_params{i},super_str_params{i},i);
    if i == 1
        Elapsed_Time = toc;
    else
        Elapsed_Time = Elapsed_Time + toc;
    end
    if Npar == 0
        Display_Timer('Comparison',Elapsed_Time,num_of_scenarios,i,super_params);
    end
    if ~Comparison_flag && (Npar == 0)
        fprintf('\n');
    end
end
%% Comparison
if Comparison_flag && (num_of_scenarios > 1)
    disp('Now plotting the comparison figures...');
    Comparison(num_of_scenarios,params1);
    fprintf('\n');
end
%% Fisher analysis
if Npar~=0
   Fisher_Analysis(f_det_M1M2z,M_vec,z_vec,free_params,params1,str_params1,Elapsed_Time);
end