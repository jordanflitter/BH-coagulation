% This function plots comparison figures for preveiously simulated
% scenarios. The function assumes that all the simulations' data is saved
% within the current directory at folders named "Run 1", "Run 2", etc...
% The function plots the following figures:
% 1. Final BHMFs (of the dynamic channel).
% 2. number of BHs (of the dynamic channel) vs. time.
% 3. total BHs mass (of the dynamic channel)vs. time.
% 4. probability to form Intermediate Mass Black Hole (IMBH) vs. time.
% 5. total merger rate vs. time.
% 6. heavier mass distribution.
% 7. lighter mass distribution.
% 8. mass ratio distribution.
% 9. redshift distribution.
% 
% The function receives the following inputs:
% num_of_scenarios: number of previously simulated scenarios to be
%                   compared.
% params:   a dictionary containing numerical parameters. See documantation
%           at Main.m for more details.
%
function Comparison(num_of_scenarios,params)
% Extract the parameters from the dictionaries
Nbins = params('Nbins');
X = params('X');
observation_flag = params('observation_flag');
plot_coagulation = params('plot_coagulation');
plot_M1_dist = params('plot_M1_dist');
plot_M2_dist = params('plot_M2_dist');
plot_q_dist = params('plot_q_dist');
plot_z_dist = params('plot_z_dist');
save_figures = params('save_figures');
% Close all open figures and create a new comaprison folder
close all;
if ~exist('Comparison', 'dir')
    mkdir('Comparison');
end
if plot_coagulation && X~=0
%   Final BHMFs figure
    for i = 1:num_of_scenarios
        load([pwd,'/Run ',num2str(i),'/Run ',num2str(i),' - Workspace']);
        if i == 1
            F1 = figure; H = area(M,f_0/M0); set(gca,'Layer','top');
            H.FaceColor = [0.7,0.7,0.7]; H.EdgeColor = [1,1,1];
            hold on; plot(M,f_0/M0,'k'); set(gca,'ColorOrderIndex',1);
            Mgap = M(find(f_src>1,1)); Mmax = M(end);
            if ~isempty(Mgap)
                xlim([0.6*Mgap,1.2*Mmax]);
            end
            min_f = 1e-10*N_SBH/M0; max_f = max(f_0/M0);
        end
        plot(M(1:length(f_src)),f_src,'LineWidth',2);
        max_f = max([max_f, f_src]);
    end
    ylim([min_f,2*max_f]);
    set(gca,'XScale','log','YScale','log');
    set(gca,'XTick',[1,2,3,5,10,30,50,100,200,1000:1000:M(end)]);
    set(gca,'FontName','Times New Roman','FontSize',20);
    xlabel('$M\,[M_{\odot}]$','interpreter','Latex','FontSize',24);
    ylabel('$f_{\mathrm{BH}}\,[M_\odot^{-1}]$','interpreter','Latex','FontSize',24);
    title(['Black Hole Mass Function, $t$ = ',num2str(t_c),' Gyears'],'interpreter','Latex','FontSize',18);
    savefig(F1,[pwd,'/Comparison/Comparison - Coagulation']);
%   eta, P_IMBH, and total merger rate figures
    F2 = figure;
    F3 = figure;
    F4 = figure;
    F5 = figure;
    for i = 1:num_of_scenarios
        load([pwd,'/Run ',num2str(i),'/Run ',num2str(i),' - Workspace']);
        Ntot = moments(1,:);
        Mtot = moments(2,:);
        digits(100);
        P_IMBH = double(1-power(1-vpa(P_500),Ntot)); digits(32);
        figure(F2.Number); plot(t_c*(0:1/I:1),Ntot,'LineWidth',2); hold on;
        figure(F3.Number); plot(t_c*(0:1/I:1),Mtot,'LineWidth',2); hold on;
        figure(F4.Number); plot(t_c*(0:1/I:1),P_IMBH,'LineWidth',2); hold on;
        figure(F5.Number); plot(t_c*(0:(1/(Nt-1)):1),Gamma_tot,'LineWidth',2); hold on;
    end
    figure(F2.Number);
    xlim([0,t_c]);
    set(gca,'FontName','Times New Roman','FontSize',20);
    xlabel('Time [Gyears]','interpreter','Latex','FontSize',24);
    ylabel('$N_\mathrm{tot}$','interpreter','Latex','FontSize',24);
    title('Number of Black Holes Vs. Time','interpreter','Latex','FontSize',20);
    savefig(F2,[pwd,'/Comparison/Comparison - N_tot']);
    figure(F3.Number);
    xlim([0,t_c]);
    set(gca,'FontName','Times New Roman','FontSize',20);
    xlabel('Time [Gyears]','interpreter','Latex','FontSize',24);
    ylabel('$M_\mathrm{tot}\,[M_\odot]$','interpreter','Latex','FontSize',24);
    title('Total Black Holes Mass Vs. Time','interpreter','Latex','FontSize',20);
    savefig(F3,[pwd,'/Comparison/Comparison - M_tot']);
    figure(F4.Number); plot([0,t_c],[0.01,0.01],'k-.');
    set(gca,'YScale','log');
    xlim([0,t_c]); ylim([1e-40,1]);
    set(gca,'FontName','Times New Roman','FontSize',20);
    xlabel('Time [Gyears]','interpreter','Latex','FontSize',24);
    ylabel('$P_{\mathrm{IMBH}}$','interpreter','Latex','FontSize',24);
    title('Probability for IMBH formation','interpreter','Latex','FontSize',20);
    savefig(F4,[pwd,'/Comparison/Comparison - IMBH Probability']);
    figure(F5.Number);
    set(gca,'YScale','log')
    xlim([0,t_c]);
    set(gca,'FontName','Times New Roman','FontSize',20);
    xlabel('Time [Gyears]','interpreter','Latex','FontSize',24);
    ylabel('$\Gamma_{\mathrm{tot}}\,[\mathrm{year}^{-1}]$','interpreter','Latex','FontSize',24);
    title('Cluster''s total merger rate Vs. Time','interpreter','Latex','FontSize',20);
    savefig(F5,[pwd,'/Comparison/Comparison - Total Merger Rate']);
end
% M1 distribution figure
if observation_flag && plot_M1_dist
    N_obs_x = zeros(num_of_scenarios,Nbins);
    for i = 1:num_of_scenarios
        load([pwd,'/Run ',num2str(i),'/Run ',num2str(i),' - Workspace']);
        if i == 1
            N_obs_x_stat = N_obs_M1_stat;
            xbins = Mbins;
        end
        N_obs_x(i,:) = N_obs_M1;
    end
    F6 = Nobs_plotter(xbins,N_obs_x_stat,N_obs_x,'M1',LIGO_RUN);
    if save_figures
        savefig(F6,[pwd,'/Comparison/Comparison - M1 Distribution']);
    end
end
% M2 distribution figure
if observation_flag && plot_M2_dist
    N_obs_x = zeros(num_of_scenarios,Nbins);
    for i = 1:num_of_scenarios
        load([pwd,'/Run ',num2str(i),'/Run ',num2str(i),' - Workspace']);
        if i == 1
            N_obs_x_stat = N_obs_M2_stat;
            xbins = Mbins;
        end
        N_obs_x(i,:) = N_obs_M2;
    end
    F7 = Nobs_plotter(xbins,N_obs_x_stat,N_obs_x,'M2',LIGO_RUN);
    if save_figures
        savefig(F7,[pwd,'/Comparison/Comparison - M2 Distribution']);
    end
end
% q distribution figure
if observation_flag && plot_q_dist
    N_obs_x = zeros(num_of_scenarios,Nbins);
    for i = 1:num_of_scenarios
        load([pwd,'/Run ',num2str(i),'/Run ',num2str(i),' - Workspace']);
        if i == 1
            N_obs_x_stat = N_obs_q_stat;
            xbins = qbins;
        end
        N_obs_x(i,:) = N_obs_q;
    end
    F8 = Nobs_plotter(xbins,N_obs_x_stat,N_obs_x,'q',LIGO_RUN);
    if save_figures
        savefig(F8,[pwd,'/Comparison/Comparison - q Distribution']);
    end
end
% z distribution figure
if observation_flag && plot_z_dist
    N_obs_x = zeros(num_of_scenarios,Nbins);
    for i = 1:num_of_scenarios
        load([pwd,'/Run ',num2str(i),'/Run ',num2str(i),' - Workspace']);
        if i == 1
            N_obs_x_stat = N_obs_z_stat;
            xbins = zbins;
        end
        N_obs_x(i,:) = N_obs_z;
    end
    F9 = Nobs_plotter(xbins,N_obs_x_stat,N_obs_x,'z',LIGO_RUN);
    if save_figures
        savefig(F9,[pwd,'/Comparison/Comparison - z Distribution']);
    end
end