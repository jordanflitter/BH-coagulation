% This functions plots the events distributions with respect to either M1,
% M2, q, or z.
% 
% The function receives the following inputs:
% xbins:        a vector of the bins locations.
% N_obs_x_stat: a vector of the events distribution which includes 
%               the contribution from only the static channels.
% N_obs_x:      a matrix of the events distributions which include 
%               contributions from both static and dynamic channels. Each 
%               row corresponds to a different tested scenario.
% x_str:        a string which specifies the corresponding distribution. 
%               The options are 'M1', 'M2', 'q' or 'z'.
% LIGO_RUN:     a string which specifies the corresponding experiment.
% 
% The function returns the following outputs:
% F:            The figure handle of the plotted distribution.
% 
function F = Nobs_plotter(xbins,N_obs_x0,N_obs_x,x_str,LIGO_RUN)
F = figure;
% Recalculate the bins edges
if strcmp(x_str,'M1') || strcmp(x_str,'M2') || strcmp(x_str,'z')
    xbins_edges = sqrt(xbins(1:end-1).*xbins(2:end));
    xbins_edges = [xbins_edges(1)^2/xbins_edges(2),xbins_edges,xbins_edges(end)^2/xbins_edges(end-1)];
elseif strcmp(x_str,'q')
    xbins_edges = 0.5*(xbins(1:end-1) + xbins(2:end));
    xbins_edges = [2*xbins_edges(1)-xbins_edges(2),xbins_edges,2*xbins_edges(end)-xbins_edges(end-1)];
end
% For O3a only, place the data in bins 
if strcmp(LIGO_RUN,'O3a')
    load('O3a_distributions');
    if strcmp(x_str,'M1')
        f_x = f_M1;
        x_vec = M1_data_vec;
    elseif strcmp(x_str,'M2')
        f_x = f_M2;
        x_vec = M2_data_vec;
    elseif strcmp(x_str,'q')
        f_x = f_q;
        x_vec = q_data_vec;
    elseif strcmp(x_str,'z')
        f_x = f_z;
        x_vec = z_data_vec;
    end
    f_x_x = @(x) interp1(x_vec,f_x,x,'linear',0);
    N_obs_data = zeros(1,length(xbins));
    for i=1:length(xbins)
        N_obs_data(i) = integral(f_x_x,xbins_edges(i),xbins_edges(i+1),'ArrayValued',true);
        if i == 1
            N_obs_data(1) = N_obs_data(1) + integral(f_x_x,0,xbins_edges(1),'ArrayValued',true);
        end
        if i == length(xbins)
            N_obs_data(end) = N_obs_data(end) + integral(f_x_x,xbins_edges(end),Inf,'ArrayValued',true);
        end
    end
    % Plot the histogram for O3a data
    for i = 1:length(xbins)
        area([xbins_edges(i),xbins_edges(i+1)],[N_obs_data(i),N_obs_data(i)],'FaceColor',[0.7,0.7,0.7]); hold on;
    end
end
% Plot the curve of the static channel contribution
P = plot(xbins,N_obs_x0,'k','LineWidth',2); hold on;
set(gca,'ColorOrderIndex',1);
% Plot error bars (for O5 only)
if strcmp(LIGO_RUN,'O5')
    rnd = poissrnd(N_obs_x0);
    for j=1:length(xbins)
        while rnd(j)==1
            rnd(j) = poissrnd(N_obs_x0(j));
        end
    end
    errorbar(xbins,rnd,sqrt(rnd),'Color',P.Color,'LineWidth',2,'LineStyle','none');
end
% Plot the curves for the different scenarios and write the total number of
% events
annotation('textbox',[0.15,0.4,0.5,0.5],'String','# events: ',...
    'FitBoxToText','on','EdgeColor','none','FontName','Times New Roman','FontSize',18);
for i = 1:size(N_obs_x,1)
    P = plot(xbins,N_obs_x(i,:),'LineWidth',2);
    annotation('textbox',[0.32+0.07*(i-1),0.4,0.5,0.5],'String',num2str(round(sum(N_obs_x(i,:)))),...
        'FitBoxToText','on','EdgeColor','none','FontName','Times New Roman','FontSize',18,'Color',P.Color);
%   Plot error bars (for O5 only)
    if strcmp(LIGO_RUN,'O5')
        rnd = poissrnd(N_obs_x(i,:));
        for j=1:length(xbins)
            while rnd(j)==1
                rnd(j) = poissrnd(N_obs_x(i,j));
            end
        end
        errorbar(xbins,rnd,sqrt(rnd),'Color',P.Color,'LineWidth',2,'LineStyle','none');
    end
end
% Adjust figures' boundaries, scale, labels, title...
if strcmp(LIGO_RUN,'O3a')
    N_obs_TH = 0.1;
    ylim_min = 0;
    ylim_factor = 1.2;
    if strcmp(x_str,'M1') || strcmp(x_str,'M2')
        plot([0.1,300],[1,1],'k--');
    elseif strcmp(x_str,'q')
        plot([0,1],[1,1],'k--');
    elseif strcmp(x_str,'z')
        plot([1e-4,10],[1,1],'k--');
    end
    max_N = max([max(N_obs_data), max(N_obs_x0)]);
    min_x = min([xbins(find(N_obs_data>N_obs_TH,1,'first')),xbins(find(N_obs_x0>N_obs_TH,1,'first'))]);
    max_x = max([xbins(find(N_obs_data>N_obs_TH,1,'last')),xbins(find(N_obs_x0>N_obs_TH,1,'last'))]);
elseif strcmp(LIGO_RUN,'O5')
    N_obs_TH = 1;
    ylim_min = 1;
    ylim_factor = 2;
    set(gca,'YScale','log');
    set(gca,'YTick',[1,10,100,1000,10000],'YTickLabel',{'1','10','100','1000','10000'});
    max_N = max(N_obs_x0);
    min_x = xbins(find(N_obs_x0>N_obs_TH,1,'first'));
    max_x = xbins(find(N_obs_x0>N_obs_TH,1,'last'));
end
for i = 1:size(N_obs_x,1)
    max_N = max([max_N, N_obs_x(i,:)]);
    min_x = min([min_x,xbins(find(N_obs_x(i,:)>N_obs_TH,1,'first'))]);
    max_x = max([max_x,xbins(find(N_obs_x(i,:)>N_obs_TH,1,'last'))]);
end
xlim([0.5*min_x,1.5*max_x]);
ylim([ylim_min,ylim_factor*max_N]);
set(gca,'Layer','top','FontName','Times New Roman','FontSize',20);
ylabel('$N_{\mathrm{obs}}$','interpreter','Latex','FontSize',24);
if strcmp(x_str,'M1')
    set(gca,'XScale','log');
    set(gca,'XTick',[1,2,3,5,10,30,50,100:100:xbins_edges(end)]);
    xlabel('$M_1\,[M_{\odot}]$','interpreter','Latex','FontSize',24);
    title(['Heavier Mass Distribution (',LIGO_RUN,')'],'interpreter','Latex','FontSize',20);
elseif strcmp(x_str,'M2')
    set(gca,'XScale','log');
    set(gca,'XTick',[1,2,3,5,10,30,50,100:100:xbins_edges(end)]);
    xlabel('$M_2\,[M_{\odot}]$','interpreter','Latex','FontSize',24);
    title(['Lighter Mass Distribution (',LIGO_RUN,')'],'interpreter','Latex','FontSize',20);
elseif strcmp(x_str,'q')
    xlim([0,1]);
    set(gca,'XTick',[0,0.25,0.5,0.75,1]);
    xlabel('$q$','interpreter','Latex','FontSize',24);
    title(['Mass Ratio Distribution (',LIGO_RUN,')'],'interpreter','Latex','FontSize',20);
elseif strcmp(x_str,'z')
    set(gca,'XScale','log');
    set(gca,'XTick',[0.01,0.05,0.1,0.2,0.5,1,2,3]);
    xlabel('$z$','interpreter','Latex','FontSize',24);
    title(['Redshift Distribution (',LIGO_RUN,')'],'interpreter','Latex','FontSize',20);
end