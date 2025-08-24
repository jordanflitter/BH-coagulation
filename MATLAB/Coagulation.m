% This function solves the coagulation equation in the context of BHs and
% NSs in a cluster. This function also plots the BHMF evolution, the number
% of BHs vs. time, the total BHs mass vs. time, and the cluster's total
% merger rate vs. time.
% 
% The function receives the following inputs:
% M:          a vector containing linearly spaced mass values
%            (in units of solar masses).
% N0:         a vector containing samples of the IMF at specific mass
%             values which are given by the vector M (in units of 
%             1/solar masses).
% params:   a dictionary containing numerical parameters. See documantation
%           at Main.m for more details.
% str_params: a dictionary containing string parameters. See documantation
%             at Main.m for more details.
% 
% The function returns the following outputs:
% f_src:      a vector containing samples of the BHMF by the end of the 
%             coagulation process at specific mass values which are given 
%             by the vector M (in units of 1/solar masses).
% Gamma:      a 3D array of the merger rates (in units of year^-1).
%             In each page, the rows and columns corrspond to different 
%             mass values (which are given by the vector M), while 
%             different pages correspond to different time samples (at 
%             each t_c/Nt).
% Gamma_tot:  a vector containing the overall merger rate (in units of 
%             year^-1). Each cell corresponds to different time sample (at
%             each t_c/I).
% moments:    a matrix where its 3 rows correspond to the BHMF's 3 first 
%             moments and each column correspond to different time sample
%             (at each t_c/I).
% P_500:      a vector of the probability to select an IMBH from the
%             cluster. Each cell corresponds to different time sample (at
%             each t_c/I).
% F1:         a figure handle for the BHMF evolution (if plot_coagulation
%             is true).
% F2:         a figure handle for the number of BHs vs. time(if 
%             plot_coagulation is true).
% F3:         a figure handle for the total BHs mass vs. time (if 
%             plot_coagulation is true).
% F4:         a figure handle for the total merger rate evolution (if 
%             plot_coagulation is true).
% 
function [f_src,Gamma,Gamma_tot,moments,P_500,F1,F2,F3,F4] = Coagulation(M,N0,params,str_params)
% Extract parameters from the dictionaries
M0 = params('M0');
i_max = params('i_max');
I = params('I');
Nt = params('Nt');
plot_coagulation = params('plot_coagulation');
SYM_TYPE = str_params('SYM_TYPE');
Mmax = params('Mmax');
mu_NS = params('mu_NS');
sigma_NS = params('sigma_NS');
X = params('X');
R11 = params('R11');
beta = params('beta');
gamma = params('gamma');
t_c = params('t_c');
f_loss = params('f_loss');
M_esc = params('M_esc');
v_esc = params('v_esc');
N_NS = params('N_NS');
t_NS = params('t_NS');
% Calculate the rate kernel's dependence on the masses
[i_mat,j_mat] = meshgrid(1:i_max);
R = power(min(i_mat./j_mat,j_mat./i_mat),beta).*power((i_mat+j_mat)/2,gamma);
% Calculate the recoil velocity of 2-body interactions
q = i_mat./j_mat;
A = 1.2e+4; %Km/s
B = -0.93;
H = 7.3e+3; %Km/s
K = 6e+4; %Km/s
xi = 145*pi/180; %degrees
alpha1 = 0.7*ones(i_max); alpha2 = 0.7*ones(i_max); %alpha1(:,1:floor(Mmax/M0))=0; alpha2(1:floor(Mmax/M0),:)=0;
theta1 = 0*ones(i_max); theta2 = pi*ones(i_max); theta1 = 2*pi*rand(i_max); theta2 = 2*pi*rand(i_max);
v_m = A*(q.^2.*(1-q)./((1+q).^5).*(1+B*q./((1+q).^2)));
v_perp = H*q.^2./((1+q).^5).*(alpha2.*cos(theta2)-q.*alpha1.*cos(theta1));
v_paral = K/sqrt(2)*q.^2./((1+q).^5).*(alpha2.*sin(theta2)-q.*alpha1.*sin(theta1));
v_ratios_2_body = sqrt((v_m+v_perp*cos(xi)).^2 + (v_perp*sin(xi)).^2 + v_paral.^2)/v_esc;
% Calculate the recoil velocity of 3-body interactions
v_ratios_3_body = transpose(sqrt(i_mat.*(j_mat.^2)./(M_esc/M0*(i_mat+j_mat).^2)));
% Calculate the average recoil velocity and calculate the probability to
% have ejection
[~,E] = ellipke(4*v_ratios_2_body.*v_ratios_3_body./((v_ratios_2_body+v_ratios_3_body).^2));
v_ratios = 2*(v_ratios_2_body+v_ratios_3_body)/pi.*E; v_ratios(isnan(v_ratios))=0;
P_ej = (1-power(1-v_ratios.^4,3/2)).*(v_ratios<1) + (v_ratios>=1);
% Calculate the indices for the first term in the coagulation equation
if f_loss == 0.05
    load('Mask_indices_(f_loss-005)');
elseif f_loss == 0
    load('Mask_indices_(f_loss-000)');
end
GW_ind = 16*f_loss*(i_mat.*j_mat).^2./(i_mat+j_mat).^3;
% i_ind = cell(1,Nmax); j_ind = cell(1,Nmax);
R_vec = cell(1,i_max);
for i=1:i_max
%     Mask = (i_mat+j_mat == i+round(GW_ind));
%     [j_ind{i},i_ind{i}] = find(Mask);
    R_vec{i} = diag(R(i_ind{i},j_ind{i}))';
end
% Set the time resolution
dt = 1e+9*t_c/I; %years
% Plot the IMF
if plot_coagulation
    F1 = figure; H = area(M,N0/M0); set(gca,'Layer','top');
    H.FaceColor = [0.7,0.7,0.7]; H.EdgeColor = [1,1,1];
    hold on; plot(M,N0/M0,'k');
    set(gca,'XScale','log','YScale','log');
    set(gca,'XTick',[1,2,3,5,10,30,50,100,200,1000:1000:M(end)]);
    set(gca,'FontName','Times New Roman','FontSize',20);
    Mmin = M(find(N0/M0>1,1)); Mmax = M(end);
    max_N = max(N0/M0);
    xlim([0.6*Mmin,1.2*Mmax]); ylim([1e-9*max_N,2*max_N]);
    xlabel('$M\,[M_{\odot}]$','interpreter','Latex','FontSize',24);
    ylabel('$f_{\mathrm{BH}}\,[M_\odot^{-1}]$','interpreter','Latex','FontSize',24);
    title('Black Hole Mass Function, $t$ = 0','interpreter','Latex','FontSize',18);
    hold on;
end
% Set the initial conditions
N = N0; Gamma = zeros(i_max,i_max,Nt);
Ntot = zeros(1,I+1); Ntot(1) = sum(N);
Mtot = zeros(1,I+1); Mtot(1) = M0*sum(N.*(1:i_max));
Msqr_tot = zeros(1,I+1); Msqr_tot(1) = M0^2*sum(N.*(1:i_max).^2);
Theta_mask = (i_mat+j_mat-round(GW_ind) <= i_max); R_mask = R.*Theta_mask;
Gamma(:,:,1) = R11*(N'*N).*R_mask(1:i_max,1:i_max);
P_500 = zeros(1,I+1); P_500(1) = sum(N(round(500/M0):end))/sum(N);
% At each iteration, the BHMF evolves based on previous values
for t = 1:I
%   Account for NSs contribution
    if t_NS ~= 0 && t == round(t_NS*I/t_c)
        f0_NS = exp(-(M(1:i_max)-mu_NS).^2/(2*sigma_NS^2));
        f0_NS = f0_NS/sum(f0_NS);
        N = N + N_NS*f0_NS;
    end
    N_new = N;
%   Evolve the BHMF by solving the coagulation equation numerically
    for i = 1:i_max
        if i~=1
            N_new(i) = N_new(i)+ dt*0.5*R11*sum(R_vec{i}.*N(i_ind{i}).*N(j_ind{i}));
        end
        if i < i_max
            j_max = j_ind{i_max}(find(i_ind{i_max}==i,1,'last'));
            N_new(i) = N_new(i)- dt*N(i)*R11*sum(R(i,1:j_max).*(1+P_ej(i,1:j_max)).*N(1:j_max));
        end
        if N_new(i)<0
            error('Calculation yields negative values. Increase the number of iterations.');
        end
    end
    N = N_new; 
%   Calculate moments
    Ntot(t+1) = sum(N);
    Mtot(t+1) = M0*sum(N.*(1:i_max));
    Msqr_tot(t+1) = M0^2*sum(N.*(1:i_max).^2);
%   Update Nmax
    if ~strcmp(SYM_TYPE,'Test')
        i_max = round(Mtot(t+1)/M0);
        N = N(1:i_max);
    end
%   Calculate P_500
    P_500(t+1) = sum(N(round(500/M0):end))/sum(N);
%   Calculate the merger rate Gamma_i,j(t)
    if ismember(t,round((1:(Nt-1))*I/(Nt-1)))
        Theta_mask = (i_mat+j_mat-round(GW_ind) <= i_max); R_mask = R.*Theta_mask;
        Gamma(1:i_max,1:i_max,round((Nt-1)*t/I)+1) = R11*(N'*N).*R_mask(1:i_max,1:i_max);
    end
%   Plot the BHMF at intermediate times
    if ismember(t,round((1:6)*I/7)) && plot_coagulation
        M_samples = logspace(log10(Mmin),log10(Mmax),150);
        N_samples = exp(interp1(M(1:i_max),log(N/M0),M_samples,'pchip',-Inf));
        plot(M_samples,N_samples,'.','Color',[0.094,0.5,0.764],'MarkerSize',3);
        max_N = max([max_N,max(N/M0)]);
        title(['Black Hole Mass Function, $t$ = ',num2str(t*t_c/I),' Gyears'],'interpreter','Latex','FontSize',18);
        hold on;
    end
%   Plot the final BHMF
    if t == I  && plot_coagulation
        plot(M(1:i_max),N/M0,'Color',[0.094,0.5,0.764],'LineWidth',2);
        Mmin = M(find(N/M0>1,1)); 
        if ~isempty(Mmin)
            xlim([0.6*Mmin,1.2*Mmax]);
        end
        max_N = max([max_N,max(N/M0)]);
        ylim([1e-9*max_N,2*max_N]);
        title(['Black Hole Mass Function, $t$ = ',num2str(t_c),' Gyears'],'interpreter','Latex','FontSize',18);
    end
end
% Prepare outputs
f_src = (X*[N,zeros(1,length(N0)-length(N))]+(1-X)*N0)/M0;
moments = [Ntot;Mtot;Msqr_tot];
Gamma_tot = reshape(0.5*sum(sum(Gamma)),1,Nt);
% Plot BHMF's moments and cluster's total merger rate
if plot_coagulation
%   N_tot figure
    F2 = figure; plot(t_c*(0:1/I:1),Ntot,'LineWidth',2);
    xlim([0,t_c]);
    set(gca,'FontName','Times New Roman','FontSize',20);
    xlabel('Time [Gyears]','interpreter','Latex','FontSize',24);
    ylabel('$N_\mathrm{tot}$','interpreter','Latex','FontSize',24);
    title('Number of Black Holes Vs. Time','interpreter','Latex','FontSize',20);
%   M_tot figure
    F3 = figure; plot(t_c*(0:1/I:1),Mtot,'LineWidth',2);
    xlim([0,t_c]);
    set(gca,'FontName','Times New Roman','FontSize',20);
    xlabel('Time [Gyears]','interpreter','Latex','FontSize',24);
    ylabel('$M_\mathrm{tot}\,[M_\odot]$','interpreter','Latex','FontSize',24);
    title('Total Black Holes Mass Vs. Time','interpreter','Latex','FontSize',20);
%   Total merger rate figure
    F4 = figure; plot(t_c*(0:(1/(Nt-1)):1),Gamma_tot,'LineWidth',2);
    xlim([0,t_c]); ylim([0.8*min(Gamma_tot),1.2*max(Gamma_tot)]);
    set(gca,'FontName','Times New Roman','FontSize',20);
    xlabel('Time [Gyears]','interpreter','Latex','FontSize',24);
    ylabel('$\Gamma_{\mathrm{tot}}\,[\mathrm{year}^{-1}]$','interpreter','Latex','FontSize',24);
    title('Cluster''s total merger rate Vs. Time','interpreter','Latex','FontSize',20);
else
    F1 = libpointer('doublePtr',[]);
    F2 = libpointer('doublePtr',[]);
    F3 = libpointer('doublePtr',[]);
    F4 = libpointer('doublePtr',[]);
end