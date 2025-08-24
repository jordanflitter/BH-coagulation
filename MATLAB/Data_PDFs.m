% This script produces PDFs for the variables M1, M2, q and z given their
% best fit values and 1-sigma uncertainties for multiple merger events that
% were detected during an experiment. For M1, M2, and z the PDFs are
% combinations of asymmetric Gaussians. To make the PDF of q, the M1, M2
% PDFs of each event are sampled multiple times to make the event's PDF,
% and then all the PDFs are combined together
close all; clc;
% Load the O3a data and set the variables vectors
load('O3a_data');
M1_data_vec = logspace(log10(0.8*min(M1_data(:,1)-M1_data(:,3))),log10(1.2*max(M1_data(:,1)+M1_data(:,2))),1000);
M2_data_vec = logspace(log10(0.8*min(M2_data(:,1)-M2_data(:,3))),log10(1.2*max(M2_data(:,1)+M2_data(:,2))),1000);
z_data_vec = logspace(log10(0.8*min(z_data(:,1)-z_data(:,3))),log10(1.2*max(z_data(:,1)+z_data(:,2))),1000);
% Calculate dM1 and dM2, the widths of the above samples in the spetra
r_M1 = M1_data_vec(2)/M1_data_vec(1);
r_M2 = M2_data_vec(2)/M2_data_vec(1);
dM1_vec = M1_data_vec*(sqrt(r_M1)-1/sqrt(r_M1));
dM2_vec = M2_data_vec*(sqrt(r_M2)-1/sqrt(r_M2));
% Set intial conditions for the PDFs
f_M1 = 0; f_M2 = 0; f_z = 0; f_q = 0;
% At each step in the loop the PDFs of a single event are calculated and
% then added to the whole PDFs
for i = 1:size(M1_data,1)
%   Calcualate the event's PDF for M1
    f_M1_plus = (M1_data_vec>=M1_data(i,1)).*exp(-(M1_data_vec-M1_data(i,1)).^2/(2*M1_data(i,2)^2)); 
    f_M1_plus = f_M1_plus/trapz(M1_data_vec,f_M1_plus);
    f_M1_minus = (M1_data_vec<M1_data(i,1)).*exp(-(M1_data_vec-M1_data(i,1)).^2/(2*M1_data(i,3)^2));
    f_M1_minus = f_M1_minus/trapz(M1_data_vec,f_M1_minus);
    f_M1_event = (f_M1_plus+f_M1_minus)/2;
    f_M1 = f_M1 + f_M1_event;
%   Calcualate the event's PDF for M2
    f_M2_plus = (M2_data_vec>=M2_data(i,1)).*exp(-(M2_data_vec-M2_data(i,1)).^2/(2*M2_data(i,2)^2)); 
    f_M2_plus = f_M2_plus/trapz(M2_data_vec,f_M2_plus);
    f_M2_minus = (M2_data_vec<M2_data(i,1)).*exp(-(M2_data_vec-M2_data(i,1)).^2/(2*M2_data(i,3)^2));
    f_M2_minus = f_M2_minus/trapz(M2_data_vec,f_M2_minus);
    f_M2_event = (f_M2_plus+f_M2_minus)/2;
    f_M2 = f_M2 + f_M2_event;
%   Calcualate the event's PDF for z
    f_z_plus = (z_data_vec>=z_data(i,1)).*exp(-(z_data_vec-z_data(i,1)).^2/(2*z_data(i,2)^2)); 
    f_z_plus = f_z_plus/trapz(z_data_vec,f_z_plus);
    f_z_minus = (z_data_vec<z_data(i,1)).*exp(-(z_data_vec-z_data(i,1)).^2/(2*z_data(i,3)^2));
    f_z_minus = f_z_minus/trapz(z_data_vec,f_z_minus);
    f_z_event = (f_z_plus+f_z_minus)/2;
    f_z = f_z + f_z_event;
%   Assign weights to the samples of the PDFs of M1 and M2
    W_M1 = dM1_vec.*f_M1_event;
    W_M2 = dM2_vec.*f_M2_event;
%   Generate samples of the mass ratio q which are taken from the M1 and M2
%   distributions
    q_values = zeros(1,10000);
    for j = 1:length(q_values)
        M1_rand = randsample(M1_data_vec,1,true,W_M1);
        M2_rand = randsample(M2_data_vec(M2_data_vec<=M1_rand),1,true,W_M2(M2_data_vec<=M1_rand));
        q_rand = M2_rand/M1_rand;
        q_values(j) = q_rand;
    end
%   Calculate the event's PDF of q by making a histogram of the randomly 
%   sampled q values 
    [f_q_event,q_edges] = histcounts(q_values,19); f_q_event = f_q_event*length(f_q_event);
    qbins = (q_edges(1:end-1)+q_edges(2:end))/2;
%   Interpolate the PDF samples in fixed values of q
    q_data_vec = linspace(0,1,length(f_q_event)+2); q_data_vec(end) = []; q_data_vec(1) = [];
    f_q_event = interp1(qbins,f_q_event,q_data_vec,'pchip',0);
    f_q_event = f_q_event/trapz(q_data_vec,f_q_event);
    f_q = f_q + f_q_event;
end