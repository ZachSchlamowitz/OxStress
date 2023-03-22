%% simulate_selvaggio_single.m
% Author: Zach Schlamowitz, 2/14/2023

% This script runs the selvaggio model for one set of parameters 

%% Model

% Setup
global Params
Params = struct();
Params.v_sup = 10^(-5+6); % (µM/sec) (converted from M/sec) ranges 10^-7 to 10^-3 in Fig 3, based on varying levels of extracellular H202
Params.k_Alt = 79; % (sec^-1) default: 79
Params.k_Ox = 40; % (µM^-1 sec^-1) (converted from 4e7 M^-1 sec^-1)
Params.k_Srx = 3.3e-3; % (sec^-1) (converted from 3.3 10^-3 sec^-1)

Params.k_Sulf = 5.1e-3; % (µM^-1 sec^-1) (converted from 5.1e3 M^-1 sec^-1)
Params.k_Sulf_I = 1.5e-3; % (µM^-1 sec^-1) (converted from 5.1e3 M^-1 sec^-1)
Params.k_Sulf_II = 3.4e-3; % (µM^-1 sec^-1) (converted from 5.1e3 M^-1 sec^-1)

Params.k_Cond = 6.4; % (sec^-1) default: 6.4
Params.k_Cond_I = 11; % (sec^-1) default: 11 [see Armindo's email with cite]
Params.k_Cond_II = 0.5; % (sec^-1) default: 0.5 [see Armindo's email with cite]

Params.k_Red = 0.21; % (µM^-1 sec^-1) (converted from 2.1e5 M^-1 sec^-1)
Params.VAppMax = 230; % (µM/sec) (converted from 0.23 mM/sec) NOTE that is reported as VMax in Table 2; presumably interchangeable?
Params.K_M = 1.8; % (µM)
Params.PrxTotal = 92; % total concentration?? FLAG of Prx (µM)
Params.TrxTotal = 20; % total concentration?? FLAG of Trx (µM) (source: Armindo's email)
Params.PrxITotal = 110; % total concentration of PrxI (µM) updated: was from Supplement Table 6, now from Armindo's email
Params.PrxIITotal = 30; % total concentration of PrxII (µM) updated: was from Supplement Table 6, now from Armindo's email
Params.n_cells = 1000; %?
Params.kappa = 1; %?
Params.S = 1; %?
Params.V_medium = 3e4; %? µL
Params.V_cytoplasm = 1e-8; %? µL guess based on average MCF7 cell vol of 10000 m^3 from google search

% Run model
t0 = 0;
tf = 5000; % seconds?
% initvals = [1; 23; 23; 23; 11.5]; % For rough initial guess, start with 1µM H202 and divide the total Prx and Trx values equally across states; [PrxX] = 92/4 = 23. [TrxX] = 23/2 = 11.5
initvals = [10; 0.01; 0.001; 0.01;  0.01; 0.001; 0.01; 0.01];
% vars = [(1)H2O2, (2)PrxISO, (3)PrxISO2, (4)PrxISS, (5)PrxIISO, (6)PrxIISO2, (7)PrxIISS, (8)TrxSS]

%opts = odeset('RelTol',1e-3, 'AbsTol',1e-6, 'InitialStep',0.1, 'MaxStep',0.1);

tic
[time, sol] = ode23s(@selvaggio_model_2spec, [t0 tf], initvals);%, opts);
toc

% Obtain Trx-SH and Prx-S time courses; add as columns to sol
temp = sol'; 
temp = [temp; zeros(2, size(sol,1))];
% Note that we back out PrxS from total PrxTotal = PrxS + PrxSO + PrxSS + PrxSO2
temp(10,:) = Params.PrxTotal - temp(2,:) - temp(3,:) - temp(4,:);
% and similarly we back out TrxSH from total TrxTotal = TrxSH + TrxSS.
temp(11,:) = Params.TrxTotal - temp(5,:);
sol = temp';

%% Plotting

% Time Series Plot
plot(time, sol, LineWidth=2)
title('Species Concentrations Over Time')
xlabel('t')
ylabel('Species')
legend('H2O2', 'PrxI-SO', 'PrxI-SO2', 'PrxI-SS', 'PrxII-SO', 'PrxII-SO2', 'PrxII-SS', 'TrxSS', '-', 'PrxS','TrxSH', 'Location', 'North', 'FontSize',12 )

figure

beginning = ceil(.75*length(time));
plot(time(1:beginning), sol(1:beginning, :), LineWidth=2)
title('Beginning of Species Concentrations Over Time')
xlabel('t')
% ylabel('Population')
legend('H2O2', 'PrxI-SO', 'PrxI-SO2', 'PrxI-SS', 'PrxII-SO', 'PrxII-SO2', 'PrxII-SS', 'TrxSS', '-', 'PrxS','TrxSH', 'Location', 'North', 'FontSize',12 )
