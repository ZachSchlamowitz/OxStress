%% simulate_selvaggio.m
% Author: Zach Schlamowitz, 2/14/2023

%% Model
function [steady_states] = simulate_selvaggio(v_sup)
% This function gets the steady state values of PrxS and PrxSO from a
% simulation of the selvaggio ODE model with H2O2 supply rate v_sup

% Setup
global Params
Params = struct();
Params.v_sup = v_sup * 10^6; % v_sup is in M/sec so convert to µM/sec
% Params.v_sup = 10^(-5+6); % (µM/sec) (converted from M/sec) ranges 10^-7 to 10^-3 in Fig 3, based on varying levels of extracellular H202
Params.k_Alt = 79; % (sec^-1) default: 79
Params.k_Ox = 40; % (µM^-1 sec^-1) (converted from 4e7 M^-1 sec^-1)
Params.k_Sulf = 5.1e-3; % (µM^-1 sec^-1) (converted from 5.1e3 M^-1 sec^-1)
Params.k_Srx = 3.3e-3; % (sec^-1) (converted from 3.3 10^-3 sec^-1)
Params.k_Cond = 6.4; % (sec^-1) default: 6.4
Params.k_Red = 0.21; % (µM^-1 sec^-1) (converted from 2.1e5 M^-1 sec^-1)
Params.VAppMax = 230; % (µM/sec) (converted from 0.23 mM/sec) NOTE that is reported as VMax in Table 2; presumably interchangeable?
Params.K_M = 1.8; % (µM)
Params.PrxTotal = 92; % total concentration?? FLAG of Prx (µM)
Params.TrxTotal = 23; % total concentration?? FLAG of Trx (µM)

% Run model
t0 = 0;
tf = 5000; % seconds?
% initvals = [1; 23; 23; 23; 11.5]; % For rough initial guess, start with 1µM H202 and divide the total Prx and Trx values equally across states; [PrxX] = 92/4 = 23. [TrxX] = 23/2 = 11.5
initvals = [1; 0.1; 0.1; 0.1; 0.1];
%opts = odeset('RelTol',1e-2, 'AbsTol',1e-5, 'InitialStep',0.1, 'MaxStep',0.1);

tic
[time, sol] = ode23s(@selvaggio_model, [t0 tf], initvals);%, opts);
toc

% Obtain Trx-SH and Prx-S time courses; add as columns to sol
temp = sol'; 
temp = [temp; zeros(2, size(sol,1))];
% Note that we back out PrxS from total PrxTotal = PrxS + PrxSO + PrxSS + PrxSO2
temp(6,:) = Params.PrxTotal - temp(2,:) - temp(3,:) - temp(4,:);
% and similarly we back out TrxSH from total TrxTotal = TrxSH + TrxSS.
temp(7,:) = Params.TrxTotal - temp(5,:);
sol = temp';

% Pull steady-state values of desired species
steady_states = NaN(1,2);
steady_states(1,1) = sol(end,6)/Params.PrxTotal; % PrxS
steady_states(1,2) = sol(end,2)/Params.PrxTotal; % PrxSO
steady_states(1,3) = sol(end,3)/Params.PrxTotal; % PrxSO2
end