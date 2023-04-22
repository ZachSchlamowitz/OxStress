%% simulate_selvaggio_full_single.m
% Author: Zach Schlamowitz, 2/14/2023

%% Model
function [steady_states, heatmap_vals] = simulate_selvaggio_full_single(bolus, intra_h2o2, v_sup, cell_type, Params)
% This function gets the steady state values of PrxS and PrxSO from a
% simulation of the selvaggio ODE model with H2O2 supply rate v_sup
global Params

% if Params == struct()
% Setup
Params = struct();
if isnan(v_sup)
    Params.v_sup = 10^(-5+6); % (µM/sec) (converted from M/sec) ranges 10^-7 to 10^-3 in Fig 3, based on varying levels of extracellular H202
else
    Params.v_sup = v_sup * 10^6; % v_sup is in M/sec so convert to µM/sec
end

if cell_type == 'HEK293'
    % From Table 2
    Params.k_Alt = 160; % (sec^-1) default: 1.6e2 (source: Table 2)
    Params.k_Ox = 40; % (µM^-1 sec^-1) (converted from 4e7 M^-1 sec^-1) 
    Params.k_Srx = 4.1e-4; % (sec^-1) (converted from 0.41 10^-3 sec^-1)
    Params.k_Red = 0.21; % (µM^-1 sec^-1) (converted from 2.1e5 M^-1 sec^-1)
    Params.VAppMax = 190; % (µM/sec) (converted from 0.19 mM/sec) NOTE that is reported as VMax in Table 2; presumably interchangeable?
    Params.K_M = 1.8; % (µM)
    Params.PrxTotal = 140; % total concentration of Prx (µM)
    Params.TrxTotal = 46; % total concentration of Trx (µM) (same value in ST6)
    Params.k_Sulf = 3.7e-3; % (µM^-1 sec^-1) (converted from 3.7e3 M^-1 sec^-1)
    Params.k_Cond = 7.3; % (sec^-1) default: 7.3

    % From Supplement section 3.2.3 (text)
    Params.k_Sulf_I = 1.3e-3; % (µM^-1 sec^-1) (converted from 1.3e3 M^-1 sec^-1)
    Params.k_Sulf_II = 1.2e-2; % (µM^-1 sec^-1) (converted from 1.2e4 M^-1 sec^-1)

    Params.k_Cond_I = 9; % (sec^-1) default: 11 [default from Armindo's email Table with cite]
    Params.k_Cond_II = 1.7; % (sec^-1) default: 0.5 [default from Armindo's email Table with cite]

    % From Supplementary Table 6  
    Params.PrxITotal = 110; % total concentration of PrxI (µM)
    Params.PrxIITotal = 32; % total concentration of PrxII (µM)

    % Permeation Parameters (various sources)
    Params.n_cells = 3e5; % number of cells in media / well of plate; default: 3e5 cells per well from Sobotta et al 2013
    Params.kappa = 15; % permeability coefficient (µm/sec) in erythrocytes (source: 3/2023 meeting with Armindo)
    Params.V_medium = 2e11; % volume of media (µm^3) default: 2e11; converted to µm^3 from 200 µL (source: guess from Andrew)
    Params.V_cytoplasm = 1150; % volume of cell (µm^3); default: 660 from Supplement 3.2.3, although this may come from Jurkat T cells, which are smaller; alternatively, guess 1150 based on average HEK293 cell diameter of 13 µm (source: https://bionumbers.hms.harvard.edu/bionumber.aspx?id=108893&ver=3&trm=HEK+293&org= among others)
    Params.S = 530; % surface area of cell (µm^2) guess based on average HEK293 cell diameter of 13 µm (see above)

elseif cell_type == 'MCF7  '
    % From Table 2
    Params.k_Alt = 79; % (sec^-1) default: 79
    Params.k_Ox = 40; % (µM^-1 sec^-1) (converted from 4e7 M^-1 sec^-1)
    Params.k_Srx = 3.3e-3; % (sec^-1) (converted from 3.3 10^-3 sec^-1)
    Params.k_Red = 0.21; % (µM^-1 sec^-1) (converted from 2.1e5 M^-1 sec^-1)
    Params.VAppMax = 230; % (µM/sec) (converted from 0.23 mM/sec) NOTE that is reported as VMax in Table 2; presumably interchangeable?
    Params.K_M = 1.8; % (µM)
    Params.k_Sulf = 5.1e-3; % (µM^-1 sec^-1) (converted from 5.1e3 M^-1 sec^-1)
    Params.k_Cond = 6.4; % (sec^-1) default: 6.4
    Params.PrxTotal = 92; % total concentration of Prx (µM) (NOTE: total = Params.PrxITotal + Params.PrxIITotal if use ST6 values)

    % From Table in Armindo's Email: (formerly from Supplementary Table 6 or Supplement section 3.2.3 (text))
    Params.TrxTotal = 20; % total concentration of Trx (µM) default: 20; updated: was 23 from Supplement Table 6, now from Armindo's email
    Params.PrxITotal = 110; % total concentration of PrxI (µM) default: 110; updated: was 59 from Supplement Table 6, now from Armindo's email
    Params.PrxIITotal = 30; % total concentration of PrxII (µM) default: 30; updated: was 33 from Supplement Table 6, now from Armindo's email

    Params.k_Sulf_I = 1.5e-3; % (µM^-1 sec^-1) default: 1.5e-3 [source: Table from Armindo's email]
    Params.k_Sulf_II = 3.4e-3; % (µM^-1 sec^-1) default: 3.4e-3 [source: Table from Armindo's email]

    Params.k_Cond_I = 11; % (sec^-1) default: 11 [see Armindo's email with cite]
    Params.k_Cond_II = 0.5; % (sec^-1) default: 0.5 [see Armindo's email with cite]

    % Permeation Parameters (various sources)
    Params.n_cells = 6000; % default: 6000 source: guess by Andrew
    Params.kappa = 15; % permeability coefficient (µm/sec) in erythrocytes (source: 3/2023 meeting with Armindo)
    Params.S = 1224.18; % surface area of cell (µm^2) guessed value based on SA of sphere of diamter 19.74 µm, https://www.researchgate.net/figure/Size-distribution-of-MCF-7-cells-used-in-this-study-The-histogram-was-derived-using-an_fig5_257966698
    Params.V_medium = 2e11; %? default: 2e11 µm^3; converted to µm^3 from 200 µL (source: guess from Andrew)
    Params.V_cytoplasm = 1760; %? µm^3 based on MCF7 cell vol of 1.76pL from https://bionumbers.hms.harvard.edu/bionumber.aspx?id=115154&ver=1

end
% end

% Run model
t0 = 0;
tf = 5000; % seconds

tspan = [0:1:5000]; % can specify specific timepoints for eval or just give start and end of time range and let solver choose

if isnan(intra_h2o2)
    init_h2o2 = 0.01; % default: 0.01 uM [drawn from Fig. 4 of https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5256672/]
else
    init_h2o2 = intra_h2o2;
end

% FOR ONE-PRX MODEL
% vars = [(1)H2O2, (2)PrxSO, (3)PrxSO2, (4)PrxSS, (5)TrxSS]
% initvals = [init_h2o2; 0.01; 0.001; 0.01; 0.01]; % init values for PRDXs and TRXs are per Armindo's suggestions

% FOR TWO-PRX MODEL W/O PERMEATION
% vars = [(1)H2O2, (2)PrxISO, (3)PrxISO2, (4)PrxISS, (5)PrxIISO, (6)PrxIISO2, (7)PrxIISS, (8)TrxSS]
% initvals = [init_h2o2; 0.01; 0.001; 0.01;  0.01; 0.001; 0.01; 0.01];

% FOR TWO-PRX MODEL WITH PERMEATION
initvals = [bolus; init_h2o2; 0.01; 0.001; 0.01;  0.01; 0.001; 0.01; 0.01];

% Initial values for knockout / knockdown experiments
initvals = [bolus; init_h2o2; 0.01; 0.001; 0.01;  0; 0; 0; 0.01];
Params.PrxITotal = 0;

% vars = [(1)H2O2_out, (2)H2O2, (3)PrxISO, (4)PrxISO2, (5)PrxISS, (6)PrxIISO, (7)PrxIISO2, (8)PrxIISS, (9)TrxSS]
%opts = odeset('RelTol',1e-2, 'AbsTol',1e-5, 'InitialStep',0.1, 'MaxStep',0.1);

tic
[time, sol] = ode23s(@selvaggio_model_2spec_perm, tspan, initvals);%, opts);
toc

% Obtain Trx-SH and Prx-S time courses; add as columns to sol
temp = sol'; 
temp = [temp; zeros(3, size(sol,1))];
% Note that we back out PrxS from total PrxTotal = PrxS + PrxSO + PrxSS +
% PrxSO2 for both PRX1 and PRX2
temp(10,:) = Params.PrxITotal - temp(2,:) - temp(3,:) - temp(4,:); % PRX1-S
temp(11,:) = Params.PrxIITotal - temp(5,:) - temp(6,:) - temp(7,:); % PRX2-S
% and similarly we back out TrxSH from total TrxTotal = TrxSH + TrxSS.
temp(12,:) = Params.TrxTotal - temp(8,:);
sol = temp';

% Initialize output matrix 
steady_states = NaN(12,size(time,1)); % *note that depending on the application this may or may not actually houes steady states

% Pull desired steady-state* values of all species as fracs of total
% steady_states(1,:) = sol(:,1)' ./ bolus;
% steady_states(2,:) = sol(:,2)' ./ init_h2o2;
% steady_states(3:5,:) = sol(:,3:5)' ./ Params.PrxITotal;
% steady_states(6:8,:) = sol(:,6:8)'./ Params.PrxIITotal;
% steady_states(9,:) = sol(:,9)' ./ Params.TrxTotal;
steady_states(:,:) = sol(:,:)';
% vars = [(1)H2O2_out, (2)H2O2, (3)PrxISO, (4)PrxISO2, (5)PrxISS, (6)PrxIISO, (7)PrxIISO2, (8)PrxIISS, (9)TrxSS] (10)PrxIS (11)PrxIIS (12)TrxSH 

% Compute percent hyperoxidation timecourses
prct_PrxIISO2 = steady_states(7,:)';
prct_PrxISO2  = steady_states(4,:)';

% Compute our heatmap statistic: Prx2-SS/(Prx2-SS + Prx1-SS)
prxII_over_prxSum = sol(:,8)./(sol(:,8)+sol(:,5));
prxI_over_prxSum = sol(:,5)./(sol(:,8)+sol(:,5));
heatmap_vals = [time, prxII_over_prxSum, prct_PrxIISO2, prct_PrxISO2, prxI_over_prxSum];

% Plot
prxI_trajectories = [sol(:,6:8)'; sol(:,11)'];
prxII_trajectories = [sol(:,3:5)'; sol(:,10)'];

% figure
% plot(time, prxI_trajectories)
% legend('PrxISO','PrxISO2', 'PrxISS', 'PrxIS')
% 
% figure
% plot(time, prxII_trajectories)
% legend('PrxIISO','PrxIISO2', 'PrxIISS', 'PrxIIS')


end