%% sensitivity_analysis.m
% Author: Zach Schlamowitz, 4/19/2023; updated 6/23/23
%
% This script performs sensitivity analysis on the parameters of the full 
% (two species) Selvaggio model of the Prx system. The intuition of the 
% following simple algorithm is to twiddle parameters by order of 
% magnitudes (multiply or divide by 2 or 10) and see how the resulting 
% "steady state" of the system changes. We define the steady state to be 
% just the resting state of the system after 30-40 minutes post-bolus. 
% (Note that it may be of value to consider how changing the parameters 
% affects the first 10 min or so of the simulation, even if the long-run 
% "steady states" don't change.) Then, we compute the change in the 
% position of the steady state vector under deformations of the parameters 
% by computing the Euclidean distance between the ss vector before and 
% after the deformations.
%
% The output is presented in a matrix:
% Parameter | Distance(default ss, ss with param/2) | Distance(default ss, ss with param*2) ETC.  
% ----------------------------------------------------------------------------------------- ETC.
%   k_Alt   |                 10.554                |                    8.554              ETC.
%   k_Ox    |                   0.45                |                    0.20               ETC.
%                                                  ETC. 
%
% The actual output has four columns (as opposed to the two shown above)
% corresponding to /10, /2, *2, *10.
%
% To use SENSITIVITY_ANALYSIS:
% 1. Specify user parameters as desired
% 2. Hit Run
% 3. Analyze sensitivity analysis matrix output "SAM"

%% User Options 
cell_type = 'HEK293'; % 'HEK293' or 'MCF7  '
h2o2_bolus = 1; % µM (we deem a reasonable range to be roughly >0 to 2000µM)

%% Default Parameters
default_params = struct();
default_params.v_sup = 10^(-5+6); % (µM/sec) (converted from M/sec) ranges 10^-7 to 10^-3 in Fig 3, based on varying levels of extracellular H202

if cell_type == 'HEK293'
    % From Table 2
    default_params.k_Alt = 160; % (sec^-1) default: 1.6e2; rate constant for H2O2 loss to alternate (exogenous) sinks
    default_params.k_Ox = 40; % (µM^-1 sec^-1) default: 40 (converted from 4e7 M^-1 sec^-1) rate constant for Prx first oxidation
    default_params.k_Srx = 4.1e-4; % (sec^-1) default: 4.1e-4 (converted from 0.41 10^-3 sec^-1) rate constant for Srx activity (reduction of PrxI/II-SO2/3)
    default_params.k_Red = 0.21; % (µM^-1 sec^-1) default: 0.21 (converted from 2.1e5 M^-1 sec^-1) rate constant for reduction of PrxI/II-SS to PrxI/II-SH
    default_params.VAppMax = 190; % (µM/sec) default: 190 (converted from 0.19 mM/sec) Max rate of TrxR activity; Note that is reported as VMax in Table 2; presumably interchangeable?
    default_params.K_M = 1.8; % (µM) default: 1.8; K_M for TrxR activity

    % From Supplement section 3.2.3 (text)
    default_params.k_Sulf_I = 1.3e-3; % (µM^-1 sec^-1) (converted from 1.3e3 M^-1 sec^-1) rate constant for PrxI hyperoxidation
    default_params.k_Sulf_II = 1.2e-2; % (µM^-1 sec^-1) (converted from 1.2e4 M^-1 sec^-1) rate constant for PrxII hyperoxidation
    
    default_params.k_Cond_I = 9; % (sec^-1) default: 9 [default from Armindo's email Table with cite] rate constant for PrxI condensation
    default_params.k_Cond_II = 1.7; % (sec^-1) default: 1.7 [default from Armindo's email Table with cite] rate constant for PrxII condensation
    
    % From Supplementary Table 6  
    default_params.TrxTotal = 46; % (µM) default: 46; total concentration of Trx (same value in ST6)
    default_params.PrxITotal = 110; % (µM) default: 110; total concentration of PrxI 
    default_params.PrxIITotal = 32; % (µM) default: 32; total concentration of PrxII

    % Permeation Parameters (various sources)
    default_params.n_cells = 3e5; % number of cells in media / well of plate; default: 3e5 cells per well from Sobotta et al 2013
    default_params.kappa = 15; % permeability coefficient (µm/sec) in erythrocytes (source: 3/2023 meeting with Armindo)
    default_params.V_medium = 2e11; % volume of media (µm^3) default: 2e11; converted to µm^3 from 200 µL (source: guess from Andrew)
    default_params.V_cytoplasm = 1150; % volume of cell (µm^3); default: 660 from Supplement 3.2.3, although this may come from Jurkat T cells, which are smaller; alternatively, guess 1150 based on average HEK293 cell diameter of 13 µm (source: https://bionumbers.hms.harvard.edu/bionumber.aspx?id=108893&ver=3&trm=HEK+293&org= among others)
    default_params.S = 530; % surface area of cell (µm^2) guess based on average HEK293 cell diameter of 13 µm (see above)

elseif cell_type == 'MCF7  '
    % From Table 2
    default_params.k_Alt = 79; % (sec^-1) default: 79; rate constant for H2O2 loss to alternate (exogenous) sinks
    default_params.k_Ox = 40; % (µM^-1 sec^-1) default: 40 (converted from 4e7 M^-1 sec^-1) rate constant for Prx first oxidation
    default_params.k_Srx = 3.3e-3; % (sec^-1) default: 3.3e-3 (converted from 3.3 10^-3 sec^-1) rate constant for Srx activity (reduction of PrxI/II-SO2/3)
    default_params.k_Red = 0.21; % (µM^-1 sec^-1) default: 0.21 (converted from 2.1e5 M^-1 sec^-1) rate constant for reduction of PrxI/II-SS to PrxI/II-SH
    default_params.VAppMax = 230; % (µM/sec) default: 230 (converted from 0.23 mM/sec) Max rate of TrxR activity; NOTE that is reported as VMax in Table 2; presumably interchangeable?
    default_params.K_M = 1.8; % (µM) default: 1.8; K_M for TrxR activity
    
    % From Table in Armindo's Email: (formerly from Supplementary Table 6 or Supplement section 3.2.3 (text))
    default_params.TrxTotal = 20; % total concentration of Trx (µM) default: 20; updated: was 23 from Supplement Table 6, now from Armindo's email
    default_params.PrxITotal = 110; % total concentration of PrxI (µM) default: 110; updated: was 59 from Supplement Table 6, now from Armindo's email
    default_params.PrxIITotal = 30; % total concentration of PrxII (µM) default: 30; updated: was 33 from Supplement Table 6, now from Armindo's email
    
    default_params.k_Sulf_I = 1.5e-3; % (µM^-1 sec^-1) default: 1.5e-3 [source: Table from Armindo's email] rate constant for PrxI hyperoxidation
    default_params.k_Sulf_II = 3.4e-3; % (µM^-1 sec^-1) default: 3.4e-3 [source: Table from Armindo's email] rate constant for PrxII hyperoxidation
    
    default_params.k_Cond_I = 11; % (sec^-1) default: 11 [see Armindo's email with cite] rate constant for PrxI condensation
    default_params.k_Cond_II = 0.5; % (sec^-1) default: 0.5 [see Armindo's email with cite] rate constant for PrxII condensation

    % Permeation Parameters (various sources)
    default_params.n_cells = 6000; % default: 6000 source: guess by Andrew
    default_params.kappa = 15; % permeability coefficient (µm/sec) in erythrocytes (source: 3/2023 meeting with Armindo)
    default_params.S = 1224.18; % surface area of cell (µm^2) guessed value based on SA of sphere of diamter 19.74 µm, https://www.researchgate.net/figure/Size-distribution-of-MCF-7-cells-used-in-this-study-The-histogram-was-derived-using-an_fig5_257966698
    default_params.V_medium = 2e11; %? default: 2e11 µm^3; converted to µm^3 from 200 µL (source: guess from Andrew)
    default_params.V_cytoplasm = 1760; %? µm^3 based on MCF7 cell vol of 1.76pL from https://bionumbers.hms.harvard.edu/bionumber.aspx?id=115154&ver=1
end

% Get the default state vector (under a given H2O2 bolus)
[trajectories_raw, trajectories_frac] = simulate_selvaggio(h2o2_bolus, NaN, default_params.v_sup, 'HEK293', 2); %bolus=h2o2_bolus, intra_h2o2=NaN, v_sup=NaN, cell_type='HEK293',num_prx=2)
vec_default = trajectories_frac(2100,:); % take state vector at 35min = 2100sec

% Now make new parameter struct that will be used in the sensitivity analysis sims
Params = default_params;
global Params

% Get container for all the parameters' names
fn = fieldnames(Params);

% Initialize sensitivity analysis output matrix
SAM = NaN(numel(fn),4); % 4 columns correspond to /10, /2, *2, *10 distances
for k=1:numel(fn) % loop over struct fields by fieldname
    %% Twiddle parameter down (/10)
    Params.(fn{k}) = default_params.(fn{k}) / 10;
    
    % Simulate with this parameter set
    %[ss_down_10, hm_v_down_10] = simulate_selvaggio_full_single(1,NaN,NaN,'HEK293',Params);
    [ss_down_10_raw, ss_down_10_frac] = run_2spec_sim(h2o2_bolus, NaN);

    % Store output "steady state" vector
    vec_down_10 = ss_down_10_frac(2100,:); % take state vector at 35min = 2100sec

    %% Twiddle parameter down (/2)
    Params.(fn{k}) = default_params.(fn{k}) / 2;
    
    % Simulate with this parameter set
    [ss_down_2_raw, ss_down_2_frac] = run_2spec_sim(h2o2_bolus, NaN);
    
    % Store output "steady state" vector
    vec_down_2 = ss_down_2_frac(2100,:); % take state vector at 35min = 2100sec

    %% Twiddle parameter up (*2)
    Params.(fn{k}) = default_params.(fn{k}) * 2;
    
    % Simulate with this parameter set
    [ss_up_2_raw, ss_up_2_frac] = run_2spec_sim(h2o2_bolus, NaN);
    
    % Store output "steady state" vector
    vec_up_2 = ss_up_2_frac(2100,:); % take state vector at 35min = 2100sec
    
    %% Twiddle parameter up (*10)
    Params.(fn{k}) = default_params.(fn{k}) * 10;
    
    % Simulate with this parameter set
    [ss_up_10_raw, ss_up_10_frac] = run_2spec_sim(h2o2_bolus, NaN);
    
    % Store output "steady state" vector
    vec_up_10 = ss_up_10_frac(2100,:); % take state vector at 35min = 2100sec
    
    %% Compute distances to default vector 
    SAM(k,1) = norm(vec_down_10 - vec_default);
    SAM(k,2) = norm(vec_down_2  - vec_default);
    SAM(k,3) = norm(vec_up_10   - vec_default);
    SAM(k,4) = norm(vec_up_2    - vec_default);

    %% Reset this parameter
    Params.(fn{k}) = default_params.(fn{k});

end

toc


%% Script Functions

% Function to run a single simulation of 2-species PTRS
function [traj_raw, traj_frac] = run_2spec_sim(bolus, intra_h2o2)
    % This function is stolen from simulate_selvaggio but housed here for clarity
    global Params

    %% Run Two-Species Model
    % ODE Solver Parameters
    % Time (can just give start and end of time range and let solver choose
    % specific timepoints for eval or can specify timepoints by passing in a longer vector; this is what we opt for here)
    %t0 = 0;
    %tf = 5000;
    tspan = [0:1:5000]; % Note: since model parameters use a time unit of seconds, this timescale is in seconds also
    
    % Determine initial value for intracellular H2O2 concentration
    if isnan(intra_h2o2)
        init_h2o2 = 0.01; % default: 0.01 uM [drawn from Fig. 4 of https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5256672/]
    else % use argument passed in to simulate_selvaggio function call
        init_h2o2 = intra_h2o2;
    end

    % FOR TWO-PRX MODEL W/O PERMEATION
    % vars are [(1)H2O2, (2)PrxISO, (3)PrxISO2, (4)PrxISS, (5)PrxIISO, (6)PrxIISO2, (7)PrxIISS, (8)TrxSS]
    % initvals = [init_h2o2; 0.01; 0.001; 0.01;  0.01; 0.001; 0.01; 0.01];
    
    % FOR TWO-PRX MODEL WITH PERMEATION
    initvals = [bolus; init_h2o2; 0.01; 0.001; 0.01;  0.01; 0.001; 0.01; 0.01];
    % vars are [(1)H2O2_out, (2)H2O2, (3)PrxISO, (4)PrxISO2, (5)PrxISS, (6)PrxIISO, (7)PrxIISO2, (8)PrxIISS, (9)TrxSS]
    
    % Solve ODEs (numerically)
    tic
    [time, sol] = ode23s(@selvaggio_model_2spec_perm, tspan, initvals);%, opts); % specify solver options by including opts
    toc
    
    % Obtain Trx-SH and PrxI/II-S time courses; add as columns to sol
    temp = sol'; 
    temp = [temp; zeros(3, size(sol,1))];
    % Note that we back out PrxS from total PrxTotal = PrxS + PrxSO + PrxSS +
    % PrxSO2 for both PRX1 and PRX2
    temp(10,:) = Params.PrxITotal - temp(3,:) - temp(4,:) - temp(5,:); % PRX1-S
    temp(11,:) = Params.PrxIITotal - temp(6,:) - temp(7,:) - temp(8,:); % PRX2-S
    % and similarly we back out TrxSH from total TrxTotal = TrxSH + TrxSS.
    temp(12,:) = Params.TrxTotal - temp(9,:);
    sol = temp';

    % Store raw trajectories
    traj_raw = sol;
    
    % Also store states (e.g., PrxII-SO) as proportion of total compartment (e.g., PrxIITotal)
    traj_frac = NaN(size(sol));
    traj_frac(:,1) = sol(:,1) ./ bolus;  % external H2O2
    traj_frac(:,2) = sol(:,2) ./ init_h2o2;  % internal H2O2
    traj_frac(:,3:5) = sol(:,3:5) ./ Params.PrxITotal; % PrxI states
    traj_frac(:,10) = sol(:,10) ./ Params.PrxITotal;
    traj_frac(:,6:8) = sol(:,6:8)./ Params.PrxIITotal; % PrxII states
    traj_frac(:,11) = sol(:,11) ./ Params.PrxIITotal;
    traj_frac(:,9) = sol(:,9) ./ Params.TrxTotal; % Trx States
    traj_frac(:,12) = sol(:,12) ./ Params.TrxTotal;
end
