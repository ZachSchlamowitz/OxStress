%% simulate_selvaggio_full_single.m
% Author: Zach Schlamowitz
%
% This script simulates the PTRS under a single experimental condition:
% either following a single bolus addition of H2O2 (using the two-prx 
% model variant), or under a specified H2O2 supply rate (using the one-prx
% model variant). It is the sister script to the function < simulate_selvaggio.m > 
% designed for interactive user use to run a single isolated simulation.
%
% To use...
% SIMULATE_SELVAGGIO_FULL_SINGLE
% 1. Set User Options
% 2. Set up desired plots in the Plot section
% 3. Run

%% User Options
num_prx = 2; % =2 or =1, specifies which version of model to use, dictating which associated user settings have effect
if num_prx == 1
    v_sup = 10^(-5); % default value: 10^(-5+6) (µM/sec) (converted from M/sec) ranges 10^-7 M/sec to 10^-3 M/sec in Fig 3, based on varying levels of extracellular H202
elseif num_prx == 2
    bolus = 1; % µM (we deem a reasonable range to be roughly >0 to 2000µM)
    intra_h2o2 = NaN; % (µM) Specify intracellular concentration of H2O2 here if desired (optional). This value will be used to set the initial value for the H2O2 state.
    cell_type = 'HEK293'; % 'HEK293' or 'MCF7  '
else
    error("Invalid nuber of peroxiredoxin species.")
end

%% Non-User Setup
global Params % global variable is used because cannot pass param struct into ode solver directly
Params = struct();

% Main if-elif-else to navigate model variant (either 1-species or 2-species)
if num_prx == 1
    %% Set Model Parameters (MCF7 used where cell lines vary)
    if isnan(v_sup) % note that this is just to catch a case where v_sup is not specified above
        Params.v_sup = 10^(-5+6); % (µM/sec) (converted from M/sec) ranges 10^-7 to 10^-3 in Fig 3, based on varying levels of extracellular H202
    else
        Params.v_sup = v_sup * 10^6; % v_sup is in M/sec so convert to µM/sec
    end

    Params.k_Alt = 79; % (sec^-1) default: 79; rate constant for H2O2 loss to alternate (exogenous) sinks
    Params.k_Ox = 40; % (µM^-1 sec^-1) default: 40 (converted from 4e7 M^-1 sec^-1) rate constant for Prx first oxidation
    Params.k_Sulf = 5.1e-3; % (µM^-1 sec^-1) default: 5.1e-3 (converted from 5.1e3 M^-1 sec^-1) rate constant for Prx hyperoxidation
    Params.k_Srx = 3.3e-3; % (sec^-1) default: 3.3e-3 (converted from 3.3 10^-3 sec^-1) rate constant for Srx activity (reduction of Prx-SO2/3)
    Params.k_Cond = 6.4; % (sec^-1) default: 6.4; rate constant for Prx condensation
    Params.k_Red = 0.21; % (µM^-1 sec^-1) default: 0.21 (converted from 2.1e5 M^-1 sec^-1) rate constant for reduction of Prx-SS to Prx-SH
    Params.VAppMax = 230; % (µM/sec) default: 230 (converted from 0.23 mM/sec) max rate of TrxR activity; NOTE that is reported as VMax in Table 2; presumably interchangeable? 
    Params.K_M = 1.8; % (µM) default: 1.8; K_M for TrxR activity
    Params.PrxTotal = 92; % (µM) default: 92; total concentration of Prx 
    Params.TrxTotal = 23; % (µM) default: 23; total concentration of Trx
    
    %% Run One-Species Model
    % ODE Solver Parameters
    % Time
    t0 = 0;
    tf = 5000; % Note: since model parameters use a time unit of seconds, this timescale is in seconds also
    
    % Initial values
    % Recall the variables are [(1)H2O2, (2)PrxSO, (3)PrxSO2, (4)PrxSS, (5)TrxSS]
    % Determine initial value for intracellular H2O2 concentration
    if isnan(intra_h2o2)
        init_h2o2 = 0.01; % default: 0.01 uM [drawn from Fig. 4 of https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5256672/]
    else % use argument passed in to simulate_selvaggio function call
        init_h2o2 = intra_h2o2;
    end
    initvals = [init_h2o2; 0.01; 0.001; 0.01; 0.01]; % init values for PRDXs and TRXs are per Armindo Salvador's suggestions
    
    % Solver options (to mess with solver tolerances and step size)
    %opts = odeset('RelTol',1e-2, 'AbsTol',1e-5, 'InitialStep',0.1, 'MaxStep',0.1);
    
    % Solve ODEs (numerically)
    tic
    [time, sol] = ode23s(@selvaggio_model, [t0 tf], initvals);%, opts); % specify solver options by including opts
    toc
    
    % Obtain Trx-SH and Prx-S time courses; add as columns to sol
    temp = sol'; 
    temp = [temp; zeros(2, size(sol,1))];
    % Note that we back out PrxS from total PrxTotal = PrxS + PrxSO + PrxSS + PrxSO2
    temp(6,:) = Params.PrxTotal - temp(2,:) - temp(3,:) - temp(4,:);
    % and similarly we back out TrxSH from total TrxTotal = TrxSH + TrxSS.
    temp(7,:) = Params.TrxTotal - temp(5,:);
    sol = temp';
    
    % Store raw trajectories
    trajectories_raw = sol;
    
    % Also store states  as (e.g., PrxII-SO) proportion of total compartment (e.g., PrxTotal)
    % Recall vars are now [(1)H2O2, (2)PrxSO, (3)PrxSO2, (4)PrxSS, (5)TrxSS, (6)PrxS, (7)TrxS]
    trajectories_frac = NaN(size(sol));
    trajectories_frac(:,1) = sol(:,1) ./ init_h2o2;  % H2O2
    trajectories_frac(:,2:4) = sol(:,2:4) ./ Params.PrxTotal; % Prx states
    trajectories_frac(:,6) = sol(:,6) ./ Params.PrxTotal;
    trajectories_frac(:,5) = sol(:,5) ./ Params.TrxTotal; % Trx States
    trajectories_frac(:,7) = sol(:,7) ./ Params.TrxTotal;


elseif num_prx == 2
    %% Set Model Parameters, depending on cell type used
    if cell_type == 'HEK293'
        % From Table 2
        Params.k_Alt = 160; % (sec^-1) default: 1.6e2; rate constant for H2O2 loss to alternate (exogenous) sinks
        Params.k_Ox = 40; % (µM^-1 sec^-1) default: 40 (converted from 4e7 M^-1 sec^-1) rate constant for Prx first oxidation
        Params.k_Srx = 4.1e-4; % (sec^-1) default: 4.1e-4 (converted from 0.41 10^-3 sec^-1) rate constant for Srx activity (reduction of PrxI/II-SO2/3)
        Params.k_Red = 0.21; % (µM^-1 sec^-1) default: 0.21 (converted from 2.1e5 M^-1 sec^-1) rate constant for reduction of PrxI/II-SS to PrxI/II-SH
        Params.VAppMax = 190; % (µM/sec) default: 190 (converted from 0.19 mM/sec) Max rate of TrxR activity; Note that is reported as VMax in Table 2; presumably interchangeable?
        Params.K_M = 1.8; % (µM) default: 1.8; K_M for TrxR activity
        Params.k_Sulf = 3.7e-3; % (µM^-1 sec^-1) default: 3.7e-3 (converted from 3.7e3 M^-1 sec^-1) rate constant for Prx hyperoxidation
        Params.k_Cond = 7.3; % (sec^-1) default: 7.3; rate constant for Prx condensation
    
        % From Supplement section 3.2.3 (text)
        Params.k_Sulf_I = 1.3e-3; % (µM^-1 sec^-1) (converted from 1.3e3 M^-1 sec^-1) rate constant for PrxI hyperoxidation
        Params.k_Sulf_II = 1.2e-2; % (µM^-1 sec^-1) (converted from 1.2e4 M^-1 sec^-1) rate constant for PrxII hyperoxidation
        
        Params.k_Cond_I = 9; % (sec^-1) default: 9 [default from Armindo's email Table with cite] rate constant for PrxI condensation
        Params.k_Cond_II = 1.7; % (sec^-1) default: 1.7 [default from Armindo's email Table with cite] rate constant for PrxII condensation
        
        % From Supplementary Table 6  
        Params.TrxTotal = 46; % (µM) default: 46; total concentration of Trx (same value in ST6)
        Params.PrxITotal = 110; % (µM) default: 110; total concentration of PrxI 
        Params.PrxIITotal = 32; % (µM) default: 32; total concentration of PrxII
    
        % Permeation Parameters (various sources)
        Params.n_cells = 3e5; % number of cells in media / well of plate; default: 3e5 cells per well from Sobotta et al 2013
        Params.kappa = 15; % permeability coefficient (µm/sec) in erythrocytes (source: 3/2023 meeting with Armindo)
        Params.V_medium = 2e11; % volume of media (µm^3) default: 2e11; converted to µm^3 from 200 µL (source: guess from Andrew)
        Params.V_cytoplasm = 1150; % volume of cell (µm^3); default: 660 from Supplement 3.2.3, although this may come from Jurkat T cells, which are smaller; alternatively, guess 1150 based on average HEK293 cell diameter of 13 µm (source: https://bionumbers.hms.harvard.edu/bionumber.aspx?id=108893&ver=3&trm=HEK+293&org= among others)
        Params.S = 530; % surface area of cell (µm^2) guess based on average HEK293 cell diameter of 13 µm (see above)
    
    elseif cell_type == 'MCF7  '
        % From Table 2
        Params.k_Alt = 79; % (sec^-1) default: 79; rate constant for H2O2 loss to alternate (exogenous) sinks
        Params.k_Ox = 40; % (µM^-1 sec^-1) default: 40 (converted from 4e7 M^-1 sec^-1) rate constant for Prx first oxidation
        Params.k_Srx = 3.3e-3; % (sec^-1) default: 3.3e-3 (converted from 3.3 10^-3 sec^-1) rate constant for Srx activity (reduction of PrxI/II-SO2/3)
        Params.k_Red = 0.21; % (µM^-1 sec^-1) default: 0.21 (converted from 2.1e5 M^-1 sec^-1) rate constant for reduction of PrxI/II-SS to PrxI/II-SH
        Params.VAppMax = 230; % (µM/sec) default: 230 (converted from 0.23 mM/sec) Max rate of TrxR activity; NOTE that is reported as VMax in Table 2; presumably interchangeable?
        Params.K_M = 1.8; % (µM) default: 1.8; K_M for TrxR activity
        
        % From Table in Armindo's Email: (formerly from Supplementary Table 6 or Supplement section 3.2.3 (text))
        Params.TrxTotal = 20; % total concentration of Trx (µM) default: 20; updated: was 23 from Supplement Table 6, now from Armindo's email
        Params.PrxITotal = 110; % total concentration of PrxI (µM) default: 110; updated: was 59 from Supplement Table 6, now from Armindo's email
        Params.PrxIITotal = 30; % total concentration of PrxII (µM) default: 30; updated: was 33 from Supplement Table 6, now from Armindo's email
        
        Params.k_Sulf_I = 1.5e-3; % (µM^-1 sec^-1) default: 1.5e-3 [source: Table from Armindo's email] rate constant for PrxI hyperoxidation
        Params.k_Sulf_II = 3.4e-3; % (µM^-1 sec^-1) default: 3.4e-3 [source: Table from Armindo's email] rate constant for PrxII hyperoxidation
        
        Params.k_Cond_I = 11; % (sec^-1) default: 11 [see Armindo's email with cite] rate constant for PrxI condensation
        Params.k_Cond_II = 0.5; % (sec^-1) default: 0.5 [see Armindo's email with cite] rate constant for PrxII condensation
    
        % Permeation Parameters (various sources)
        Params.n_cells = 6000; % default: 6000 source: guess by Andrew
        Params.kappa = 15; % permeability coefficient (µm/sec) in erythrocytes (source: 3/2023 meeting with Armindo)
        Params.V_medium = 2e11; % volume of media (µm^3) default: 2e11; converted to µm^3 from 200 µL (source: guess from Andrew)
        Params.V_cytoplasm = 1760; % volume of cell (µm^3) default: 1760 based on MCF7 cell vol of 1.76pL from https://bionumbers.hms.harvard.edu/bionumber.aspx?id=115154&ver=1
        Params.S = 1224.18; % surface area of cell (µm^2) guessed value based on SA of sphere of diamter 19.74 µm, https://www.researchgate.net/figure/Size-distribution-of-MCF-7-cells-used-in-this-study-The-histogram-was-derived-using-an_fig5_257966698

    end
    
    %% Run Two-Species Model
    % ODE Solver Parameters
    % Time (can just give start and end of time range and let solver choose
    % specific timepoints for eval or can specify timepoints by passing in a longer vector)
    %t0 = 0;
    %tf = 5000; % Note: since model parameters use a time unit of seconds, this timescale is in seconds also
    tspan = [0:1:5000]; 
    
    % Determine initial value for intracellular H2O2 concentration (unless specified by the user)
    if isnan(intra_h2o2)
        init_h2o2 = 0.01; % default: 0.01 uM [drawn from Fig. 4 of https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5256672/]
    else % use argument passed in to simulate_selvaggio function call
        init_h2o2 = intra_h2o2;
    end

    % FOR TWO-PRX MODEL W/O PERMEATION
    % vars = [(1)H2O2, (2)PrxISO, (3)PrxISO2, (4)PrxISS, (5)PrxIISO, (6)PrxIISO2, (7)PrxIISS, (8)TrxSS]
    %initvals = [init_h2o2; 0.01; 0.001; 0.01;  0.01; 0.001; 0.01; 0.01];
    
    % FOR TWO-PRX MODEL WITH PERMEATION
    initvals = [bolus; init_h2o2; 0.01; 0.001; 0.01;  0.01; 0.001; 0.01; 0.01];
    % vars = [(1)H2O2_out, (2)H2O2, (3)PrxISO, (4)PrxISO2, (5)PrxISS, (6)PrxIISO, (7)PrxIISO2, (8)PrxIISS, (9)TrxSS]
    
    % ---------------------------------------------------------------------
    % TO RUN A KNOCKOUT SIMULATION, UNCOMMENT THE APPROPRIATE SECTION HERE
    % Initial values for knockout / knockdown experiments
    % ----------------------
    % PrxI Knockout >>>
    % initvals = [bolus; init_h2o2;  0; 0; 0; 0.01; 0.001; 0.01; 0.01];
    % Params.PrxITotal = 0;
    % Params.PrxIITotal = 142; % OPTIONAL, correct for removal of PrxI pool
    % ----------------------
    % PrxII Knockout >>>
    % initvals = [bolus; init_h2o2; 0.01; 0.001; 0.01;  0; 0; 0; 0.01];
    % Params.PrxIITotal = 0;
    % Params.PrxITotal = 32; % OPTIONAL, reduce PrxI pool to PrxII pool size
    % -----------------------
    % ---------------------------------------------------------------------

    % Solver options (to mess with solver tolerances and step size)
    %opts = odeset('RelTol',1e-2, 'AbsTol',1e-5, 'InitialStep',0.1, 'MaxStep',0.1);
    
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
    trajectories_raw = sol;
    
    % Also store states (e.g., PrxII-SO) as proportion of total compartment (e.g., PrxIITotal)
    trajectories_frac = NaN(size(sol));
    trajectories_frac(:,1) = sol(:,1) ./ bolus;  % external H2O2
    trajectories_frac(:,2) = sol(:,2) ./ init_h2o2;  % internal H2O2
    trajectories_frac(:,3:5) = sol(:,3:5) ./ Params.PrxITotal; % PrxI states
    trajectories_frac(:,10) = sol(:,10) ./ Params.PrxITotal;
    trajectories_frac(:,6:8) = sol(:,6:8)./ Params.PrxIITotal; % PrxII states
    trajectories_frac(:,11) = sol(:,11) ./ Params.PrxIITotal;
    trajectories_frac(:,9) = sol(:,9) ./ Params.TrxTotal; % Trx States
    trajectories_frac(:,12) = sol(:,12) ./ Params.TrxTotal;
    
    
else

    error("ERROR: Model only supports one or two peroxiredoxin species. Please use num_prx =1 or =2 to specify model variant (see README for more).")

end

%% Plot

if num_prx == 1
    % Get trajectories of prx species
    prx_trajectories = [trajectories_frac(:,2:4)'; trajectories_frac(:,6)'];
    
    % Plots of potential interest: trajectories of each Prx state, by species
    figure
    plot(time, prx_trajectories)
    title("Prx Trajectories")
    xlabel("Time (sec)")
    ylabel("Proportion of Prx in Each Oxidation State")
    legend('PrxSO','PrxSO2', 'PrxSS', 'PrxS')
    
    
    % ADD WHATEVER PLOTS YOU LIKE HERE %
    % vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv %
    
    
    
    % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %

elseif num_prx == 2 
    % Get trajectoriees for each prx
    prxII_trajectories = [trajectories_frac(:,6:8)'; trajectories_frac(:,11)'];
    prxI_trajectories = [trajectories_frac(:,3:5)'; trajectories_frac(:,10)'];
    
    % Plots of potential interest: trajectories of each Prx state, by species
    figure
    plot(time, prxI_trajectories)
    title("PrxI Trajectories")
    xlabel("Time (sec)")
    ylabel("Proportion of Prx in Each Oxidation State")
    legend('PrxISO','PrxISO2', 'PrxISS', 'PrxIS')
    
    figure
    plot(time, prxII_trajectories)
    title("PrxII Trajectories")
    xlabel("Time (sec)")
    ylabel("Proportion of Prx in Each Oxidation State")
    legend('PrxIISO','PrxIISO2', 'PrxIISS', 'PrxIIS')
    
    % ADD WHATEVER PLOTS YOU LIKE HERE %
    % vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv %
    
    
    
    % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
end

