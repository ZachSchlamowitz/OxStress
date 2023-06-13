%% plot_ptrs.m
% Author: Zach Schlamowitz
% Last Updated: 6/12/23
%
% This file is the main script for Zach Schlamowitz's 2023 Senior Honors 
% Thesis project at the University of Arizona. The script runs multi-condition 
% simulations of the peroxiredoxin-thioredoxin system (PTRS) as described by the 
% mathematical model put forth by Selvaggio et al. (2018). The model uses 
% ordinary differential equations (ODEs) to specify two species of
% peroxiredoxin (prx) acting inside a single cell. Please see the
% README.txt file for this project (accessible via:
% https://github.com/ZachSchlamowitz/OxStress). Clerical note: this script 
% was formerly named plot_dimeric.m
%
% The core experiment we simulate is bolus addition of hydrogen peroxide
% into the extracellular matrix, reminiscent of an in vitro experiment. The
% resulting PTRS protein concentrations in a single cell are simulated over
% time, producing discrete timecourses as simulation outputs. Various doses
% for the bolus addition are simulated. Among user-specified settings are
% options to toggle between the 1-prx and 2-prx models, change cell types
% (for the 2-prx model), or change which bolus doses are simulated.
%
% To use this script, simply set desired parameters in the User Options
% section below, and hit run! The script will generate two data structs 
% containing the aforementoined timecourses, showing, for each dose, the 
% concentrations of each PTRS species every second for 5000 seconds. These 
% gross values are contained in one struct (gross_traj), while proportional
% values giving each species' concentration as a fraction of the total 
% compartment (e.g., PrxII-SS / PrxIITotal) are contained in another (frac_traj).
% See the plotting section below to examine and, if desired, modify figures 
% that are generated from this core simulation. 

%% User Options
% vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

% Initialize simulation settings
num_prx = 1; % =2 or =1, specifies which version of model to use
% if num_prx is =1, ignore all other settings in this section and instead
% look at intracellular_val in 1-Prx Logic section
boluses = [0.001 2.5 5 10 25 50 100 250 500 1000 2000]'; % extracellular H2O2 boluses to simulate (values are micromolar)
bolus_list = {'dose1', 'dose2', 'dose3', 'dose4', 'dose5', 'dose6', 'dose7', 'dose8', 'dose9', 'dose10', 'dose11'};%, 'dose12', 'dose13'}; % Modify this to match length of boluses
cell_type = 'HEK293'; % ='HEK293' or 'MCF7  '; specifies which cell line's parameters to use

% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

%% Internal Initialization
gross_traj = struct;
frac_traj = struct;

%% 1-Prx Logic
if num_prx == 1
    
    % Specify intracellular concentration of H2O2 here if desired (optional).
    % This value will be used to set the initial value for the H2O2 state.
    intracellular_val = NaN; % default: NaN

    % Use plot_phenotypes function to loop over values of v_sup and call
    % appropriate model scripts for the 1-prx model. 
    [PrxS_fig, PrxSO_fig, PrxSO2_fig] = plot_phenotypes(intracellular_val, num_prx);

%% 2-Prx Logic
elseif num_prx == 2
    v_sup = NaN; % v_sup not used in 2-prx model
    
    % Loop over bolus doses and simulate PTRS under each
    for i = 1:size(boluses,1)
        bolus = boluses(i);
        intracellular_val = min(0.01*bolus, 0.01); % intracellular is ~1% of bolus concentration, but...
        % for small bolus we assume a resting level of h2o2 in the cell of 0.01 uM;
        % for large bolus this is negligible
    
        % Run simulation
        [trajectories_raw, trajectories_frac] = simulate_selvaggio(bolus, intracellular_val, v_sup, cell_type, num_prx);
    
        % Store simulation results
        gross_traj.(bolus_list{i}) = trajectories_raw;
        frac_traj.(bolus_list{i}) = trajectories_frac;

    end
    
    %% Pull desired data for plotting
    
    % Initialize matrices for plotting data (we re-arrange simulation data
    % into "plotting data" for clarity). The name of each matrix specifies
    % the timecourse it contains.
    frac_prxII_hyperox_matrix = NaN(size(trajectories_raw,1),size(boluses,1)); % PrxII-SO2 (frac)
    frac_prxI_hyperox_matrix = NaN(size(frac_prxII_hyperox_matrix)); % PrxI-SO2 (frac)
    frac_prxII_disulfide_matrix = NaN(size(frac_prxII_hyperox_matrix)); % PrxII-SS (frac)
    frac_prxI_disulfide_matrix = NaN(size(frac_prxII_hyperox_matrix)); % PrxI-SS (frac)
    prxII_disulfide_matrix = NaN(size(frac_prxII_hyperox_matrix)); % PrxII-SS (gross)
    prxI_disulfide_matrix = NaN(size(frac_prxII_hyperox_matrix)); % PrxI-SS (gross)
    
    % Sort data into the plotting matrices
    for j = 1:size(boluses,1)
        % Get desired trajectories under each dose j
        frac_prxII_hyperox_matrix(:,j) = frac_traj.(bolus_list{j})(:,7); % PrxII-SO2 (frac)
        frac_prxI_hyperox_matrix(:,j) = frac_traj.(bolus_list{j})(:,4); % PrxI-SO2 (frac)
        frac_prxII_disulfide_matrix(:,j) = frac_traj.(bolus_list{j})(:,8); % PrxII-SS (frac)
        frac_prxI_disulfide_matrix(:,j) = frac_traj.(bolus_list{j})(:,5); % PrxI-SS (frac)
        prxII_disulfide_matrix(:,j) = gross_traj.(bolus_list{j})(:,8); % PrxII-SS (gross)
        prxI_disulfide_matrix(:,j) = gross_traj.(bolus_list{j})(:,5); % PrxI-SS (gross)
    end
    
    % If desired, uncomment this to: 
    % % Rescale prct and raw disulfide Prx to max seen in that timecourse
    % [prxII_prct_disulfide_maxes, prxII_prct_max_idxs] = max(frac_prxII_disulfide_matrix(:,:));
    % [prxI_prct_disulfide_maxes, prxI_prct_max_idxs] = max(frac_prxI_disulfide_matrix(:,:));
    % [prxII_disulfide_maxes, prxII_max_idxs] = max(prxII_disulfide_matrix(:,:));
    % [prxI_disulfide_maxes, prxI_max_idxs] = max(prxI_disulfide_matrix(:,:));
    % 
    % prxII_prct_disulfide_matrix_scaled = frac_prxII_disulfide_matrix./prxII_prct_disulfide_maxes;
    % prxI_prct_disulfide_matrix_scaled = frac_prxI_disulfide_matrix./prxI_prct_disulfide_maxes;
    % prxII_disulfide_matrix_scaled = prxII_disulfide_matrix./prxII_disulfide_maxes;
    % prxI_disulfide_matrix_scaled = prxI_disulfide_matrix./prxI_disulfide_maxes;
    
    
    %% Plots
    % Plot 1 ------------------------------------
    % Replicating Selvaggio Supplemental Figure 9
    % -------------------------------------------
    prct_prxII_disulfide_vs_dose = frac_prxII_disulfide_matrix(300,:)'; % get prxII disulfide proportion @ 5 min post-bolus
    figure
    bar(categorical(boluses), prct_prxII_disulfide_vs_dose)
    ylim([0 0.1]) % NOTE this is the zoomed version, shown in the thesis document, not the scale shown in Selvaggio/Sobotta.
    title('Proportion PrxII-SS vs Dose at 5min Post-Bolus')
    xlabel('H2O2 Bolus (µM)')
    ylabel('Proportion Dimeric PrxII')
    % -------------------------------------------
    
    
    % Plot 2 ----------------------------------
    % Hyperoxidaiton Heatmaps 
    % -----------------------------------------
    
    % To plot the heatmaps so that time is increasing along y axis, we flip the
    % relevant data vertically:
    frac_prxII_hyperox_matrix_plotting = frac_prxII_hyperox_matrix(end:-1:1,:);
    frac_prxI_hyperox_matrix_plotting = frac_prxI_hyperox_matrix(end:-1:1,:);
    
    % PrxII Hyperoxidation Heatmap
    figure
    h2 = heatmap(frac_prxII_hyperox_matrix_plotting(2501:5000,:)); % We plot the first half of the time sereies only to focus on first ~40min
    h2.GridVisible = 'off'; % gets rid of gridlines (because at this scale they obscure the colors)
    h2.XLabel = 'Dose (μM)';
    h2.YLabel = 'Time (sec)';
    h2.Title = 'Fraction PrxII-SO2 vs Time and Bolus Dose';
    h2.caxis([0 1]); % coloraxis scale
    h2.XData = {'0', '2.5', '5', '10', '25', '50', '100', '250', '500', '1000', '2000'};
    % Clean up y ticks:
    Y_labels = [2500:-1:1]';
    % Convert each number in the array into a string
    CustomYLabels = string(Y_labels);
    % Replace all but the fifth elements by spaces (so they won't show)
    CustomYLabels(mod(Y_labels,600) ~= 0) = " ";
    h2.YDisplayLabels = CustomYLabels; % Clean up y ticks
    
    % PrxI Hyperoxidation Heatmap
    figure
    h1 = heatmap(frac_prxI_hyperox_matrix_plotting(2501:5000,:));
    h1.GridVisible = 'off';
    h1.XLabel = 'Dose (μM)';
    h1.YLabel = 'Time (sec)';
    h1.Title = 'Fraction PrxI-SO2 vs Time and Bolus Dose';
    h1.caxis([0 1]);
    h1.XData = {'0', '2.5', '5', '10', '25', '50', '100', '250', '500', '1000', '2000'};
    h1.YDisplayLabels = CustomYLabels; % Clean up y ticks
    
    % ---------------------------------
    
    
    % Plot 3 ------------------------------
    %  Hyperoxidation Dose Response Curve
    %     ("Andrew's Favorite Figure")
    % -------------------------------------
    
    % Obtain percent hyperoxidized across different doses at a post-30min 
    % timepoint, compile into vector for plotting
    prxII_hyperox_vs_dose = frac_prxII_hyperox_matrix(2100,:)'; % 2100sec = 35min
    prxI_hyperox_vs_dose = frac_prxI_hyperox_matrix(2100,:)';
  
    % % (optional) Do the same for percent disulfide: 
    % prxII_disulfide_vs_dose = frac_prxII_disulfide_matrix(4700,:)';
    % prxI_disulfide_vs_dose = frac_prxI_disulfide_matrix(4700,:)';
    
    % % Plot figure using linear axis scaling:
    % figure
    % p1 = plot(boluses, prxII_hyperox_vs_dose, '-o', MarkerSize=10, LineWidth=2);
    % p1.Color(4) = 0.75;
    % title("Fraction Hyperoxidized Prxs vs Dose at 35min Post-Bolus")
    % xlabel('Dose (μM)')
    % ylabel('% Hyperoxidzed Prx')
    % hold on
    % p2 = plot(boluses, prxI_hyperox_vs_dose, '-s', MarkerSize=10,LineWidth=2);
    % p2.Color(4) = 0.5;
    % legend('PrxII-SO2', 'PrxI-SO2')
    
    % Plot dose response of each Prx on log-scaled dose (x) axis
    figure
    % Plot prxII
    s1 = semilogx(boluses, prxII_hyperox_vs_dose, '-o', MarkerSize=10, LineWidth=2.5);%, Color="#77AC30");%  MarkerFaceColor="#0072BD") % color tags were used when buliding knockout plots which require running this code block multiple times, in which case default colors can cause confusion
    hold on
    % Plot prxI
    s2 = semilogx(boluses, prxI_hyperox_vs_dose, '-s', MarkerSize=10,LineWidth=2.5);%, Color="#EDB120");% , MarkerFaceColor="#D95319")
    % For clarity, we increase the transparency by reducing the opacity of the curves
    s1.Color(4) = 0.75;
    s2.Color(4) = 0.5; % reduce opacity for the darker color more
    xlabel('Dose (μM)')
    ylim([-0.01 0.75])
    xlim([10^0 10^4])
    ylabel('Fraction Hyperoxidzed Prx')
    legend('PrxII-SO2', 'PrxI-SO2')
    title("Fraction Hyperoxidized Prxs vs Dose at 35min Post-Bolus")
    
    % % If wish to add disulfide dose-response curves as well, uncomment:
    % hold on
    % figure
    % s1 = semilogx(boluses, prxII_disulfide_vs_dose, '--o', MarkerSize=10, LineWidth=2, Color="#0072BD");%  MarkerFaceColor="#0072BD")
    % hold on
    % s2 = semilogx(boluses, prxI_disulfide_vs_dose, '--s', MarkerSize=10,LineWidth=2, Color="#D95319");% , MarkerFaceColor="#D95319")
    % s1.Color(4) = 0.75;
    % s2.Color(4) = 0.5;
    % xlabel('Dose (μM)')
    % ylim([-0.01 0.75])
    % xlim([10^0 10^4])
    % ylabel('Prx Fraction')
    % % legend('PrxII-SO2', 'PrxI-SO2', 'PrxII-SS', 'PrxI-SS')
    % legend('PrxII-SS', 'PrxI-SS')
    % title("Percent Disulfide Prx vs Dose at 10 min Post-Bolus")
    %title("Percent Hyperoxidized and Disulfide Prx vs Dose at 35min Post-Bolus")
    % ------------------------------------
    
    
    % % Plot 4 --------------------------
    % %     (%) PrxI/II-SS vs Dose
    % % ---------------------------------
    %
    % NOTE: the following are variations on dose-response curves for the
    % disulfide prx states. They did not provide a clear takeaway nor use 
    % upon first examination, and were left untouched thereafter. Hence,
    % there are two versions. Feel free to play with them as you see fit. 
    %
    %
    % % Plot attempt 1
    % % Obtain percent disulfide scaled to max percent disulfide across different 
    % % doses at post-30min timepoint, compile into vector for plotting
    % prxII_prct_disulfide_vs_dose_fixedtpt = prxII_prct_disulfide_matrix_scaled(4400,:)';
    % prxI_prct_disulfide_vs_dose_fixedtpt = prxI_prct_disulfide_matrix_scaled(4400,:)';
    % 
    % prxII_disulfide_vs_dose_fixedtpt = prxII_disulfide_matrix_scaled(4940,:)';
    % prxI_disulfide_vs_dose_fixedtpt = prxI_disulfide_matrix_scaled(4940,:)';
    % 
    % figure
    % s1 = semilogx(boluses, prxII_disulfide_vs_dose_fixedtpt, '-o', MarkerSize=10, LineWidth=2);%, Color="#0072BD");%  MarkerFaceColor="#0072BD")
    % hold on
    % s2 = semilogx(boluses, prxI_disulfide_vs_dose_fixedtpt, '-s', MarkerSize=10,LineWidth=2);%, Color="#EDB120");% , MarkerFaceColor="#D95319")
    % s1.Color(4) = 0.75;
    % s2.Color(4) = 0.5;
    % xlabel('Dose (μM)')
    % ylim([-0.01 1.01])
    % xlim([10^0 10^4])
    % ylabel('% Disulfide Prx (scaled to max)')
    % legend('PrxII-SO2', 'PrxI-SO2')
    % title("Percent Disulfide Prx's vs Dose at 10min Post-Bolus")
    % 
    % 
    % % Plot attempt 2
    % % Just use prxII_disulfide_maxes, that is, plot the max values achieved
    % % figure
    % % s1 = semilogx(boluses, prxII_prct_disulfide_maxes, '-o', MarkerSize=10, LineWidth=2);%, Color="#0072BD");%  MarkerFaceColor="#0072BD")
    % % hold on
    % % s2 = semilogx(boluses, prxI_prct_disulfide_maxes, '-s', MarkerSize=10,LineWidth=2);%, Color="#EDB120");% , MarkerFaceColor="#D95319")
    % % s1.Color(4) = 0.75;
    % % s2.Color(4) = 0.5;
    % % xlabel('Dose (μM)')
    % % % ylim([-0.01 0.75])
    % % xlim([10^0 10^4])
    % % ylabel('Max % Disulfide Prx')
    % % legend('PrxII-SO2', 'PrxI-SO2')
    % % title("Maximum Percent Disulfide Prx's vs Dose")
    % -------------------------------------

    

    % Plot 5 ------------------------------
    %        Disulfide Timecourses
    % -------------------------------------
    % Prx-SS vs Time (for each dose)
    for j = 1:size(boluses,1)
        figure
        plot(frac_prxII_disulfide_matrix(:,j), LineWidth=2)
        hold on
        plot(frac_prxI_disulfide_matrix(:,j), LineWidth=2)
        xlim([0 600])
        xlabel('Time (sec)')
        ylabel('Fraction Prx-SS')
        legend('PrxII-SS', 'PrxI-SS')
        titlestr = strcat('Prx Disulfide Activity Time Course Following ', num2str(boluses(j)), 'uM Bolus');
        title(titlestr)
    end
    % -------------------------------------

else
    error("num_prx must be either =1 or =2. No other number of prxs supported.")
end


