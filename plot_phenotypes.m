%% plot_phenotypes.m
% Author: Zach Schlamowitz
%
% This function creates (and returns handles to) three figures which 
% replicate results from in Selvaggio et al. (2018). The figures show 
% steady state concentrations of species of the peroxiredoxin-thioredoxin 
% system (PTRS) at various values of one model parameter: v_sup, the H2O2 
% supply rate. The name "plot_phenotypes" is used in reference to the 
% various "phenotypic regions" Selvaggio et al. (2018) identify as 
% qualitative outcomes of their ODE model of the PTRS. The replicated 
% figures generated here are figures Fig.4a-c of our manuscript, 
% Schlamowitz (2023), which replicate curves seen in Figure 4i of Selvaggio
% et al. (2018).
%
% PLOT_PHENOTYPES
% Parameters: 
% --> intracellular_val: intracellular concentration of H2O2 (uM)
% --> num_prx: number of peroxiredoxin species to model (=1 or =2), sets model version
% Returns:
% <-- PrxS_fig: figure handle for steady state values of PrxS vs v_sup values
% <-- PrxSO_fig: figure handle for steady state values of PrxSO vs v_sup values
% <-- PrxSO2_fig: figure handle for steady state values of PrxSO2 vs v_sup values
% Other notes:
% - plots each of the above figures as well (toggle on/off via commenting)

function [PrxS_fig, PrxSO_fig, PrxSO2_fig] = plot_phenotypes(intracellular_val, num_prx)
    %% Obtain species trajectories vs v_sup
    % Initialize 
    v_sups = zeros(61,1);
    steady_states = NaN(61,3);

    % Loop over v_sup values logarithmically and simulate PTRS under each
    for i = 1:61
        v_sup_exp = -8+0.1*(i-1); % exponent for v_sup
        v_sup_val = 10^(v_sup_exp); % compute v_sup value (M/sec)
        v_sups(i,1) = v_sup_val; % store v_sup value
        
        if num_prx ~= 1
            error("Expected 1-prx model (num_prx=1) but was passed different value for num_prx.")
        end

        % Simulate PTRS model with the inputted intracellular initial value
        % for H2O2 and the current H2O2 supply value (v_sup)
        [trajectories_raw, trajectories_frac] = simulate_selvaggio(NaN, intracellular_val, v_sup_val, NaN, num_prx);     
        
        % Pull steady-state values of desired species 
        % (NOTE: that here we assume the simulation has run long enough to 
        % reach a steady state so are just pulling the last values of the 
        % timecourses. Check this assumption if change duration of simulation.) 
        ss = NaN(1,3);
        ss(1,1) = trajectories_frac(end,6); % PrxS
        ss(1,2) = trajectories_frac(end,2); % PrxSO
        ss(1,3) = trajectories_frac(end,3); % PrxSO2
        
        % Store steady states for this v_sup val
        steady_states(i,:) = ss;
    end
    
    
    %% Plot
    PrxS_fig = figure;
    loglog(v_sups, steady_states(:,1), 'o')
    title("Steady State Fraction of PrxS versus v_{sup}", FontSize=12)
    ylim([10^-7 1])
    xlim([10^-8 10^-2])
    xlabel("v_{sup} (M/sec)", FontSize=12)
    ylabel("Fraction PrxS", FontSize=12)
    
    PrxSO_fig = figure;
    loglog(v_sups, steady_states(:,2), 'o')
    title("Steady State Fraction of PrxSO versus v_{sup}", FontSize=12)
    ylim([10^-7 1])
    xlim([10^-8 10^-2])
    xlabel("v_{sup} (M/sec)", FontSize=12)
    ylabel("Fraction PrxSO", FontSize=12)
    
    PrxSO2_fig = figure;
    loglog(v_sups, steady_states(:,3), 'o')
    title("Steady State Fraction of PrxSO2 versus v_{sup}", FontSize=12)
    ylim([10^-7 1])
    xlim([10^-8 10^-2])
    xlabel("v_{sup} (M/sec)", FontSize=12)
    ylabel("Fraction PrxSO2", FontSize=12)

end