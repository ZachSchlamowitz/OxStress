%% plot_phenotypes.m
% Script to plot various "phenotypic regions" arising from subsets of the "design space"

function [PrxS_fig, PrxSO_fig, PrxSO2_fig] = plot_phenotypes(intracellular_val, num_prx)
    %% Obtain species trajectories vs v_sup
    % Initialize 
    v_sups = zeros(61,1);
    steady_states = NaN(61,3);
    for i = 1:61
        v_sup_exp = -8+0.1*(i-1);
        v_sup_val = 10^(v_sup_exp);
        v_sups(i,1) = v_sup_val;
    
        [trajectories_raw, trajectories_frac] = simulate_selvaggio(NaN, intracellular_val, v_sup_val, NaN, num_prx);
        
        if num_prx ~= 1
            error("Expected 1-prx model (num_prx=1) but was passed different value for num_prx.")
        end
        
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
    xlabel("v_{sup} (μM/sec)", FontSize=12)
    ylabel("Fraction PrxS", FontSize=12)
    
    PrxSO_fig = figure;
    loglog(v_sups, steady_states(:,2), 'o')
    title("Steady State Fraction of PrxSO versus v_{sup}", FontSize=12)
    ylim([10^-7 1])
    xlim([10^-8 10^-2])
    xlabel("v_{sup} (μM/sec)", FontSize=12)
    ylabel("Fraction PrxSO", FontSize=12)
    
    PrxSO2_fig = figure;
    loglog(v_sups, steady_states(:,3), 'o')
    title("Steady State Fraction of PrxSO2 versus v_{sup}", FontSize=12)
    ylim([10^-7 1])
    xlim([10^-8 10^-2])
    xlabel("v_{sup} (μM/sec)", FontSize=12)
    ylabel("Fraction PrxSO2", FontSize=12)

end