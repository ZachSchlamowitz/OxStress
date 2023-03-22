%% plot_dimeric.m
% Script to recreate Supplementary Figure S9 from Selvaggio supplement. 
% Plots dimeric PRXII concentrations after 5min of simulation for 2 species
% model against different H2O2 boluses


% Initialize 
boluses = [0 2.5 5 10 25 50 100 250 1000]'; % micromolar values
dimeric_levels = NaN(9,1);
for i = 1:size(boluses,1)
    bolus = boluses(i);
    intracellular_val = 0.01*bolus; % intracellular is ~1% of bolus concentration
    
    ss = simulate_selvaggio(intracellular_val, NaN, 'HEK293');
    dimeric_levels(i,:) = ss;
end


%% Plot
bar(categorical(boluses), dimeric_levels)
ylim([0 1])
