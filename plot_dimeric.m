%% plot_dimeric.m
% Script to recreate Supplementary Figure S9 from Selvaggio supplement. 
% Plots dimeric PRXII concentrations after 5min of simulation for 2 species
% model against different H2O2 boluses


% Initialize 
boluses = [0.001 2.5 5 10 25 50 100 250 500 1000 2000]'; % micromolar values
outputs = NaN(9,3,size(boluses,2)); % this is just 1!
dimeric_levels = NaN(size(boluses,1),1); % this is for Plot 1 -- replicating % dimeric PrxII at 5min post-bolus
heatmap_data = struct;
bolus_list = {'dose1', 'dose2', 'dose3', 'dose4', 'dose5', 'dose6', 'dose7', 'dose8', 'dose9', 'dose10', 'dose11'};
cell_type = 'HEK293'; % 'HEK293' or 'MCF7  '
for i = 1:size(boluses,1)
    bolus = boluses(i);
    intracellular_val = min(0.01*bolus, 0.01); % intracellular is ~1% of bolus concentration, but...
    % for small bolus we assume a resting level of h2o2 in the cell of 0.01 uM;
    % for large bolus this is negligible

    [ss, heatmap_trajectories] = simulate_selvaggio(bolus, intracellular_val, NaN, cell_type);
    
    outputs(:,:,i) = ss;
    dimeric_levels(i,:) = ss(8,2); % 8=PRXIISS, 2=@5min

    heatmap_data.(bolus_list{i}) = heatmap_trajectories;

end


%% Plots
% Plot 1 ------------------------------------
% Replicating Selvaggio Supplemental Figure 9
% -------------------------------------------
%
% bar(categorical(boluses), dimeric_levels)
% ylim([0 1])
% -------------------------------------------



% Plot 2 ----------------------------------
% PRXI/II-S02/SS vs boluses at 3 timepoints
% -----------------------------------------
% Plot bar chart vs boluses at each time point for dimeric PRX 2 and PRX 1
% for j=1:3
%     figure
%     % get plotting data
%     data1 = [outputs(5,j,:)];
%     d = NaN(size(data1,3));
%     for i = 1:size(data1,3)
%         d(i) = data1(1,1,i);
%     end
%     bar(categorical(boluses), d)
%     titlestr = strcat('PRXI-SS at timepoint ', num2str(j));
%     title(titlestr)
%     ylim([0 1])
% 
%     figure
%     % get plotting data
%     data2 = outputs(8,j,:);
%     d = NaN(size(data2,3));
%     for i = 1:size(data2,3)
%         d(i) = data2(1,1,i);
%     end
% 
%     bar(categorical(boluses), d)
%     titlestr = strcat('PRXII-SS at timepoint ', num2str(j));
%     title(titlestr)
%     ylim([0 1])
% 
%     figure
%     % get plotting data
%     data3 = outputs(4,j,:);
%     d = NaN(size(data3,3));
%     for i = 1:size(data3,3)
%         d(i) = data3(1,1,i);
%     end
% 
%     bar(categorical(boluses), d)
%     titlestr = strcat('PRXI-SO2 at timepoint ', num2str(j));
%     title(titlestr)
%     ylim([0 1])
% 
%     figure
%     % get plotting data
%     data4 = outputs(7,j,:);
%     d = NaN(size(data4,3));
%     for i = 1:size(data4,3)
%         d(i) = data4(1,1,i);
%     end
% 
%     bar(categorical(boluses), d)
%     titlestr = strcat('PRXII-SO2 at timepoint ', num2str(j));
%     title(titlestr)
%     ylim([0 1])
% 
% end
% ---------------------------------


% Plot 3 ----------------------------------
% Heatmaps of "Our Stat" and Hyperoxidaiton
% -----------------------------------------
prxII_over_prxsum_matrix = NaN(size(heatmap_trajectories,1),size(boluses,2)); % matrix housing "our stat"
prxII_hyperox_matrix = NaN(size(prxII_over_prxsum_matrix));
prxI_hyperox_matrix = NaN(size(prxII_over_prxsum_matrix));
prxI_over_prxsum_matrix = NaN(size(prxII_over_prxsum_matrix));

for j = 1:size(boluses,1)
    % column 1 of heatmap_data is time
    prxII_over_prxsum_matrix(:,j) = heatmap_data.(bolus_list{j})(end:-1:1,2);
    prxII_hyperox_matrix(:,j) = heatmap_data.(bolus_list{j})(end:-1:1,3);
    prxI_hyperox_matrix(:,j) = heatmap_data.(bolus_list{j})(end:-1:1,4);
    prxI_over_prxsum_matrix(:,j) = heatmap_data.(bolus_list{j})(end:-1:1,5);
end

% PrxII/PrxSum Heatmap
% figure
% h = heatmap(prxII_over_prxsum_matrix(2501:5000,:));
% h.GridVisible = 'off';
% h.XLabel = 'Dose (μM)';
% h.YLabel = 'Time (sec)';
% h.Title = '(PrxIISS)/(PrxISS + PrxIISS) Heatmap';
% h.caxis([0 1]);
% h.XData = {'0', '2.5', '5', '10', '25', '50', '100', '250', '500', '1000', '2000'};
% % Clean up y ticks:
% Y_labels = [2500:-1:1]';
% % Convert each number in the array into a string
% CustomYLabels = string(Y_labels);
% % Replace all but the fifth elements by spaces
% CustomYLabels(mod(Y_labels,600) ~= 0) = " ";
% h.YDisplayLabels = CustomYLabels;
% 
% % PrxI/PrxSum Heatmap
% figure
% h3 = heatmap(prxI_over_prxsum_matrix(2501:5000,:));
% h3.GridVisible = 'off';
% h3.XLabel = 'Dose (μM)';
% h3.YLabel = 'Time (sec)';
% h3.Title = '(PrxISS)/(PrxISS + PrxIISS) Heatmap';
% h3.caxis([0 1]);
% h3.XData = {'0', '2.5', '5', '10', '25', '50', '100', '250', '500', '1000', '2000'};
% h3.YDisplayLabels = CustomYLabels; % Clean up y ticks
% 
% % PrxII Hyperoxidation Heatmap
% figure
% h2 = heatmap(prxII_hyperox_matrix(2501:5000,:));
% h2.GridVisible = 'off';
% h2.XLabel = 'Dose (μM)';
% h2.YLabel = 'Time (sec)';
% h2.Title = 'PrxII-SO2 Heatmap';
% h2.caxis([0 1]);
% h2.XData = {'0', '2.5', '5', '10', '25', '50', '100', '250', '500', '1000', '2000'};
% h2.YDisplayLabels = CustomYLabels; % Clean up y ticks
% 
% % PrxI Hyperoxidation Heatmap
% figure
% h1 = heatmap(prxI_hyperox_matrix(2501:5000,:));
% h1.GridVisible = 'off';
% h1.XLabel = 'Dose (μM)';
% h1.YLabel = 'Time (sec)';
% h1.Title = 'PrxI-SO2 Heatmap';
% h1.caxis([0 1]);
% h1.XData = {'0', '2.5', '5', '10', '25', '50', '100', '250', '500', '1000', '2000'};
% h1.YDisplayLabels = CustomYLabels; % Clean up y ticks

% ---------------------------------




% Plot 4 --------------------------
%     Andrew's Favorite Plot
% ---------------------------------

% Obtain percent hyperoxidized across different doses at post-30min 
% timepoint, compile into vector for plotting
prxII_hyperox_vs_dose = prxII_hyperox_matrix(400,:)';
prxI_hyperox_vs_dose = prxI_hyperox_matrix(400,:)';

figure
p1 = plot(boluses, prxII_hyperox_vs_dose, '-o', MarkerSize=10, LineWidth=2);
p1.Color(4) = 0.75;
title("Andrew's Favorite Figure")
xlabel('Dose (μM)')
ylabel('% Hyperoxidzed Prx')
hold on
p2 = plot(boluses, prxI_hyperox_vs_dose, '-s', MarkerSize=10,LineWidth=2);
p2.Color(4) = 0.5;
legend('PrxII-SO2', 'PrxI-SO2')

figure
s1 = semilogx(boluses, prxII_hyperox_vs_dose, '-o', MarkerSize=10, LineWidth=2);%  MarkerFaceColor="#0072BD")
hold on
s2 = semilogx(boluses, prxI_hyperox_vs_dose, '-s', MarkerSize=10,LineWidth=2);% , MarkerFaceColor="#D95319")
s1.Color(4) = 0.75;
s2.Color(4) = 0.5;
xlabel('Dose (μM)')
ylim([-0.01 0.25])
xlim([10^0 10^4])
ylabel('% Hyperoxidzed Prx')
legend('PrxII-SO2', 'PrxI-SO2')
title("Percent Hyperoxidized Prx's vs Dose at 35min Post-Bolus")
% ---------------------------------