%% plot_dimeric.m
% Script to recreate Supplementary Figure S9 from Selvaggio supplement. 
% Plots dimeric PRXII concentrations after 5min of simulation for 2 species
% model against different H2O2 boluses


% Initialize 
boluses = [0 2.5 5 10 25 50 100 250 500 1000 2000]'; % micromolar values
outputs = NaN(9,3,size(boluses,2)); % this is just 1!
dimeric_levels = NaN(size(boluses,1),1);
heatmap_data = struct;
bolus_list = {'dose1', 'dose2', 'dose3', 'dose4', 'dose5', 'dose6', 'dose7', 'dose8', 'dose9', 'dose10', 'dose11'};
for i = 1:size(boluses,1)
    bolus = boluses(i);
    intracellular_val = 0.01*bolus; % intracellular is ~1% of bolus concentration

    [ss,heatmap_vals] = simulate_selvaggio(bolus, intracellular_val, NaN, 'HEK293');
    outputs(:,:,i) = ss;
    dimeric_levels(i,:) = ss(8,2); % 8=PRXIISS, 2=@5min

    heatmap_data.(bolus_list{i}) = heatmap_vals;

end


%% Plots
% Plot 1 --------------------------
% bar(categorical(boluses), dimeric_levels)
% ylim([0 1])
% ---------------------------------



% Plot 2 --------------------------
% 1. PRXI/II-S02/SS vs boluses at 3 timepoints
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


% Plot 3 --------------------------
% heatmap(heatmap_data)
dose_time_stat_matrix = NaN(size(heatmap_vals,1),size(boluses,2));
prxII_hyperox_matrix = NaN(size(dose_time_stat_matrix));
prxI_hyperox_matrix = NaN(size(dose_time_stat_matrix));
prxI_over_prxsum_matrix = NaN(size(dose_time_stat_matrix));

for j = 1:size(boluses,1)
    dose_time_stat_matrix(:,j) = heatmap_data.(bolus_list{j})(:,2);
    prxII_hyperox_matrix(:,j) = heatmap_data.(bolus_list{j})(:,3);
    prxI_hyperox_matrix(:,j) = heatmap_data.(bolus_list{j})(:,4);
    prxI_over_prxsum_matrix(:,j) = heatmap_data.(bolus_list{j})(:,5);
end

figure
h = heatmap(dose_time_stat_matrix(1:2500,:));
h.GridVisible = 'off';
h.XLabel = 'Dose (μM)';
h.YLabel = 'Time (sec)';
h.Title = '(PrxIISS)/(PrxISS + PrxIISS) Heatmap';
h.caxis([0 1]);
h.XData = {'0', '2.5', '5', '10', '25', '50', '100', '250', '500', '1000', '2000'};

figure
h3 = heatmap(prxI_over_prxsum_matrix(1:2500,:));
h3.GridVisible = 'off';
h3.XLabel = 'Dose (μM)';
h3.YLabel = 'Time (sec)';
h3.Title = '(PrxISS)/(PrxISS + PrxIISS) Heatmap';
h3.caxis([0 1]);
h3.XData = {'0', '2.5', '5', '10', '25', '50', '100', '250', '500', '1000', '2000'};

figure
h2 = heatmap(prxII_hyperox_matrix(1:2500,:));
h2.GridVisible = 'off';
h2.XLabel = 'Dose (μM)';
h2.YLabel = 'Time (sec)';
h2.Title = 'PrxII-SO2 Heatmap';
h2.caxis([0 1]);
h2.XData = {'0', '2.5', '5', '10', '25', '50', '100', '250', '500', '1000', '2000'};

figure
h1 = heatmap(prxI_hyperox_matrix(1:2500,:));
h1.GridVisible = 'off';
h1.XLabel = 'Dose (μM)';
h1.YLabel = 'Time (sec)';
h1.Title = 'PrxI-SO2 Heatmap';
h1.caxis([0 1]);
h1.XData = {'0', '2.5', '5', '10', '25', '50', '100', '250', '500', '1000', '2000'};


% ---------------------------------




% Plot 4 --------------------------
%     Andrew's Favorite Plot
% ---------------------------------

% Obtain percent hyperoxidized across different doses at post-30min 
% timepoint, compile into vector for plotting
prxII_hyperox_vs_dose = prxII_hyperox_matrix(2000,:)';
prxI_hyperox_vs_dose = prxI_hyperox_matrix(2000,:)';

figure
plot(boluses, prxII_hyperox_vs_dose, 'o')
title("Andrew's Favorite Figure")
xlabel('Dose (μM)')
ylabel('% Hyperoxidzed Prx')
hold on
plot(boluses, prxI_hyperox_vs_dose, 'o')
legend('PrxII-SO2', 'PrxI-SO2')

figure
semilogx(boluses, prxII_hyperox_vs_dose, 'o', MarkerSize=12, MarkerFaceColor="#D95319")
hold on
semilogx(boluses, prxI_hyperox_vs_dose, 'o', MarkerSize=12)
xlabel('Dose (μM)')
ylabel('% Hyperoxidzed Prx')
legend('PrxII-SO2', 'PrxI-SO2')
% ---------------------------------