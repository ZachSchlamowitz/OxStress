%% plot_dimeric.m
% Script to recreate Supplementary Figure S9 from Selvaggio supplement. 
% Plots dimeric PRXII concentrations after 5min of simulation for 2 species
% model against different H2O2 boluses


% Initialize 
boluses = [0 2.5 5 10 25 50 100 250 1000]'; % micromolar values
outputs = NaN(9,3,size(boluses,2));
dimeric_levels = NaN(9,1);
for i = 1:size(boluses,1)
    bolus = boluses(i);
    intracellular_val = 0.01*bolus; % intracellular is ~1% of bolus concentration
    
    ss = simulate_selvaggio(bolus, intracellular_val, NaN, 'HEK293');
    outputs(:,:,i) = ss;
    dimeric_levels(i,:) = ss(8,2); % 8=PRXIISS, 2=@5min
end


%% Plots
% bar(categorical(boluses), dimeric_levels)
% ylim([0 1])

% Plot bar chart vs boluses at each time point for dimeric PRX 2 and PRX 1
for j=1:3
    figure
    % get plotting data
    data1 = [outputs(5,j,:)];
    d = NaN(size(data1,3));
    for i = 1:size(data1,3)
        d(i) = data1(1,1,i);
    end
    bar(categorical(boluses), d)
    titlestr = strcat('PRXI-SS at timepoint ', num2str(j));
    title(titlestr)
    ylim([0 1])

    figure
    % get plotting data
    data2 = outputs(8,j,:);
    d = NaN(size(data2,3));
    for i = 1:size(data2,3)
        d(i) = data2(1,1,i);
    end

    bar(categorical(boluses), d)
    titlestr = strcat('PRXII-SS at timepoint ', num2str(j));
    title(titlestr)
    ylim([0 1])

    figure
    % get plotting data
    data3 = outputs(4,j,:);
    d = NaN(size(data3,3));
    for i = 1:size(data3,3)
        d(i) = data3(1,1,i);
    end

    bar(categorical(boluses), d)
    titlestr = strcat('PRXI-SO2 at timepoint ', num2str(j));
    title(titlestr)
    ylim([0 1])

    figure
    % get plotting data
    data4 = outputs(7,j,:);
    d = NaN(size(data4,3));
    for i = 1:size(data4,3)
        d(i) = data4(1,1,i);
    end

    bar(categorical(boluses), d)
    titlestr = strcat('PRXII-SO2 at timepoint ', num2str(j));
    title(titlestr)
    ylim([0 1])
    
end