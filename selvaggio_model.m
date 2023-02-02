% selvaggio_model.m

% This is an implementation of the Selvaggio (2018) PRX system ODE model
%% Template
% x' = x - axy
% y' = -y + Bxy
% 
% t0 = 0;
% tf = 15;
% initvals = [20; 20];
% [time, sol] = ode23(@lotka, [t0 tf], initvals);
% 
% % Time Series Plot
% plot(time, sol)
% title('Predator, Prey Populations Over Time')
% xlabel('t')
% ylabel('Population')
% legend('Prey', 'Predators', 'Location', 'North')
% 
% % Phase Plane Plot
% plot(sol(:,1), sol(:,2))
% title('Phase Space')
% xlabel('Prey Population')
% ylabel('Predator Population')
% 
% % Add in version with ode45 as well
% [T,S] = ode45(@lotka, [t0 tf], initvals);
% plot(sol(:,1), sol(:,2), '-', S(:,1), S(:,2), '-');
% title('Comparing Phase Plane Portraits')
% legend('ode23', 'ode45')
% 
% function eqs = lotka(t,initvals)
%     a = 0.01;
%     B = 0.02;
%     eqs = zeros(2,1);
%     eqs(1) = initvals(1) - a*initvals(1)*initvals(2);
%     eqs(2) = -initvals(2) + B*initvals(1)*initvals(2);
% end

%% Equations

% Initialize parameters
global Params
Params = struct();
Params.v_sup = 10^(-5+6); % (µM/sec) (converted from M/sec) ranges 10^-7 to 10^-3 in Fig 3, based on varying levels of extracellular H202
Params.k_Alt = 79; % (sec^-1)
Params.k_Ox = 40; % (µM^-1 sec^-1) (converted from 4e7 M^-1 sec^-1)
Params.k_Sulf = 5.1e-3; % (µM^-1 sec^-1) (converted from M^-1 sec^-1)
Params.k_Srx = 3.3e-3; % (sec^-1) (converted from 10^-3 sec^-1)
Params.k_Cond = 6.4; % (sec^-1)
Params.k_Red = 0.21; % (µM^-1 sec^-1) (converted from M^-1 sec^-1)
Params.VAppMax = 230; % (µM/sec) (converted from 0.23 mM/sec) NOTE that is reported as VMax in Table 2; presumably interchangeable?
Params.K_M = 1.8; % (µM)
Params.PrxTotal = 92; % total concentration?? FLAG of Prx (µM)
Params.TrxTotal = 23; % total concentration?? FLAG of Trx (µM)

% Run Model
t0 = 0;
tf = 50; % seconds?
% initvals = [1; 23; 23; 23; 11.5]; % For rough initial guess, start with 1µM H202 and divide the total Prx and Trx values equally across states; [PrxX] = 92/4 = 23. [TrxX] = 23/2 = 11.5
initvals = [10; 0.1; 0.1; 0.1; 0.1];
%opts = odeset('RelTol',1e-2, 'AbsTol',1e-5, 'InitialStep',0.1, 'MaxStep',0.1);
tic
[time, sol] = ode23s(@selvaggio, [t0 tf], initvals);%, opts);
toc
% disp(strcat('Time elapsed: ', num2str(toc)))

% Obtain Trx-SH and Prx-S time courses; add as columns to sol
temp = sol'; 
temp = [temp; zeros(2, size(sol,1))];
% Note that we back out PrxS from total PrxTotal = PrxS + PrxSO + PrxSS + PrxSO2
temp(6,:) = Params.PrxTotal - temp(2,:) - temp(3,:) - temp(4,:);
% and similarly we back out TrxSH from total TrxTotal = TrxSH + TrxSS.
temp(7,:) = Params.TrxTotal - temp(5,:);
sol = temp';

% Time Series Plot
plot(time, sol, LineWidth=2)
% title('Predator, Prey Populations Over Time')
xlabel('t')
% ylabel('Population')
legend('H2O2', 'PrxSO', 'PrxSO2', 'PrxSS', 'TrxSS', 'PrxS','TrxSH', 'Location', 'North', 'FontSize',14 )

figure

beginning = ceil(.75*length(time));
plot(time(1:beginning), sol(1:beginning, :), LineWidth=2)
% title('Predator, Prey Populations Over Time')
xlabel('t')
% ylabel('Population')
legend('H2O2', 'PrxSO', 'PrxSO2', 'PrxSS', 'TrxSS','PrxS','TrxSH','Location', 'North', 'FontSize',14 )

function eqs = selvaggio(t,vars)
    global Params
    % Key:
    % vars = [(1)H2O2, (2)PrxSO, (3)PrxSO2, (4)PrxSS, (5)TrxSS]
    % Note that we back out PrxS from total PrxTotal = PrxS + PrxSO + PrxSS + PrxSO2
    PrxS = Params.PrxTotal - vars(2) - vars(3) - vars(4);
    % and similarly we back out TrxSH from total TrxTotal = TrxSH + TrxSS.
    TrxSH = Params.TrxTotal - vars(5);

    eqs = zeros(5,1);
    % dH202/dt
    % eqs(1) = v_sup - k_Alt*H2O2 - k_Ox*PrxS*H2O2 - k_Sulf*PrxSO*H2O2;
    eqs(1) = Params.v_sup - Params.k_Alt*vars(1) - Params.k_Ox*PrxS*vars(1) - Params.k_Sulf*vars(2)*vars(1);
    % dPrxSO/dt
    % eqs(2) = k_Ox*PrxS*H2O2 + k_Srx*PrxSO2 - k_Sulf*PrxSO*H2O2 - k_Cond*PrxSO;
    eqs(2) = Params.k_Ox*PrxS*vars(1) + Params.k_Srx*vars(3) - Params.k_Sulf*vars(2)*vars(1) - Params.k_Cond*vars(2);
    % dPrxSO2/dt
    % eqs(3) = k_Sulf*PrxSO*H2O2 - k_Srx*PrxSO2;
    eqs(3) = Params.k_Sulf*vars(2)*vars(1) - Params.k_Srx*vars(3);
    % dPrxSS/dt
    % eqs(4) = k_Cond*PrxSO - k_Red*TrxSH*PrxSS;
    eqs(4) = Params.k_Cond*vars(2) - Params.k_Red*TrxSH*vars(4);
    % dTrxSS/dt
    % eqs(5) = k_Red*TrxSH*PrxSS - VAppMax*(TrxSS)/(K_M+TrxSS);
    eqs(5) = Params.k_Red*TrxSH*vars(4) - Params.VAppMax*(vars(5))/(Params.K_M+vars(5));


end