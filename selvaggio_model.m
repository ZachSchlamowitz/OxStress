selvaggio_model.m

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




function eqs = selvaggio(t,vars)
    % Key:
    % vars = [(1)H2O2, (2)PrxSO, (3)PrxSO2, (4)PrxSS, (5)TrxSS]
    % Note that we back out PrxS from total Prx = PrxS + PrxSO + PrxSS + PrxSO2
    % and similarly we back out TrxSH from total Trx = TrxSH + TrxSS.
    
    eqs = zeros(5,1);
    % dH202/dt
    % eqs(1) = v_sup - k_Alt*H2O2 - k_Ox*PrxS*H2O2 - k_Sulf*PrxSO*H2O2;
    eqs(1) = v_sup - k_Alt*vars(1) - k_Ox*PrxS*vars(1) - k_Sulf*vars(2)*vars(1);
    % dPrxSO/dt
    % eqs(2) = k_Ox*PrxS*H2O2 + k_Srx*PrxSO2 - k_Sulf*PrxSO*H2O2 - k_Cond*PrxSO;
    eqs(2) = k_Ox*PrxS*vars(1) + k_Srx*vars(3) - k_Sulf*vars(2)*vars(1) - k_Cond*vars(2);
    % dPrxSO2/dt
    % eqs(3) = k_Sulf*PrxSO*H2O2 - k_Srx*PrxSO2;
    eqs(3) = k_Sulf**vars(2)*vars(1) - k_Srx*vars(3);
    % dPrxSS/dt
    % eqs(4) = k_Cond*PrxSO - k_Red*TrxSH*PrxSS;
    eqs(4) = k_Cond*vars(2) - k_Red*TrxSH*vars(4);
    % dTrxSS/dt
    % eqs(5) = k_Red*TrxSH*PrxSS - VAppMax*(TrxSS)/(K_M+TrxSS);
    eqs(5) = k_Red*TrxSH*vars(4) - VAppMax*(vars(5))/(K_M+vars(5));


end