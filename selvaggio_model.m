% selvaggio_model.m

% This is an implementation of the Selvaggio (2018) PRX system ODE model

%% Params, Commands, Etc.

% Initialize parameters
% global Params


% disp(strcat('Time elapsed: ', num2str(toc)))


%% Plotting

% % Time Series Plot
% plot(time, sol, LineWidth=2)
% % title('Predator, Prey Populations Over Time')
% xlabel('t')
% % ylabel('Population')
% legend('H2O2', 'PrxSO', 'PrxSO2', 'PrxSS', 'TrxSS', 'PrxS','TrxSH', 'Location', 'North', 'FontSize',14 )
% 
% figure
% 
% beginning = ceil(.75*length(time));
% plot(time(1:beginning), sol(1:beginning, :), LineWidth=2)
% % title('Predator, Prey Populations Over Time')
% xlabel('t')
% % ylabel('Population')
% legend('H2O2', 'PrxSO', 'PrxSO2', 'PrxSS', 'TrxSS','PrxS','TrxSH','Location', 'North', 'FontSize',14 )

%% Model
function eqs = selvaggio_model(t,vars, Params)
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