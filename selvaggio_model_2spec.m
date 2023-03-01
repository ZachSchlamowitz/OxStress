% selvaggio_model_2spec.m

% This is an implementation of the Selvaggio (2018) PRX system ODE model
% from their supplement that additionally models two species of PRX and
% accounts for permeability

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
function eqs = selvaggio_model_2spec(t,vars, Params)
    global Params
    % Key:
    % vars = [(1)H2O2_out, (2)H2O2, (3)PrxISO, (4)PrxISO2, (5)PrxISS, (6)PrxIISO, (7)PrxIISO2, (8)PrxIISS, (9)TrxSS]
    % Note that we back out PrxS from total PrxTotal = PrxS + PrxSO + PrxSS + PrxSO2
    PrxIS = Params.PrxITotal - vars(3) - vars(4) - vars(5);
    PrxIIS = Params.PrxIITotal - vars(6) - vars(7) - vars(8);
    % and similarly we back out TrxSH from total TrxTotal = TrxSH + TrxSS.
    TrxSH = Params.TrxTotal - vars(9);

    eqs = zeros(9,1);
    % dH202_out/dt
    eqs(1) = (Params.kappa*Params.n_cells*Params.S/Params.V_medium)*(vars(2)-vars(1));
    
    % dH2O2/dt
    eqs(2) = (Params.kappa*Params.S/Params.V_cytoplasm)*()
    
    Params.v_sup - Params.k_Alt*vars(1) - Params.k_Ox*PrxS*vars(1) - Params.k_Sulf*vars(2)*vars(1);
    
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