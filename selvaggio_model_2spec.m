%% selvaggio_model_2spec.m
% Author: Zach Schlamowitz
%
% SELVAGGIO_MODEL_2SPEC is a modified implementation of the Selvaggio et 
% al. (2018) Prx system ODE model that models two species of PRX (but
% does not account for permeability). This function houses the ODEs themselves in an 
% isolated file that can be passed in (via a function handle) to an ODE 
% solver. Such usage is core to the simulation process run using either 
% plot_ptrs.m > simulate_selvaggio.m or simulate_selvaggio_full_single.m.
%
% Source equations for the model are modified from the supplementary 
% material to Selvaggio et al. (2018); in this implementation, there is 
% only one H2O2 compartment (cytoplasmic), resulting in removal of the 
% first equation and modification of the second from their presentation. 
% For the two Prx species model with H2O2 membrane permation, see 
% selvaggio_model_2spec_perm.m, and for the one Prx species model, see 
% selvaggio_model.m.

%% Model
function eqs = selvaggio_model_2spec(t,vars, Params)
    global Params

    % Key:
    % vars = [(1)H2O2_out, (2)H2O2, (3)PrxISO, (4)PrxISO2, (5)PrxISS, (6)PrxIISO, (7)PrxIISO2, (8)PrxIISS, (9)TrxSS]
    % Note that we back out PrxS from total PrxTotal = PrxS + PrxSO + PrxSS + PrxSO2
    PrxIS = Params.PrxITotal - vars(2) - vars(3) - vars(4);
    PrxIIS = Params.PrxIITotal - vars(5) - vars(6) - vars(7);
    % and similarly we back out TrxSH from total TrxTotal = TrxSH + TrxSS.
    TrxSH = Params.TrxTotal - vars(8);

    eqs = zeros(8,1);    
    % dH2O2/dt
    eqs(1) = Params.v_sup - Params.k_Alt*vars(1) - Params.k_Ox*(PrxIS + PrxIIS)*vars(1)...
        - Params.k_Sulf_I*vars(2)*vars(1) - Params.k_Sulf_II*vars(5)*vars(1);
        
    % dPrxI-SO/dt
    eqs(2) = Params.k_Ox*PrxIS*vars(1) + Params.k_Srx*vars(3) - Params.k_Sulf_I*vars(2)*vars(1) - Params.k_Cond_I*vars(2);

    % dPrxI-SO2/dt
    eqs(3) = Params.k_Sulf_I*vars(2)*vars(1) - Params.k_Srx*vars(3);

    % dPrxI-SS/dt
    eqs(4) = Params.k_Cond_I*vars(2) - Params.k_Red*TrxSH*vars(4);
    
    % dPrxII-SO/dt
    eqs(5) = Params.k_Ox*PrxIIS*vars(1) + Params.k_Srx*vars(6) - Params.k_Sulf_II*vars(5)*vars(1) - Params.k_Cond_II*vars(5);

    % dPrxII-SO2/dt
    eqs(6) = Params.k_Sulf_II*vars(5)*vars(1) - Params.k_Srx*vars(6);

    % dPrxII-SS/dt
    eqs(7) = Params.k_Cond_II*vars(5) - Params.k_Red*TrxSH*vars(7);
    
    % dTrxSS/dt
    eqs(8) = Params.k_Red*TrxSH*(vars(4) + vars(7)) - Params.VAppMax*(vars(8))/(Params.K_M+vars(8));

end