%% plot_phenotypes.m
% Script to plot various "phenotypic regions" arising from subsets of the "design space"

%% Obtain species trajectories vs parameter val: v_sup
% global Params
% Params = struct();

% % Initialize 
% v_sups = zeros(7,1);
% steady_states = NaN(7,2);
% for v_sup_exp = -8:-2
%     v_sup_val = 10^(v_sup_exp);
%     v_sups(v_sup_exp+9,1) = v_sup_val;
% 
%     ss = simulate_selvaggio(v_sup_val);
%     steady_states(v_sup_exp+9,:) = ss;
% end

% Initialize 
v_sups = zeros(61,1);
steady_states = NaN(61,3);
for i = 1:61
    v_sup_exp = -8+0.1*(i-1);
    v_sup_val = 10^(v_sup_exp);
    v_sups(i,1) = v_sup_val;

    ss = simulate_selvaggio(v_sup_val);
    steady_states(i,:) = ss;
end


%% v_sup UB for TTPU: 
% a = Params.k_Cond*Params.PrxTotal;
% b = sqrt(Params.k_Cond*Params.k_Ox*Params.k_Srx/Params.k_Sulf)*Params.PrxTotal;
% c = Params.k_Cond*Params.k_Ox*Params.PrxTotal/Params.k_Sulf;
% d = Params.k_Red*Params.PrxTotal*Params.TrxTotal;
% e=Params.VAppMax;
% f=Params.VAppMax*Params.TrxTotal/Params.K_M;
% base = min([a b c d e f]);
% UB = base^Params.k_Alt;
% REMARK: I potentially did this computation wrong but for TTPU with
% default param values I'm getting UB = 3.7e186 which seems unrealisitcally
% high?


%% Plot
figure
loglog(v_sups, steady_states(:,1), 'o')
title("PrxS Steady State versus v_{sup}")
ylim([10^-7 1])
xlim([10^-8 10^-2])

figure
loglog(v_sups, steady_states(:,2), 'o')
title("PrxSO Steady State versus v_{sup}")
ylim([10^-7 1])
xlim([10^-8 10^-2])

figure
loglog(v_sups, steady_states(:,3), 'o')
title("PrxSO2 Steady State versus v_{sup}")
ylim([10^-7 1])
xlim([10^-8 10^-2])
