%% plot_phenotypes.m
% Script to plot various "phenotypic regions" arising from subsets of the "design space"

%% Obtain species trajectories vs parameter val: v_sup
% global Params
% Params = struct();

% Initialize 
v_sups = zeros(16,1);
for v_sup_exp = -5:10
    v_sup_val = exp(v_sup_exp);
%     Params.v_sup = v_sup_val;
%     selvaggio_model.m
    v_sups(v_sup_exp+6,1) = v_sup_val;
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

