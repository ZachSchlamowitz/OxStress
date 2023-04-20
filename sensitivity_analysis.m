%% sensitivity_analysis.m
% Author: Zach Schlamowitz, 4/19/2023
% This script performs sensitivity analysis on the parameters of the full 
% Selvaggio model of the Prdx System. The intuition of the following simple
% algorithm is to twiddle parameters by order of magnitudes (* or / by 2 or
% 10) and see how the resulting "steady state" of the system changes. We
% define the steady state to be just the resting state of the system after
% 30-40 minutes post-bolus. (Note that it may be of value to consider how
% changing the parameters affects the first 10 min or so of the simulation,
% even if the long-run "steady states" don't change.) Then, we compute the
% change in the position of the steady state vector under deformations of
% the parameters by computing the Euclidean distance between the ss vector
% before and after the deformations.

% The output is presented in a matrix:
% Parameter | Distance(default ss, ss with param/2) | Distance(default ss, ss with param*2)  
% -----------------------------------------------------------------------------------------
%   k_Alt   |                 10.554                |                   8.554  
%   k_Ox    |                   0.45                |                    0.20   
%                                                  ETC. 
tic 
cell_type = 'HEK293';

default_params = struct();

default_params.v_sup = 10^(-5+6); % (µM/sec) (converted from M/sec) ranges 10^-7 to 10^-3 in Fig 3, based on varying levels of extracellular H202

if cell_type == 'HEK293'
    % From Table 2
    default_params.k_Alt = 160; % (sec^-1) default: 1.6e2 (source: Table 2)
    default_params.k_Ox = 40; % (µM^-1 sec^-1) (converted from 4e7 M^-1 sec^-1) 
    default_params.k_Srx = 4.1e-4; % (sec^-1) (converted from 0.41 10^-3 sec^-1)
    default_params.k_Red = 0.21; % (µM^-1 sec^-1) (converted from 2.1e5 M^-1 sec^-1)
    default_params.VAppMax = 190; % (µM/sec) (converted from 0.19 mM/sec) NOTE that is reported as VMax in Table 2; presumably interchangeable?
    default_params.K_M = 1.8; % (µM)
    default_params.PrxTotal = 140; % total concentration of Prx (µM)
    default_params.TrxTotal = 46; % total concentration of Trx (µM) (same value in ST6)
    default_params.k_Sulf = 3.7e-3; % (µM^-1 sec^-1) (converted from 3.7e3 M^-1 sec^-1)
    default_params.k_Cond = 7.3; % (sec^-1) default: 7.3

    % From Supplement section 3.2.3 (text)
    default_params.k_Sulf_I = 1.3e-3; % (µM^-1 sec^-1) (converted from 1.3e3 M^-1 sec^-1)
    default_params.k_Sulf_II = 1.2e-2; % (µM^-1 sec^-1) (converted from 1.2e4 M^-1 sec^-1)
    
    default_params.k_Cond_I = 9; % (sec^-1) default: 11 [default from Armindo's email Table with cite]
    default_params.k_Cond_II = 1.7; % (sec^-1) default: 0.5 [default from Armindo's email Table with cite]
    
    % From Supplementary Table 6  
    default_params.PrxITotal = 110; % total concentration of PrxI (µM)
    default_params.PrxIITotal = 32; % total concentration of PrxII (µM)

    % Permeation Parameters (various sources)
    default_params.n_cells = 3e5; % number of cells in media / well of plate; default: 3e5 cells per well from Sobotta et al 2013
    default_params.kappa = 15; % permeability coefficient (µm/sec) in erythrocytes (source: 3/2023 meeting with Armindo)
    default_params.V_medium = 2e11; % volume of media (µm^3) default: 2e11; converted to µm^3 from 200 µL (source: guess from Andrew)
    default_params.V_cytoplasm = 1150; % volume of cell (µm^3); default: 660 from Supplement 3.2.3, although this may come from Jurkat T cells, which are smaller; alternatively, guess 1150 based on average HEK293 cell diameter of 13 µm (source: https://bionumbers.hms.harvard.edu/bionumber.aspx?id=108893&ver=3&trm=HEK+293&org= among others)
    default_params.S = 530; % surface area of cell (µm^2) guess based on average HEK293 cell diameter of 13 µm (see above)

elseif cell_type == 'MCF7  '
    % From Table 2
    default_params.k_Alt = 79; % (sec^-1) default: 79
    default_params.k_Ox = 40; % (µM^-1 sec^-1) (converted from 4e7 M^-1 sec^-1)
    default_params.k_Srx = 3.3e-3; % (sec^-1) (converted from 3.3 10^-3 sec^-1)
    default_params.k_Red = 0.21; % (µM^-1 sec^-1) (converted from 2.1e5 M^-1 sec^-1)
    default_params.VAppMax = 230; % (µM/sec) (converted from 0.23 mM/sec) NOTE that is reported as VMax in Table 2; presumably interchangeable?
    default_params.K_M = 1.8; % (µM)
    default_params.k_Sulf = 5.1e-3; % (µM^-1 sec^-1) (converted from 5.1e3 M^-1 sec^-1)
    default_params.k_Cond = 6.4; % (sec^-1) default: 6.4
    default_params.PrxTotal = 92; % total concentration of Prx (µM) (NOTE: total = Params.PrxITotal + Params.PrxIITotal if use ST6 values)
    
    % From Table in Armindo's Email: (formerly from Supplementary Table 6 or Supplement section 3.2.3 (text))
    default_params.TrxTotal = 20; % total concentration of Trx (µM) default: 20; updated: was 23 from Supplement Table 6, now from Armindo's email
    default_params.PrxITotal = 110; % total concentration of PrxI (µM) default: 110; updated: was 59 from Supplement Table 6, now from Armindo's email
    default_params.PrxIITotal = 30; % total concentration of PrxII (µM) default: 30; updated: was 33 from Supplement Table 6, now from Armindo's email
    
    default_params.k_Sulf_I = 1.5e-3; % (µM^-1 sec^-1) default: 1.5e-3 [source: Table from Armindo's email]
    default_params.k_Sulf_II = 3.4e-3; % (µM^-1 sec^-1) default: 3.4e-3 [source: Table from Armindo's email]
    
    default_params.k_Cond_I = 11; % (sec^-1) default: 11 [see Armindo's email with cite]
    default_params.k_Cond_II = 0.5; % (sec^-1) default: 0.5 [see Armindo's email with cite]

    % Permeation Parameters (various sources)
    default_params.n_cells = 6000; % default: 6000 source: guess by Andrew
    default_params.kappa = 15; % permeability coefficient (µm/sec) in erythrocytes (source: 3/2023 meeting with Armindo)
    default_params.S = 1224.18; % surface area of cell (µm^2) guessed value based on SA of sphere of diamter 19.74 µm, https://www.researchgate.net/figure/Size-distribution-of-MCF-7-cells-used-in-this-study-The-histogram-was-derived-using-an_fig5_257966698
    default_params.V_medium = 2e11; %? default: 2e11 µm^3; converted to µm^3 from 200 µL (source: guess from Andrew)
    default_params.V_cytoplasm = 1760; %? µm^3 based on MCF7 cell vol of 1.76pL from https://bionumbers.hms.harvard.edu/bionumber.aspx?id=115154&ver=1
end

% Get the default state vector
[ss_default, hm_v_default] = simulate_selvaggio_full_single(1,NaN,NaN,'HEK293',default_params);
vec_default = ss_default(:,2100); % take state vector at 35min = 2100sec

% Now make new parameter struct that will be used in the sensitivity analysis sims
Params = default_params;
global Params

% Get container for all the parameters' names
fn = fieldnames(Params);

% Initialize sensitivity analysis output matrix
SAM = NaN(numel(fn),4); % 4 columns correspond to /10, /2, *2, *10 distances
for k=1:numel(fn) % loop over struct fields by fieldname
    %% Twiddle parameter down (/10)
    Params.(fn{k}) = default_params.(fn{k}) / 10;
    
    % Simulate with this parameter set
    [ss_down_10, hm_v_down_10] = simulate_selvaggio_full_single(1,NaN,NaN,'HEK293',Params);
    
    % Store output "steady state" vector
    vec_down_10 = ss_down_10(:,2100); % take state vector at 35min = 2100sec

    %% Twiddle parameter down (/2)
    Params.(fn{k}) = default_params.(fn{k}) / 2;
    
    % Simulate with this parameter set
    [ss_down_2, hm_v_down_2] = simulate_selvaggio_full_single(1,NaN,NaN,'HEK293',Params);
    
    % Store output "steady state" vector
    vec_down_2 = ss_down_2(:,2100); % take state vector at 35min = 2100sec

    %% Twiddle parameter up (*2)
    Params.(fn{k}) = default_params.(fn{k}) * 2;
    
    % Simulate with this parameter set
    [ss_up_2, hm_v_up_2] = simulate_selvaggio_full_single(1,NaN,NaN,'HEK293',Params);
    
    % Store output "steady state" vector
    vec_up_2 = ss_up_2(:,2100); % take state vector at 35min = 2100sec
    
    %% Twiddle parameter up (*10)
    Params.(fn{k}) = default_params.(fn{k}) * 10;
    
    % Simulate with this parameter set
    [ss_up_10, hm_v_up_10] = simulate_selvaggio_full_single(1,NaN,NaN,'HEK293',Params);
    
    % Store output "steady state" vector
    vec_up_10 = ss_up_10(:,2100); % take state vector at 35min = 2100sec
    
    %% Compute distances to default vector 
    SAM(k,1) = norm(vec_down_10 - vec_default);
    SAM(k,2) = norm(vec_down_2  - vec_default);
    SAM(k,3) = norm(vec_up_10   - vec_default);
    SAM(k,4) = norm(vec_up_2    - vec_default);

    %% Reset this parameter
    Params.(fn{k}) = default_params.(fn{k});
end
toc
