%% selvaggio_LSA.m 
% Author: Zach Schlamowitz, 2/15/23
%
% SELVAGGIO_LSA: Linear Stability Analysis of Selvaggio ODEs
% This script implements a simple linear stability analyis on the single
% Prx species variant of the PTRS model (see text of Selvaggio et al.
% (2018) for the source equations). The script outputs eigenvectors and
% eignevalues corresponding to a given steady state of the PTRS (these
% outputs can be used to identify the stability of said steady state and
% rates of convergence to it). Obtaining the steady state itself can be
% approached in two ways: 
% (1) identification of equilbiria is included as 
% a functionality in this script using symbolic computation; however, this
% may fail to identify all equilibria and we did not thoroughly exlpore its
% accuracy. 
% (2) simply inputting a steady state vector identified elsewhere. For
% instance, one can run a simulation of the PTRS under the desired
% experimental condition and look at state trajectories to identify a
% steady state which the system approaches as time progresses. Saving this
% steady state in the workspace (as a vector) and then inputting it in the
% appropriate place below (see comments) allows the same LSA as does 
% computing the steady state using this script. 
%
% Note that option (1) requires installation of the MATLAB Symbolic 
% Computation Toolbox. 

%% Setup
% Store default parameters for 1-prx model
Params = struct();
Params.v_sup = 10^(-5+6); % (µM/sec) (converted from M/sec) ranges 10^-7 to 10^-3 in Fig 3, based on varying levels of extracellular H202
Params.k_Alt = 79; % (sec^-1) default: 79
Params.k_Ox = 40; % (µM^-1 sec^-1) (converted from 4e7 M^-1 sec^-1)
Params.k_Sulf = 5.1e-3; % (µM^-1 sec^-1) (converted from 5.1e3 M^-1 sec^-1)
Params.k_Srx = 3.3e-3; % (sec^-1) (converted from 3.3 10^-3 sec^-1)
Params.k_Cond = 6.4; % (sec^-1) default: 6.4
Params.k_Red = 0.21; % (µM^-1 sec^-1) (converted from 2.1e5 M^-1 sec^-1)
Params.VAppMax = 230; % (µM/sec) (converted from 0.23 mM/sec) NOTE that is reported as VMax in Table 2; presumably interchangeable?
Params.K_M = 1.8; % (µM)
Params.PrxTotal = 92; % total concentration?? FLAG of Prx (µM)
Params.TrxTotal = 23; % total concentration?? FLAG of Trx (µM)

%% Symbolic Calculation of Steady State(s)
% Declare the state variables of the system as symbolic variables
syms H202 PrxSO PrxSO2 PrxSS TrxSS

% Note that we back out PrxS from total PrxTotal = PrxS + PrxSO + PrxSS + PrxSO2
PrxS = Params.PrxTotal - PrxSO - PrxSO2 - PrxSS;
% and similarly we back out TrxSH from total TrxTotal = TrxSH + TrxSS.
TrxSH = Params.TrxTotal - TrxSS;

% Store system equations symbolically...
F = [Params.v_sup - Params.k_Alt*H202 - Params.k_Ox*PrxS*H202 - Params.k_Sulf*PrxSO*H202, ...
    Params.k_Ox*PrxS*H202 + Params.k_Srx*PrxSO2 - Params.k_Sulf*PrxSO*H202 - Params.k_Cond*PrxSO, ...
    Params.k_Sulf*PrxSO*H202 - Params.k_Srx*PrxSO2, ...
    Params.k_Cond*PrxSO - Params.k_Red*TrxSH*PrxSS, ...
    Params.k_Red*TrxSH*PrxSS - Params.VAppMax*(TrxSS)/(Params.K_M+TrxSS)];

% ... and state vector
v = [H202 PrxSO PrxSO2 PrxSS TrxSS];

% Compute Jacobian matrix of system
tic
J = jacobian(F,v)
toc

x0 = [0.01; 0.01; 0.001; 0.01; 0.01]; % initial values
fun = @(vars)selvaggio(vars,Params); % anonymous function for the model equations
EQa = fsolve(fun, x0, Params); % Note that there's a way to get it to output jacobian here also, obviating need for jacobian computation above
% NOTE: For option (2), simply replace this computation with an upload of
% your pre-computed steady state equilibrium vector.

%% Perform Linear Stability Analysis using EQa and Jacobian
% Evaluate jacobian at the steady state equilibrium (EQa)
JatEQa = subs(J,v,EQa');
JatEQa = double(JatEQa);

% Get eigenvalues and eigenvectors
[evec, eval] = eig(JatEQa)

% Now go analyze!!
% -----------------------------------------------------------------------

%% Script Functions
% Model
function eqs = selvaggio(vars, Params)
    % global Params
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