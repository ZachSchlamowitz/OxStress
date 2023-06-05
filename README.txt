README.txt
Project: EXAMINING THE LINK BETWEEN PEROXIREDOXIN PROTEINS AND MUTUALLY EXCLUSIVE 
         TRANSCRIPTION FACTOR ACTIVATION WITH A MATHEMATICAL MODEL
Author: Zach Schlamowitz, University of Arizona, zschlamowitz at gmail dot com
Date: June 2023

This is the README file for Zach Schlamowitz's implementation of a mathematical model of the 
peroxiredoxin-thioredoxin system (PTRS) presented by Selvaggio et al. (2018) in their paper
"Mapping the phenotypic repertoire of the cytoplasmic 2-Cys peroxiredoxin-thioredoxin system 
The model is a system of nonlinear ordinary differential equations (ODEs) and allows description 
of one or two species of peroxiredoxin (prx). (In theory, any number of species of prxs could 
be modeled by just adding copies of corresponding sets of equations.) This implementation was 
originally part of Zach Schlamowitz's Senior Honors Thesis project at the University of Arizona 
(2023), entitled "Examining the Link Between Peroxiredoxin Proteins and Mutually Exclusive 
Transcription Factor Activation with a Mathematical Model," and was used to examine connections 
between cellular oxidative stress response (namelly the peroxiredoxin signaling / stress response 
system) and the dynamics of related transcription factor activity. A copy of this manuscript can
be found here: [LINK].

Below, we outline the key code structure and provide instructions for using the model to simulate
the PTRS under different H2O2 conditions.

-------------------------------------------
CODE STRUCTURE:
A brief note about how this code works:
To simulate the PTRS system described by the ODEs of the Selvaggio et al. model, we use a standard 
built-in MATLAB numerical ODE solver (e.g., ode23s). After passing this solver our system of ODEs 
(among other function parameters), the solver returns a numerical approximation to a solution of the ODEs 
in discrete form as a vector of timepoints and a corresponding vector of solution values. In this way, 
we obtain a system of discrete mathematical equations that describe the changes in concentrations of 
components of the PTRS as functions of time. Then, we can observe the simulated effects of an
experimental condition by looking at the trajectories of these solution values. In other words, we use
the following pipeline:

experimental conditions (e.g., H2O2 bolus)  -->  ODE u'(t) = f(u(t))  -->  ode23s  --> u(t) in discrete 
form (u(1)=x1, u(2)=x2, etc.) -->  observe predicted outcomes

To see the 1-prx and 2-prx ODE equations that are passed into these solvers, please refer to the
following sources: 
- for 1-prx equations: see the text of Selvaggio et al. (2018)
- for 2-prx equations: see either the text of Schlamowitz (2023) thesis document or 
                       the supplementary material for Selvaggio et al. (2018)

However, the way MATLAB designed these ODE solvers (e.g., ode23s) is such that when calling them, one 
has to pass in a function handle (e.g., "@myfunc") that refers to a separate script file. This separate
script file must *only* contain the system of ODEs to be solved, written in prpoper format for the solver
(for e.g., ode23s), and continaing nothing else. As a result, we did not see a way to directly pass
model parameters (e.g., k_Red) into the solver (since we did not wish to hard-code parameter values in the
model equations in their separate script). To troubleshoot this problem, we used the following workaround; 
perhaps a better one is known to better programmers. (If you are this programmer, please feel free to contact 
Zach Schlamowitz). We make a parameters struct, Params, which we declare as a global variable; this way, 
parameters do not have to be passed in to the equations script diectly. 

Now, to run one simulated condiiton requires calling the solver (e.g., ode23s) on our PTRS equations
once. However, for plotting purposes, we often wish to examine simulated timecourses under various 
conditions side by side, e.g., following various doses of bolus admisssion. As such, we wrap another
script around that which contains the equations; this wrapper file evolved into simulate_ptrs.m, the 
main file for running this project. 

-------------------------------------------
HOW TO USE:
Usage of this code should almost exclusively require simply running either < plot_ptrs.m > or 
< simulate_selvaggio_full_single.m >. For simulations and figure replication, run < plot_ptrs.m >. The code 
has been designed so that all figures shown in the thesis document can be gerenated by running this file with 
little to no modification. See below for details on how to generate figures. If you instead wish to simply 
simulate the model under one condition (e.g., one specific value of a parameter instead of a range of values), 
use the single-simulation script < simulate_selvaggio_full_single.m > (see below for more instructions).

FULL SIMULATION (plot_ptrs): This will generate the figures from Schlamowitz (2023). 
***************
< plot_ptrs > has a designated user options section at the top, containing: 
(1) a switch for 1-prx vs 2-prx model
(2) the vector of extracellular H2O2 bolus values that will be used in simulations of the 2-prx model
(3) corresonding bolus labels (make sure there is one label for each bolus value)
(4) a switch for setting the cell type (MCF7 or HEK293)
Note that (2)-(4) only have effect if (1) has been set to use the 2-prx model. 
[Note also that other cell types could be added by adding their parameters to < simulate_selvaggio.m  >
and adding a corresponding value to the swtich and if-else logic.]

For 1-prx model: set num_prx =1.
The user may modify the parameter intracellular_val if they wish to change the starting concentration of H2O2.
Note that although in this model variant an H2O2 supply rate, v_sup, is necesary, the user does not need to specify
it. This is because the loop in < plot_ptrs > iterates over values of v_sup. The < selvaggio_model > parameters
"bolus" and "cell_type" are superfluous to this model variant and can be left as NaNs. 

For 2-prx model: set num_prx =2. 
No further modifications are necessary beyond the User Options section, as the code will loop over the inputted
bolus values and simulate the 2-prx model following extracellular bolus addition of those concentrations of H2O2. 
In this model variant, the v_sup parameter is not used, and may be left as a NaN. Contrastingly, all the parameters
set in the User Options section will have effect. 

FIGURES: 
Running < plot_ptrs > with num_prx = 1 will generate three replication figures, recapitulating 
results from Selvaggio et al. (2018) by plotting steady state values of model states (e.g., PrxS)
vs values of v_sup. These are Fig.4a-c of Schlamowitz (2023).

Running < plot_ptrs > with num_prx = 2 will generate (up to) fifteen figures:
(1) replication figure Fig.4D of Schlamowitz (2023)
(2-3) heatmaps showing prx hyperoxidation vs dose and time, Fig.5 of Schlamowitz (2023)
(4) dose response curve showing prx hyperoxidation vs dose only, Fig.6 of Schlamowitz (2023)
(5-15) disulfide prx timecourses, used in Fig.7 of Schlamowitz (2023)

To generate the prx hyperoxidation curves corresponding to prxI/II knockouts shown in Fig.8
of Schlamowitz (2023), we essentially generated copies of Fig.6 of Schlamowitz (2023), using
the command [>> hold on ] to overlay one on top of the next. This also requires commenting out 
/ uncommenting initial values and model parameters corresponding to those settings. See the code 
comments to find these sections.

***************
SINGLE SIMULATIONS (simulate_selvaggio_full_single): This will run a simulation under a single experimental condition.
This script is largely the same as the simulation script < simulate_selvaggio > called in the full simulation,
but with the user input capability of < plot_ptrs >. Simply set the user-specified parameters in the "User Options" 
section as you like (they are the same as in < plot_ptrs >), set up desired plots in the "Plotting" section at the 
bottom of the file and hit run! This script is largely designed for playful exploration and experimentation, so 
the included plots are very minimal.

-------------------------------------------
OTHER FUNCTIONALITY:

The main uses we expect users to seek are available in the two main sets of file(s) described above, denoted
below with symbols (*) and (o). However, additional functionality was included in the design of the project, 
and is available in the project directory on github as well. Here we provide a brief description of each of 
the files in the github repository to faciliate navigation of them. 

Main functions:
- (*) plot_ptrs: main function to run full (multi-condition) simulations
- (o) simulate_selvaggio_full_single: main function to run a single condition simulation
Susidiary functions:
- (*) simulate_selvaggio: internal function called by < plot_ptrs > to simulate a single condition
- (*/o) selvaggio_model_2spec_perm: function containing ODEs for the two peroxiredoxin model with membrane H2O2 permeation
- (*/o) selvaggio_model: function containing ODEs for the one peroxiredoxin model (no H2O2 permeation)
Other Files:
- selvaggio_model_2spec: function containing ODEs for a variant of the two peroxiredoxin model with no H2O2 permeation 
    (use with caution; we did not explore this model variant to confirm soundness or accuracy)
- selvaggio_LSA: performs linear stability analysis on the 1-prx model using symbolic logic 
    (requires MATLAB symbolic math toolbox)
- sensitivity_analysis: script to test sensitivity of dependence of model states on parameter values

This concludes the README for Schlamowitz's (2023) Honors Thesis Project. 