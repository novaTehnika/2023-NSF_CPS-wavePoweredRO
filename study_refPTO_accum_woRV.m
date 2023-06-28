% study_refPTO_accum_woRV.m script m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 6/28/2023
%
% PURPOSE/DESCRIPTION:
% This script performs parameter variation studies
% using the model contained in sys_refPTO.m and solved by
% sim_refPTO.m.
% The parameter initiallization functions are called within this
% script before the sim_refPTO.m script is called.
%
% This specific script studies the size of high-pressure accumulator.
%
% This script is set up to be run as part of a SLURM job array. The
% following lines are required before this script is called:
%   iVar = ${SLURM_ARRAY_TASK_ID};
%   addpath('Utilities'); 
%   startParPool($SLURM_JOB_CPUS_PER_NODE);
%
% FILE DEPENDENCY:
% ./Reference PTO/
%   initialConditionDefault_refPTO
%   parameters_refPTO.m
%   sim_refPTO.m
%   stateIndex_refPTO.m
%   sys_refPTO.m
% ./WEC model/
%   flapModel.m
%   hydroStaticTorque.m
%   parameters_WECmodel.m
% ./WEC model/WECdata
%   nemohResults_vantHoff2009_20180802.mat
%   vantHoffTFCoeff.mat
% ./Solvers/
%   deltaE_NI.m
%   deltaV_NI.m
%   ode1.m
% ./Components/
%   areaFracPWM.m
%   capAccum.m
%   deadVCap.m
%   flowCV.m
%   flowPRV.m
% ./Utilities/
%   statsTimeVar_cdf.m
%
% UPDATES:
% 6/28/2023 - Created from study_refPTO_accum_wActiveRV.m
%
% Copyright (C) 2023  Jeremy W. Simmons II
% 
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program. If not, see <https://www.gnu.org/licenses/>.
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear
% clc
addpath('WEC model') 
addpath(['WEC model' filesep 'WECdata']) 
addpath('Reference PTO')
addpath('Components')
addpath('Sea States')
addpath('Solvers')
addpath('Utilities')
%% %%%%%%%%%%%%   SIMULATION PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation timeframe
par.Tramp = 250; % [s] excitation force ramp period
par.tstart = 0; %[s] start time of simulation
par.tend = 500; %[s] end time of simulation

% Solver parameters
% par.odeSolverRelTol = 1e-4; % Rel. error tolerance parameter for ODE solver
% par.odeSolverAbsTol = 1e-4; % Abs. error tolerance parameter for ODE solver
par.MaxStep = 5e-5;             % [s] for fixed time solver, this is the step size for solver
par.downSampledStepSize = 1e-2; % [s] specifies time step for data output
if mod(par.downSampledStepSize,par.MaxStep)
    warning('down-sampled time step is not an integer multiple of the maximum step size')
end

% Sea State and Wave construction parameters
% Hs = [2.34 2.64 5.36 2.05 5.84 3.25];
% Tp = [7.31 9.86 11.52 12.71 15.23 16.5];
par.wave.Hs = 2.64;
par.wave.Tp = 9.86;
par.WEC.nw = 1000; % num. of frequency components for harmonic superposition
par.wave.rngSeedPhase = 3; % seed for the random number generator

% load parameters
par = parameters_refPTO(par,...
    'nemohResults_vantHoff2009_20180802.mat','vantHoffTFCoeff.mat');

% Define initial conditions
stateIndex_refPTO % load state indices, provides 'iy_...'
initialConditionDefault_refPTO % default ICs, provides 'y0'

%% Special modifications to base parameters
% par.Sro = 3000; % [m^3]
% par.D_WEC = 0.3;         % [m^3/rad] flap pump displacement
% par.control.p_ro_nom = 7e6; % [Pa]

% par.ERUconfig = 1;

par.rvIncluded = 0; % RO inlet valve is 1 - present, 0 - absent
par.rvConfig = (0)*par.rvIncluded; % RO inlet valve is 1 - active, 0 - passive
dp_rated = 1e5; % [Pa] 
q_rated = 1000e-3; % [(lpm) -> m^3/s]
par.kv_rv = q_rated/dp_rated;

par.duty_sv = 0;

%% %%%%%%%%%%%%   Study Variables  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% total accumulator volume

nVar = 15;
Vtotal = logspace(log10(1e3),log10(50e3),nVar);% [m^3] total accumulator volume

saveSimData = 1; % save simulation data (1) or just output variables (0)

%% %%%%%%%%%%%%   COLLECT DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% change design parameter
par.Vc_h = 0.5*Vtotal(iVar);
par.Vc_ro = 0.5*Vtotal(iVar);

% run simulation
ticSIM = tic;
out = sim_refPTO(y0,par);
toc(ticSIM)

% Calculate metrics
it_vec = find(out.t>=par.tstart);
% max rate of change in pressure
% 97th percentile ratof change
% power loss from valve
% power loss through PRVs
% permeate production
% 
q_permMean = mean(out.q_perm(it_vec));
PP_WEC = mean(out.power.P_WEC(it_vec));
PP_wp = mean(out.power.P_wp(it_vec));
PP_rv = mean(out.power.P_rv(it_vec));
PP_hPRV = mean(out.power.P_hPRV(it_vec));
PP_roPRV = mean(out.power.P_roPRV(it_vec));
PP_aPRV = mean(out.power.P_bPRV(it_vec));
PP_bPRV = mean(out.power.P_bPRV(it_vec));
dpdt_max = max(abs(out.y(it_vec,iyp_ro)));

dist_dpdt = statsTimeVar_cdf(out.t(it_vec),abs(out.y(it_vec,iyp_ro)));
dpdt_97 = dist_dpdt.xi(find(dist_dpdt.f > 0.97,1,'first'));

if ~saveSimData
    clear out
end



%% %%%%%%%%%%%%   Save Data  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeStamp = datestr(now,'yyyymmddTHHMMSS'); % time in ISO8601

% Save data
filename = ['data_refPTO_accum_wPassiveRV_', ...
            datestr(now,'yyyymmdd'),'_', ...
            num2str(iVar)];
save(filename)

%% %%%%%%%%%%%%   End Computations  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return

%% %%%%%%%%%%%%   PLOTTING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% WEC pump efficiency
X = kv*sqrt(1000);
Y = Tsw;
Z = eff_wecPump_2D; 
inc = 250;
levels = round(max(eff_wecPump_2D(:)) - logspace(log10(0.01),log(0.9),8),2);

figure
C1 = contour(log10(X),Y,Z,levels,'-k','ShowText','on');
xt = xticks;
xticklabels("10^{"+string(xt')+"}");

xlabel('flow coefficient (L/s/kPa^{1/2})');
ylabel('switching period (s)');
title('WEC-driven Pump Efficiency')

% Absorbed power
X = kv*sqrt(1000);
Y = Tsw;
Z = PP_WEC_2D*1e-3;


figure
C2 = contour(log10(X),Y,Z,'-k','ShowText','on');
xt = xticks;
xticklabels("10^{"+string(xt')+"}");

xlabel('flow coefficient (L/s/kPa^{1/2})');
ylabel('switching period (s)');
title('WEC Power Absorption (kW)')

% Pump output power
X = kv*sqrt(1000);
Y = Tsw;
Z = PP_wp_2D*1e-3;


figure
C2 = contour(log10(X),Y,Z,'-k','ShowText','on');
xt = xticks;
xticklabels("10^{"+string(xt')+"}");

xlabel('flow coefficient (L/s/kPa^{1/2})');
ylabel('switching period (s)');
title('WEC-Driven Pump Power Output (kW)')
