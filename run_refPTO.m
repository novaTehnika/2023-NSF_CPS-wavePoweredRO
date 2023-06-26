% run_refPTO.m script m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 6/12/2023
%
% PURPOSE/DESCRIPTION:
% This script serves as a shell for running a single simulation
% using the model contained in sys_refPTO.m and solved by 
% sim_refPTO.m.
% The parameter initiallization functions are called within this
% script before the sim_refPTO.m script is called.
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
%
% UPDATES:
% 6/12/2023 - Created from run_parallelPTO.m.
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
clear
% clc
addpath('WEC model') 
addpath(['WEC model' filesep 'WECdata']) 
addpath('Reference PTO')
addpath('Components')
addpath('Sea States')
addpath('Solvers')
%% %%%%%%%%%%%%   SIMULATION PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation timeframe
par.Tramp = 250; % [s] excitation force ramp period
par.tstart = 0; %[s] start time of simulation
par.tend = 2000; %[s] end time of simulation

% Solver parameters
% par.odeSolverRelTol = 1e-4; % Rel. error tolerance parameter for ODE solver
% par.odeSolverAbsTol = 1e-4; % Abs. error tolerance parameter for ODE solver
par.MaxStep = 5e-5;             % [s] for fixed time solver, this is the step size for solver
par.downSampledStepSize = 1e-2; % [s] specifies time step for data output
if mod(par.downSampledStepSize,par.MaxStep)
    warning('down-sampled time step is not an integer multiple of the maximum step size')
end

% Sea State and Wave construction parameters
par.wave.Hs = 2.75;
par.wave.Tp = 12;
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

% par.rvIncluded = 1; % RO inlet valve is 1 - present, 0 - absent
% par.rvConfig = (1)*par.rvIncluded; % RO inlet valve is 1 - active, 0 - passive
% dp_rated = 1e5; % [Pa] 
% q_rated = (100)*60/1e3; % [(lpm) -> m^3/s]
% par.kv_rv = q_rated/dp_rated;

par.duty_sv = 0.0;

%% %%%%%%%%%%%%   COLLECT DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ticSIM = tic;
out = sim_refPTO(y0,par);
toc(ticSIM)

%% %%%%%%%%%%%%   PLOTTING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% WEC pos., vel., and Torque
bottomEdge = 1;
leftEdge = 2;
width = 7.5; % one column: 3+9/16, two column: 7.5
height = 6;
fontSize = 8;
lineWidth = 1;

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];


ax(1) = subplot(3,1,1);
plot(out.t,out.waveElev)
xlabel('time (s)')
ylabel('elevation (m)')
title('Wave Elevation')

ax(2) = subplot(3,1,2);
xlabel('time (s)')
hold on

yyaxis left
plot(out.t,out.theta)
xlabel('time (s)')
ylabel('position (rad)')
% ylim([-pi/2 pi/2])

yyaxis right
plot(out.t,out.theta_dot)
ylabel('angular velocity (rad/s)')
% ylim(10*[-pi/2 pi/2])

ax(3) = subplot(3,1,3);
hold on
plot(out.t,1e-6*out.T_wave)
plot(out.t,1e-6*out.T_pto)
plot(out.t,1e-6*out.T_rad)
plot(out.t,1e-6*out.T_hydroStatic)
ylabel('Torque (MNm)')

legend('T_{w}','T_{PTO}','T_{r}','T_{h}')
      
linkaxes(ax,'x');

sgtitle('WEC Behaviour')

%% Pressure
bottomEdge = 1;
leftEdge = 3;
width = 7.5; % one column: 3+9/16, two column: 7.5
height = 6;
fontSize = 8;
lineWidth = 1;

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];

ax(1) = subplot(2,1,1);
plot(out.t,out.p_h)
hold on
plot(out.t,out.p_ro)
xlabel('Time (s)')
ylabel('Pressure (Pa)')
legend('p_{h}','p_{ro}')

ax(2) = subplot(2,1,2);
plot(out.t,out.p_l)
hold on
xlabel('Time (s)')
ylabel('Pressure (Pa)')
legend('p_{l}')

linkaxes(ax,'x');

sgtitle('Pressures')

%% Flow rates
bottomEdge = 1;
leftEdge = 3;
width = 7.5; % one column: 3+9/16, two column: 7.5
height = 6;
fontSize = 8;
lineWidth = 1;

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];

ax(1) = subplot(3,1,1);
plot(out.t,out.q_h(:))
hold on
plot(out.t,out.q_pm(:))
plot(out.t,out.q_rv(:))
xlabel('Time (s)')
ylabel('Flow rate (m^3/s)')
legend('q_{h}','q_{pm}','q_{rv}')

ax(2) = subplot(3,1,2);
plot(out.t,out.q_rv(:))
hold on
plot(out.t,out.q_ERUfeed(:))
plot(out.t,out.q_feed(:))
plot(out.t,out.q_perm(:))
xlabel('Time (s)')
ylabel('Flow rate (m^3/s)')
legend('q_{rv}','q_{ERUfeed}','q_{feed}','q_{perm}')

ax(3) = subplot(3,1,3);
plot(out.t,out.q_l(:))
hold on
plot(out.t,out.q_pm(:))
plot(out.t,out.q_c(:))
plot(out.t,out.q_ERUfeed(:))
xlabel('Time (s)')
ylabel('Flow rate (m^3/s)')
legend('q_{l}','q_{pm}','q_{c}','q_{ERUfeed}')

linkaxes(ax,'x')

sgtitle('Flow rates')


%% Controller behavior
bottomEdge = 1;
leftEdge = 3;
width = 7.5; % one column: 3+9/16, two column: 7.5
height = 8;
fontSize = 8;
lineWidth = 1;

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];

ax(1) = subplot(4,1,1);
plot(out.t,1e-6*out.p_h)
hold on
plot(out.t,1e-6*out.p_ro)
plot(out.t,1e-6*out.control.p_filt)
xlabel('Time (s)')
ylabel('Pressure (MPa)')
legend('p_{h}','p_{ro}','p_{filt}')

ax(2) = subplot(4,1,2);
hold on
yyaxis left
plot(out.t,60/(2*pi)*out.control.w_pm)
plot(out.t,60/(2*pi)*out.w_pm)
ylabel('Shaft speed (rpm)')
yyaxis right
plot(out.t,1e-3*out.Tgen)
ylabel('Torque (kNm)')
legend('nominal','actual','Generator')

ax(3) = subplot(4,1,3);
hold on
plot(out.t,1e3*out.q_pm)
plot(out.t,1e3*out.q_h)
ylabel('Flow rate (Lpm)')
legend('q_{pm}','q_{h}')

ax(4) = subplot(4,1,4);
hold on
yyaxis left
plot(out.t,out.control.errInt_p_filt)
ylabel('Error integral')
legend('pressure control')
ylabel('Error integral')
xlabel('Time (s)')

linkaxes(ax,'x')

sgtitle('Controller behaviour')

%% WEC-driven pump
bottomEdge = 1;
leftEdge = 3;
width = 7.5; % one column: 3+9/16, two column: 7.5
height = 6;
fontSize = 8;
lineWidth = 1;

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];

ax(1) = subplot(2,1,1);
plot(out.t,out.p_a)
hold on
plot(out.t,out.p_b)
plot(out.t,out.p_h,'k-')
plot(out.t,out.p_l,'k--')
xlabel('Time (s)')
ylabel('Pressure (Pa)')
legend('p_{a}','p_{b}','p_{h}','p_{l}')

ax(2) = subplot(2,1,2);
plot(out.t,out.q_h,'k-')
hold on
plot(out.t,out.q_l,'k--')
plot(out.t,out.q_sv,'b')
plot(out.t,out.q_aPRV+out.q_bPRV,'r')
xlabel('Time (s)')
ylabel('Flow rate (Lpm)')
legend('q_{h}','q_{l}','q_{sv}','q_{aPRV} + q_{bPRV}')

linkaxes(ax,'x');

sgtitle('WEC-driven Pump Dynamics')
