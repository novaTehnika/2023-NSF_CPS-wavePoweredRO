function par = parameters_refPTO(par,filenameCoeff,filenameRadSS)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters_refPTO.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 6/12/2023
%
% PURPOSE/DESCRIPTION:
% This function loads default parameters into the "par" structure. This
% includes calling a simular function to load parameters for the WEC model.
% File names for the hydrodynamic data for the WEC are passed through that
% function. Not all parameters can be modified after this function runs
% without issue. Changing parameters in any design study script should be
% check for affecting the other parameters.
%
% This function is for a use with sys_refPTO.m.
%
% FILE DEPENDENCY:
% parameters_WECmodel.m
%
% UPDATES:
% 6/12/2023 - created from parameters_parallelPTO.m.
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

    % fluid and entrianed gas properties
    par.rho = 1023; % [kg/m3] density of air
    par.mu = 9.4e-4; % [Pa-s]  Dynamic (absolute) viscosity
    par.beta = 2.2e9; % [Pa]  Bulk Modulus of air free fluid
    par.p_vap = 0.037e5; % [Pa] vapour pressure of seawater
    par.R = 0.0001; % [-] fraction  Baseline fraction of air by volume entrained in the hydraulic fluid at atm
    par.p_o = 101.3e3; % [Pa]  Atmospheric pressure (reference)
    par.gamma = 1.4; % [-]ratio of specific heats for air
    
    % WEC parameters
    par = parameters_WECmodel(par,filenameCoeff,filenameRadSS);

    % WEC-pump
     % pumping chamber
    par.D_WEC = 0.3;         % [m^3/rad] flap pump displacement
    par.eta_v_WEC = 1;
    par.eta_m_WEC = 0.9;                % Flap pump mechanical efficiency

     % switching valve
    par.duty_sv = 0;
    par.T_sv = 0.5; % [s] switching period (1/switching freq)
    par.tr_sv = 0.05; % transition ratio (frac. of cyc. transitioning 1way)

    dp_rated = 1e5; % [Pa] 
    q_rated = (100)*60/1e3; % [(lpm) -> m^3/s]
    par.kv_sv = q_rated/dp_rated;
    par.dp_svlin = 1e4; % [Pa] linearization margin
    par.kv_svlin = par.kv_sv*sqrt(par.dp_svlin)/par.dp_svlin;

     % check valves
      % inlet check valves (low-pressure)
    dp_rated = 1e5; % [Pa] 
    q_rated = (100)*60/1e3; % [(lpm) -> m^3/s]
    par.kvWECin = q_rated/dp_rated;
    par.pc_WECin = 1e5; % [Pa] Cracking pressure
    par.dp_WECin = 1e5; % [Pa] margin between cracking pressure and fully open condition
    
     % outlet check valves (high-pressure)
    dp_rated = 2e5; % [Pa] 
    q_rated = (100)*60/1e3; % [(lpm) -> m^3/s]
    par.kvWECout = q_rated/dp_rated;
    par.pc_WECout = 1e5; % [Pa] Cracking pressure
    par.dp_WECout = 1e5; % [Pa] margin between cracking pressure and fully open condition

    % RO module parameters
    par.p_perm = par.p_o;
    par.p_osm = 2.7e6;
    par.Sro = 3000; % [m^3]
    par.Aperm = 2.57e-12; % [m^3/(N-s)] permeabiity coefficient (Yu and Jenne,2018)
    par.Y = 0.25;

    % ERU
    par.eta_ERUv = 0.95;
    par.eta_ERUm = 0.95;
    
    % power control unit
      % pump/motor
    par.D_pm = (6500)*1e-6/(2*pi); % [(cc/rev) -> m^3/rad]  Motor displacement
    par.w_pm_max = (3600)/60*2*pi; % [(rpm) -> rad/s] maximum speed of motor
    par.w_pm_min = (1)/60*2*pi; % [(rpm) -> rad/s] minimum speed of motor

      % efficiency model coeffs. (McCandlish and Dory model)
       % Axial piston pump (data from Danfoss APP 43/1700)
    par.APP.C_s = 3.0554e-10;
    par.APP.V_r = 1.103;
    par.APP.C_v = 7.1755e5;
    par.APP.C_f = 0.0259;
    
      % generator
    par.eta_g = 0.9;  % efficiency of electric generator
       % rotor inertia estimated by NEMA MG-1 (14.46)
    % n_poles = 3; % [qty.] number of poles of induction motor
    % genPowerRating = 100; % [hp] power rating of induction motor
    % par.Jpm = ( 0.02*2^n_poles*genPowerRating^(1.35-.05*n_poles/2) ) ...
    %             * 0.0421401101; % [lb ft^2 --> kg m^2] rotor inertia
    
    % Accumulators
    par.f = 1e-2; % fraction of dead volume of working fluid compared to 
                  % charge volume
     % LPA
    par.Vc_l = (3000)*1e-3; % [(L) -> m^3] gas volume at charge pressure
    par.pc_l = 0.2e6; % [Pa] charge pressure
     % HPA at outlet of WEC-driven pump
    par.Vc_h = (6000)*1e-3; % [(L) -> m^3] gas volume at charge pressure
    par.pc_h = 4e6; % [Pa] charge pressure
     % HPA at inlet to RO module
    par.Vc_ro = (6000)*1e-3; % [(L) -> m^3] gas volume at charge pressure
    par.pc_ro = 4e6; % [Pa] charge pressure

    % Contoller Parameters
    par.control.p_h_nom = 6e6; % [Pa]
    par.control.p_RO_nom = par.control.p_h_nom; % [Pa] (not actually used as control ref.)
    par.control.p_l_nom = 0.5e6; % [Pa] (not actually used as control ref.)
    par.control.p_RO_max = 8.2e6; % [Pa]
    par.control.p_RO_min = max(3e6,par.pc_hout); % [Pa]

     % Signal filtering
    par.control.tau_pfilt = 0.01; % [s] time constant for LPF for pressure signal
    
     % PI control of w_pm using T_gen
    par.control.wpm_ctrl.kp = 1e3;
    par.control.wpm_ctrl.ki = 2e5;
    
    % Charging system (Intake & Boost pump)
    par.cn = 7;
    par.cq = -1e6;
    par.w_c = (3600)*2*pi/60; % [(rpm) -> rad/s]
    par.eta_c = 0.7;  % pumping efficiency of pressure boost pump
    par.eta_m = 0.9;  % efficiency of charge pump motor
    par.p_tank = .65e6;

    % Pressure relief valves
         % WEC-driven pump chamber 'a'
    maxPressure = 10e6; % [Pa]
    margin = 5e4; % [Pa]
    maxFlow = 100e-3; % [m^3/s]
    par.aPRV.p_crack = maxPressure - margin;
    par.aPRV.C = (maxPressure^(3/2) ...
                 - (maxPressure-margin)*maxPressure^(1/2))/maxFlow;

         % WEC-driven pump chamber 'b'
    maxPressure = 10e6; % [Pa]
    margin = 5e4; % [Pa]
    maxFlow = 100e-3; % [m^3/s]
    par.bPRV.p_crack = maxPressure - margin;
    par.bPRV.C = (maxPressure^(3/2) ...
                 - (maxPressure-margin)*maxPressure^(1/2))/maxFlow;

     % low-pressure inlet to WEC-driven pump/outlet of charge pump
    maxPressure = 10e6; % [Pa]
    margin = 5e4; % [Pa]
    maxFlow = 100e-3; % [m^3/s]
    par.lPRV.p_crack = maxPressure - margin;
    par.lPRV.C = (maxPressure^(3/2) ...
                 - (maxPressure-margin)*maxPressure^(1/2))/maxFlow;

     % high-pressure outlet of WEC-driven pump
    maxPressure = 300e6; % [Pa]
    margin = 5e4; % [Pa]
    maxFlow = (100)*1e-3; % [(L/s) -> m^3/s]
    par.hPRV.p_crack = maxPressure - margin;
    par.hPRV.C = (maxPressure^(3/2) ...
                 - (maxPressure-margin)*maxPressure^(1/2))/maxFlow;
    
     % high-pressure inlet to RO module
    maxPressure = 8.3e6; % [Pa]
    margin = 5e4; % [Pa]
    maxFlow = 100e-3; % [m^3/s]
    par.roPRV.p_crack = maxPressure - margin;
    par.roPRV.C = (maxPressure^(3/2) ...
                 - (maxPressure-margin)*maxPressure^(1/2))/maxFlow;


end