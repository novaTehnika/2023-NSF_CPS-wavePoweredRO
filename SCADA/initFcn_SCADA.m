% loadModelParams.m script m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 5/21/2025
%
% PURPOSE/DESCRIPTION:
% This script loads all required parameters as Simulink.Parameter objects
% for use in the Simulink Model "SCADA.slx".
%
% FILE DEPENDENCY:
%
% UPDATES:
% 5/21/2025 - Created.
%
% Copyright (C) 2025  Jeremy W. Simmons II
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

addpath('WEC model') 
addpath(['WEC model' filesep 'WECdata'])  
addpath('Sea States')

%%
par.wave.Hs = 2;
par.wave.Tp = 10;

par.WEC.nw = 1000; % num. of frequency components for harmonic superposition 
par.wave.rngSeedPhase = 3; % seed for the random number generator

par = parameters_WECmodel(par,...
    'nemohResults_vantHoff2009_20180802.mat','vantHoffTFCoeff.mat');

par = parameters_HIL(par);

%% System parameters for simulink
%% General
rho = 1023; % [kg/m3] density of seawater
g = 9.81;   % [m/s^2] gravitational acceleration at sea level

%% WEC    
WEC_w = Simulink.Parameter(par.WEC.w);
WEC_w.CoderInfo.StorageClass = 'ExportedGlobal';

WEC_dw = Simulink.Parameter(par.WEC.dw);
WEC_dw.CoderInfo.StorageClass = 'ExportedGlobal';

WEC_phi_e = Simulink.Parameter(par.WEC.phi_e);
WEC_phi_e.CoderInfo.StorageClass = 'ExportedGlobal';

WEC_F_amp = Simulink.Parameter(par.WEC.F_amp);
WEC_F_amp.CoderInfo.StorageClass = 'ExportedGlobal'; 

WEC_A_rad = Simulink.Parameter(par.WEC.A_rad);
WEC_A_rad.CoderInfo.StorageClass = 'ExportedGlobal';

WEC_B_rad = Simulink.Parameter(par.WEC.B_rad);
WEC_B_rad.CoderInfo.StorageClass = 'ExportedGlobal';

WEC_L_flap = Simulink.Parameter(par.WEC.L_flap);
WEC_L_flap.CoderInfo.StorageClass = 'ExportedGlobal';

WEC_T = Simulink.Parameter(par.WEC.T);
WEC_T.CoderInfo.StorageClass = 'ExportedGlobal';

WEC_h = Simulink.Parameter(par.WEC.h);
WEC_h.CoderInfo.StorageClass = 'ExportedGlobal';

WEC_W = Simulink.Parameter(par.WEC.W);
WEC_W.CoderInfo.StorageClass = 'ExportedGlobal';

WEC_L_cm = Simulink.Parameter(par.WEC.L_cm);
WEC_L_cm.CoderInfo.StorageClass = 'ExportedGlobal';

WEC_m = Simulink.Parameter(par.WEC.m);
WEC_m.CoderInfo.StorageClass = 'ExportedGlobal';

WEC_I = Simulink.Parameter(par.WEC.I);
WEC_I.CoderInfo.StorageClass = 'ExportedGlobal';

WEC_I_inf = Simulink.Parameter(par.WEC.I_inf);
WEC_I_inf.CoderInfo.StorageClass = 'ExportedGlobal';

WEC_rho = Simulink.Parameter(par.WEC.rho);
WEC_rho.CoderInfo.StorageClass = 'ExportedGlobal';

WEC_g = Simulink.Parameter(par.WEC.g);
WEC_g.CoderInfo.StorageClass = 'ExportedGlobal';

WEC_ICs = Simulink.Parameter(zeros(1,par.WEC.ny_rad));
WEC_ICs.CoderInfo.StorageClass = 'ExportedGlobal';

%% Wave
wave_phi = Simulink.Parameter(par.wave.phi);
wave_phi.CoderInfo.StorageClass = 'ExportedGlobal';

wave_S_w = Simulink.Parameter(par.wave.S_w);
wave_S_w.CoderInfo.StorageClass = 'ExportedGlobal';

%% HIL system
HIL_scale_sys = Simulink.Parameter(par.scale_sys);
HIL_scale_sys.CoderInfo.StorageClass = 'ExportedGlobal';

 % WEC-driven Pump
HIL_WECpump_ang2lin = Simulink.Parameter(par.WECpump.ang2lin);
HIL_ang2lin.CoderInfo.StorageClass = 'ExportedGlobal';

HIL_WECpump_A_cap = Simulink.Parameter(par.WECpump.A_cap);
HIL_WECpump_A_cap.CoderInfo.StorageClass = 'ExportedGlobal';

HIL_WECpump_A_rod = Simulink.Parameter(par.WECpump.A_rod);
HIL_WECpump_A_rod.CoderInfo.StorageClass = 'ExportedGlobal';

HIL_WECpump_positiveLimit = Simulink.Parameter(par.WECpump.positiveLimit);
HIL_WECpump_positiveLimit.CoderInfo.StorageClass = 'ExportedGlobal';

HIL_WECpump_negativeLimit = Simulink.Parameter(par.WECpump.negativeLimit);
HIL_WECpump_negativeLimit.CoderInfo.StorageClass = 'ExportedGlobal';

HIL_WECpump_positiveCtrlLimit = Simulink.Parameter(par.WECpump.positiveCtrlLimit);
HIL_WECpump_positiveCtrlLimit.CoderInfo.StorageClass = 'ExportedGlobal';

HIL_WECpump_negativeCtrlLimit = Simulink.Parameter(par.WECpump.negativeCtrlLimit);
HIL_WECpump_negativeCtrlLimit.CoderInfo.StorageClass = 'ExportedGlobal';

 % Actuator Control
  % Actuator/Actuator Control
HIL_act_A_cap = Simulink.Parameter(par.actuator.A_cap);
HIL_act_A_cap.CoderInfo.StorageClass = 'ExportedGlobal';

HIL_act_A_rod = Simulink.Parameter(par.actuator.A_rod);
HIL_act_A_rod.CoderInfo.StorageClass = 'ExportedGlobal';

  % Danfoss H1T/Actuator Control
HIL_H1T_D = Simulink.Parameter(par.H1T.D);
HIL_H1T_D.CoderInfo.StorageClass = 'ExportedGlobal';

HIL_H1T_P1posCutIn = Simulink.Parameter(par.H1T.P1.posCutIn);
HIL_H1T_P1posCutIn.CoderInfo.StorageClass = 'ExportedGlobal';

HIL_H1T_P1negCutIn = Simulink.Parameter(par.H1T.P1.negCutIn);
HIL_H1T_P1negCutIn.CoderInfo.StorageClass = 'ExportedGlobal';

HIL_H1T_P1fBias = Simulink.Parameter(par.H1T.P1.fBias);
HIL_H1T_P1fBias.CoderInfo.StorageClass = 'ExportedGlobal';

HIL_H1T_P2posCutIn = Simulink.Parameter(par.H1T.P2.posCutIn);
HIL_H1T_P2posCutIn.CoderInfo.StorageClass = 'ExportedGlobal';

HIL_H1T_P2negCutIn = Simulink.Parameter(par.H1T.P2.negCutIn);
HIL_H1T_P2negCutIn.CoderInfo.StorageClass = 'ExportedGlobal';

HIL_H1T_P2fBias = Simulink.Parameter(par.H1T.P2.fBias);
HIL_H1T_P2fBias.CoderInfo.StorageClass = 'ExportedGlobal';

HIL_H1T_kp = Simulink.Parameter(par.H1T.kp);
HIL_H1T_kp.CoderInfo.StorageClass = 'ExportedGlobal';

HIL_H1T_ki = Simulink.Parameter(par.H1T.ki);
HIL_H1T_ki.CoderInfo.StorageClass = 'ExportedGlobal';

HIL_H1T_kd = Simulink.Parameter(par.H1T.kd);
HIL_H1T_kd.CoderInfo.StorageClass = 'ExportedGlobal';


 % Charge Pump Control
HIL_charge_nMax = Simulink.Parameter(par.charge.nMax);
HIL_charge_nMax.CoderInfo.StorageClass = 'ExportedGlobal';

HIL_charge_nMin = Simulink.Parameter(par.charge.nMin);
HIL_charge_nMin.CoderInfo.StorageClass = 'ExportedGlobal';

HIL_charge_cn = Simulink.Parameter(par.charge.cn);
HIL_charge_cn.CoderInfo.StorageClass = 'ExportedGlobal';

HIL_charge_kp = Simulink.Parameter(par.charge.kp);
HIL_charge_kp.CoderInfo.StorageClass = 'ExportedGlobal';

HIL_charge_ki = Simulink.Parameter(par.charge.ki);
HIL_charge_ki.CoderInfo.StorageClass = 'ExportedGlobal';

HIL_charge_kd = Simulink.Parameter(par.charge.kd);
HIL_charge_kd.CoderInfo.StorageClass = 'ExportedGlobal';


 % Motor/Gen. Control
HIL_gen_nMax = Simulink.Parameter(par.gen.nMax);
HIL_gen_nMax.CoderInfo.StorageClass = 'ExportedGlobal';

HIL_gen_nMin = Simulink.Parameter(par.gen.nMin);
HIL_gen_nMin.CoderInfo.StorageClass = 'ExportedGlobal';

HIL_gen_kp = Simulink.Parameter(par.gen.kp);
HIL_gen_kp.CoderInfo.StorageClass = 'ExportedGlobal';

HIL_gen_ki = Simulink.Parameter(par.gen.ki);
HIL_gen_ki.CoderInfo.StorageClass = 'ExportedGlobal';

HIL_gen_kd = Simulink.Parameter(par.gen.kd);
HIL_gen_kd.CoderInfo.StorageClass = 'ExportedGlobal';


 % Motor/Gen. Control
HIL_ERU_nMax = Simulink.Parameter(par.ERU.nMax);
HIL_ERU_nMax.CoderInfo.StorageClass = 'ExportedGlobal';

HIL_ERU_nMin = Simulink.Parameter(par.ERU.nMin);
HIL_ERU_nMin.CoderInfo.StorageClass = 'ExportedGlobal';

HIL_ERU_kp = Simulink.Parameter(par.ERU.kp);
HIL_ERU_kp.CoderInfo.StorageClass = 'ExportedGlobal';

HIL_ERU_ki = Simulink.Parameter(par.ERU.ki);
HIL_ERU_ki.CoderInfo.StorageClass = 'ExportedGlobal';

HIL_ERU_kd = Simulink.Parameter(par.ERU.kd);
HIL_ERU_kd.CoderInfo.StorageClass = 'ExportedGlobal';

%% runtime tunable parameters
% motionOnOff = Simulink.Parameter(0);
% motionOnOff.CoderInfo.StorageClass = 'ExportedGlobal';

    function S_w = PiersonSpec(w,par)
        % Based on Falnes (2002) "Ocean Waves and Oscillating Systems:..."
        S_w = 10*pi^5*par.wave.Hs^2/par.wave.Tp^4./w.^5 ... 
            .*exp(-20*pi^4/par.wave.Tp^4./w.^4)/(2*pi);
    end