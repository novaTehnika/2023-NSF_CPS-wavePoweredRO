function [dydt, nonState, control] = sys_refPTO(t,y,par)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sys_refPTO.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 6/12/2023
%
% PURPOSE/DESCRIPTION:
% Calculate the state derivatives for a wave energy PTO specified in the 
% NSF-Cyber-Physical Systems project proposal.
%
% FILE DEPENDENCY:
%
% UPDATES:
% 6/12/2023 - created from sys_parallelPTO.m.
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
% Define indicies of state vector
iyp_a = [];
iyp_b = [];

iyp_l = [];
iyp_h = [];
iyp_ro = [];

iyp_filt = [];
iy_errInt_p_filt = [];
iy_errInt_w_pm = [];
iycontrol = [];

iytheta = [];
iytheta_dot = [];
iyrad = [];

stateIndex_refPTO

iWEC = [iytheta iytheta_dot iyrad];


% Calculate control input and the change in controller states (if any)
[dydt_control, control] = controller(t,y,par);

% Calculate non-state variables like coeffs., forces, and flow rates
nonState = nonStateVars(t,y,par);

% Calculate the hydrodynamics WEC state derivatives and output nonstate
% variables like forces and wave elevation
[dydt_WEC, nonState.torqueWEC, nonState.waveElev] = ...
                    flapModel(t,y(iWEC),nonState.T_pto,par);

% State derivatives
dydt = zeros(11 + par.WEC.ny_rad,1);

dydt(iyp_a) = 1/nonState.C_a*(par.D_WEC*y(iytheta_dot) - nonState.q_sv ...
                + nonState.q_ain - nonState.q_aout ...
                - nonState.q_aPRV + nonState.q_bPRV);
dydt(iyp_b) = 1/nonState.C_b*(-par.D_WEC*y(iytheta_dot) + nonState.q_sv ...
                + nonState.q_bin - nonState.q_bout ...
                + nonState.q_aPRV - nonState.q_bPRV);

dydt(iyp_l) = 1/nonState.C_l*(nonState.q_c ...
                - nonState.q_ain - nonState.q_bin ...
                + nonState.q_pm - nonState.q_ERUfeed ...
                - nonState.q_lPRV + nonState.q_hPRV + nonState.q_roPRV);
dydt(iyp_h) = 1/nonState.C_h*(nonState.q_aout - nonState.q_bout ...
                - nonState.q_pm - nonState.q_rv - nonState.q_hPRV);
dydt(iyp_ro) = 1/nonState.C_ro*(nonState.q_rov - nonState.q_feed ...
                + nonState.q_ERUfeed - nonState.q_roPRV);

dydt(iycontrol) = dydt_control;

dydt(iytheta) = dydt_WEC(1); % angular velocity
dydt(iytheta_dot) = dydt_WEC(2); % angular acceleration
dydt(iyrad) = dydt_WEC(3:end); % radiation damping states for WEC model

%% %%%%%%%%%%%%   FUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [dydt_control, control] = controller(t,y,par)
        %% PI control of P_Hout using w_pm
         % Error
        err_p_filt = y(iyp_filt) - par.control.p_hout_nom;
         % Feedforward
        w_ff = 0*par.control.w_pm_ctrl.min;
         % Control signal
        w_pm_nom = w_ff ...
            + (par.control.w_pm_ctrl.kp*err_p_filt ...
            + par.control.w_pm_ctrl.ki*y(iy_errInt_p_filt));
        control.w_pm_nom = (w_pm_nom <= par.control.w_pm_ctrl.max ...
                          & w_pm_nom >= par.control.w_pm_ctrl.min) ...
                          * w_pm_nom ...
                          + (w_pm_nom > par.control.w_pm_ctrl.max) ...
                          * par.control.w_pm_ctrl.max ...
                          + (w_pm_nom < par.control.w_pm_ctrl.min) ...
                          * par.control.w_pm_ctrl.min;
                      
        %% deriviatives for filtered signal and error integrals (w/
         % anti-wind-up)
        dydt_control = [    % filtered signal
                        (y(iyp_hout) - y(iyp_filt))...
                        /par.control.tau_pfilt;

                            % error integral for pressure control
                        (y(iy_errInt_p_filt) < 1e12 ...
                        & y(iy_errInt_p_filt) > -1e12) ...
                        *err_p_filt];
                        
    end

    function nonState = nonStateVars(t,y,par)
        % Accumulator capacitance
        nonState.C_l = capAccum(y(iyp_l),par.pc_l,par.Vc_l,par.f,par);
        nonState.C_lout = capAccum(y(iyp_lout),par.pc_lout,par.Vc_lout,par.f,par) ...
                        + lineCap(y(iyp_lout),1,par);
        nonState.C_hin = capAccum(y(iyp_hin),par.pc_hin,par.Vc_hin,par.f,par) ...
                        + lineCap(y(iyp_hin),2,par);
        nonState.C_hout = capAccum(y(iyp_hout),par.pc_hout,par.Vc_hout,par.f,par) ...
                        + lineCap(y(iyp_hout),2,par);

        % WEC-driven pump
         % pumping chamber capacitance
        nonState.C_a = capWEC(y(iytheta),y(iyp_a),par) ...
                        + deadVCap(y(iyp_a),par.VwecDead,par);
        nonState.C_b = capWEC(-y(iytheta),y(iyp_b),par) ...
                        + deadVCap(y(iyp_b),par.VwecDead,par);
         % Switching valve flow
        dp = y(iyp_a) - y(iyp_b);
        nonState.q_sv = areaFracPWM(t,par.duty_sv,par.T_sv) ...
                        *((abs(dp) <= par.dp_svlin)*par.kv_svlin*dp ...
                        + (abs(dp) > par.dp_svlin)*par.kv_sv*sqrt(dp));
         % check valves
        nonState.q_ain = flowCV(y(iyp_l) - y(iyp_a), ...
                         par.kvWECin,par.pc_WECin,par.dp_WECin);
        nonState.q_aout = flowCV(y(iyp_a) - y(iyp_h), ...
                         par.kvWECout,par.pc_WECout,par.dp_WECout);
        nonState.q_bin = flowCV(y(iyp_l) - y(iyp_b), ...
                         par.kvWECin,par.pc_WECin,par.dp_WECin);
        nonState.q_bout = flowCV(y(iyp_b) - y(iyp_h), ...
                         par.kvWECout,par.pc_WECout,par.dp_WECout);
        


        delta_p_wp = y(iyp_a)-y(iyp_b);
        WECpumpPumping = y(iytheta_dot)*delta_p_wp < 0;
        WECpumpMotoring = y(iytheta_dot)*delta_p_wp >= 0;
        nonState.q_w = (WECpumpPumping*par.eta_v_WEC ...
                    + WECpumpMotoring/par.eta_v_WEC)...
                    * par.D_WEC*abs(y(iytheta_dot));

        nonState.T_pto = -sign(y(iytheta_dot))...
                    * (WECpumpPumping/par.eta_m_WEC ...
                    + WECpumpMotoring*par.eta_m_WEC)...
                    * (abs(y(iytheta_dot)) > 5e-3)...
                    * par.D_WEC*delta_p_wp;

        % nonState.pLinePPfricHP = pLsoln.PPfric;

        % house power pump/motor
        delta_p_pm = y(iyp_lin) - y(iyp_hout);
        nonState.pmPumping = y(iyw_pm)*(delta_p_pm) >= 0;
        nonState.pmMotoring = y(iyw_pm)*(delta_p_pm) < 0;

        volLoss_pm = par.APP.C_s*(delta_p_pm/(par.mu*abs(y(iyw_pm)))) ...
                           + (delta_p_pm/par.beta)*(par.APP.V_r + 1);
        nonState.eta_v_pm = nonState.pmPumping*(1 - volLoss_pm) ...
                          + nonState.pmMotoring/(1 + volLoss_pm);
        nonState.q_pm = (nonState.pmPumping*nonState.eta_v_pm ...
                        + nonState.pmMotoring/nonState.eta_v_pm) ...
                        *par.D_pm*y(iyw_pm);

        mechLoss_pm = par.APP.C_v*par.mu*abs(y(iyw_pm))/delta_p_pm + par.APP.C_f;
        nonState.eta_m_pm = nonState.pmPumping/(1 + mechLoss_pm) ...
                          + nonState.pmMotoring*(1 - mechLoss_pm);
        nonState.Tpm = (nonState.pmPumping/nonState.eta_m_pm ...
                        + nonState.pmMotoring*nonState.eta_m_pm)...
                        *par.D_pm*delta_p_pm;

        % Reverse osmosis module
        nonState.q_perm = par.Sro*par.Aperm*(y(iyp_hout) - par.p_perm - par.p_osm);
        nonState.q_feed = 1/par.Y*nonState.q_perm;

        % ERU
        nonState.q_brine = nonState.q_feed - nonState.q_perm;
        nonState.q_ERUfeed = (par.eta_ERUv)^2*nonState.q_brine;

        % Charge Pump
        dP_SO = (y(iyp_lin) - par.p_o) - par.cn*par.w_c^2; % difference between shut-off pressure and current pressure differential
        nonState.q_c = (dP_SO < 0)* sqrt(dP_SO/par.cq);

        % Pressure relief valves
        nonState.q_aPRV = prv(y(iyp_a),par.aPRV.p_crack,par.aPRV.C);
        nonState.q_bPRV = prv(y(iyp_b),par.bPRV.p_crack,par.bPRV.C);
        nonState.q_lPRV = prv(y(iyp_l),par.lPRV.p_crack,par.lPRV.C);
        nonState.q_hPRV = prv(y(iyp_h),par.hPRV.p_crack,par.hPRV.C);
        nonState.q_roPRV = prv(y(iyp_ro),par.roPRV.p_crack,par.roPRV.C);

    end

end
