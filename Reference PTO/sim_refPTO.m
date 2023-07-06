function out = sim_refPTO(y0,par)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sim_refPTO.m function m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 6/12/2023
%
% PURPOSE/DESCRIPTION:
% This script executes the set-up, solution, and basic post-processing for
% the refPTO model.
%
% FILE DEPENDENCY:
% ../Reference PTO/
%   sys_refPTO.m
%   stateIndex_refPTO.m
% ../WEC model/
%   flapModel.m
%   hydroStaticTorque.m
% ../Solvers/
%   deltaE_NI.m
%   deltaV_NI.m
%   ode1.m
% ../Components/
%   areaFracPWM.m
%   capAccum.m
%   deadVCap.m
%   flowCV.m
%   flowPRV.m
%
% UPDATES:
% 6/12/2023 - Created from sim_parallelPTO.m.
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
iycontrol = [];

iytheta = [];
iytheta_dot = [];
iyrad = [];

stateIndex_refPTO % load state indices

%% %%%%%%%%%%%%   SOLUTION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Solve for states of the dynamic system
 % Set-up solver
    tspan = [par.tstart-par.Tramp par.tend];   % time interval
    
 % Solver options
%     options = odeset('RelTol',par.odeSolverRelTol,...
%                      'AbsTol',par.odeSolverAbsTol,...
%                      'MaxOrder',3);
%     options.MaxStep = par.MaxStep;
    dt = par.MaxStep;
    downSampleRate = floor(par.downSampledStepSize/dt);

 % Run solver
    ticODE = tic;
    [t, y] = ode1(@(t,y) sys(t,y,par)', ...
                                tspan(1),dt,tspan(2),y0,downSampleRate);
    toc(ticODE)

%% %%%%%%%%%%%%   POST-PROCESS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameters
    out.par = par;
    
    % Select desired time indices
    itVec = find(t >= par.tstart);
    
    % Extract system states from simulation results post ramp
    out.t = t(itVec);
    out.y = y(itVec,:);
    
    out.p_a = y(itVec,iyp_a);
    out.p_b = y(itVec,iyp_b);
    
    out.p_l = y(itVec,iyp_l);
    out.p_h = y(itVec,iyp_h);
    out.p_ro = y(itVec,iyp_ro);
    
    out.control.p_filt = y(itVec,iyp_filt);
    out.control.errInt_p_filt = y(itVec,iy_errInt_p_filt);
    
    out.theta = y(itVec,iytheta); % [rad] position of the WEC
    out.theta_dot = y(itVec,iytheta_dot); % [rad/s] angular velocity of the WEC
            
    % Post-process non-state variables and state derivatives
    syspost = @(t,y,par) sysPost(t,y,par);

    nt_ramp = itVec(1)-1;
    startParPool
    parfor it = 1:length(itVec)
        
        [dydt(it,:), nonState(it), control(it)] = ...
                            syspost(t(it+nt_ramp),y(it+nt_ramp,:),par);
        
        % Move WEC torque results up a level in nonState stucture so that
        % they can be used like arrays in assiging to the output structure
        temp(it).T_hydroStatic = nonState(it).torqueWEC.hydroStatic;
        temp(it).T_wave = nonState(it).torqueWEC.wave;
        temp(it).T_rad = nonState(it).torqueWEC.radiation;    
    end
    
    % State derivatives
    out.dydt = dydt;
    
     % System input: wave elevation
    out.waveElev = [nonState(:).waveElev]';
    
     % forces on WEC
    out.T_pto = [nonState(:).T_pto]';
    out.T_hydroStatic = [temp(:).T_hydroStatic]';
    out.T_wave = [temp(:).T_wave]';
    out.T_rad = [temp(:).T_rad]';
          
     % Control signals
    out.control.w_pm = [control(:).w_pm]';
    out.control.kv_rv = [control(:).kv_rv]';
    out.control.kv_ideal = [control(:).kv_ideal]';

     % torque on pump/motor and generator shafts
    out.Tpm = [nonState(:).Tpm]';
    out.Tgen = out.Tpm;
    
     % WEC-driven pump flow
    out.q_ain = [nonState(:).q_ain]';
    out.q_aout = [nonState(:).q_aout]';
    out.q_bin = [nonState(:).q_bin]';
    out.q_bout = [nonState(:).q_bout]';
    out.q_sv = [nonState(:).q_sv]';

    out.q_h = [nonState(:).q_aout]' + [nonState(:).q_bout]';
    out.q_l = [nonState(:).q_ain]' + [nonState(:).q_bin]';

     % pump/motor flow
    out.w_pm = out.control.w_pm;
    out.q_pm = [nonState(:).q_pm]';
    
     % RO performance
    out.q_perm = [nonState(:).q_perm]';
    out.q_brine = [nonState(:).q_brine]';
    out.q_feed = out.q_perm + out.q_brine;

     % RO inlet valve
    out.q_rv = [nonState(:).q_rv]';

     % ERU flow rate
    out.q_ERUfeed = [nonState(:).q_ERUfeed]';

     % Charge Pump
    out.q_c = [nonState(:).q_c]';

     % Pressure relief valves at system pressure nodes
    out.q_lPRV = [nonState(:).q_lPRV]';
    out.q_hPRV = [nonState(:).q_hPRV]';
    out.q_roPRV = [nonState(:).q_roPRV]';

%% Post-process analysis

    %% Energy analysis
    % WEC-driven pump
    out.power.P_WEC = -out.T_pto.*out.theta_dot;
    out.power.P_wp = out.p_h.*out.q_h - out.p_l.*out.q_l;
    out.power.P_wpLoss = out.power.P_WEC - out.power.P_wp;
    
     % switching valve
    out.power.P_sv = out.q_sv.*(out.p_a-out.p_b);
    
     % check valve rectifier
    out.power.P_ain = out.q_ain.*(out.p_l-out.p_a);
    out.power.P_aout = out.q_aout.*(out.p_a-out.p_h);
    out.power.P_bin = out.q_bin.*(out.p_l-out.p_b);
    out.power.P_bout = out.q_bout.*(out.p_b-out.p_h);

    % Pump/motor and generator
    out.power.P_pmLoss = (out.p_l - out.p_h).*out.q_pm ...
                       - out.Tpm.*out.w_pm;
    out.power.P_gen = -((-out.Tgen.*out.w_pm > 0).*par.eta_g ...
                    + (-out.Tgen.*out.w_pm < 0)./par.eta_g) ...
                    .*out.Tgen.*out.w_pm;
    out.power.P_genLoss = -out.Tgen.*out.w_pm - out.power.P_gen;

    % RO ripple control valve
    out.power.P_rv = out.q_rv.*(out.p_h-out.p_ro);

    % Charge pump
    out.power.P_cElec = 1/(par.eta_c*par.eta_m)*out.q_c.*(out.p_l-out.par.p_o);
    out.power.P_cLoss = out.power.P_cElec - out.q_c.*(out.p_l-out.par.p_o);

    % ERU
    dp_ERUfeed = par.ERUconfig.outlet*out.p_ro ...
                + (~par.ERUconfig.outlet)*out.p_h ...
                - out.p_l;
    P_ERUfeed = out.q_ERUfeed.*dp_ERUfeed;
    P_ERUbrine = out.q_brine.*(out.p_ro - out.par.p_o);
    PbalERU = par.ERUconfig.present ...
                *(1/(par.eta_ERUv*par.eta_ERUm)*P_ERUfeed ...
                - par.eta_ERUv*par.eta_ERUm*P_ERUbrine);
    out.power.P_ERULoss = PbalERU + P_ERUbrine - P_ERUfeed;
    out.power.P_ERUelec = ((PbalERU > 0)./par.eta_m ...
                        + (PbalERU < 0).*par.eta_m) ...
                        .*PbalERU;
    out.power.P_ERUelecLoss = out.power.P_ERUelec - PbalERU;

    % Pressure relief valves
    out.power.P_lPRV = out.q_lPRV.*(out.p_l-out.par.p_o);
    out.power.P_hPRV = out.q_hPRV.*(out.p_h-out.p_l);
    out.power.P_roPRV = out.q_roPRV.*(out.p_ro-out.p_l);

    % Electrical Energy Storage
    out.power.P_battery = out.power.P_gen ...
                        - out.power.P_cElec - out.power.P_ERUelec;
    out.power.deltaE_battery = trapz(out.t,out.power.P_battery);

    %% Energy Balance
    % Change in available potential energy in WEC-driven pump chambers
    dp = 1e2;
    cap = @(p,V) deadVCap(p,V,par);
    Va = @(theta) par.V_wecDead + par.D_WEC*(par.theta_max - theta);
    Vb = @(theta) par.V_wecDead + par.D_WEC*(par.theta_max - theta);

     % Chamber 'a'
    cap1 = @(p) cap(p,Va(out.theta(1)));
    capend = @(p) cap(p,Va(out.theta(end)));
    deltaE_a = deltaE_NI(out.par.p_o,out.p_a(end),capend,dp) ...
                - deltaE_NI(out.par.p_o,out.p_a(1),cap1,dp);

     % Chamber 'b'
    cap1 = @(p) cap(p,Vb(out.theta(1)));
    capend = @(p) cap(p,Vb(out.theta(end)));
    deltaE_b = deltaE_NI(out.par.p_o,out.p_a(end),capend,dp) ...
                - deltaE_NI(out.par.p_o,out.p_a(1),cap1,dp);

    % Change in available potential energy in accumulators
    dp = 1e2;
     % Low-pressure
    cap = @(p) capAccum(p,par.pc_l,par.Vc_l,par.f,par);
    deltaE_l = deltaE_NI(out.p_l(1),out.p_l(end),cap,dp);

     % High-pressure outlet of WEC-driven pump
    cap = @(p) capAccum(p,par.pc_h,par.Vc_h,par.f,par);
    deltaE_h = deltaE_NI(out.p_h(1),out.p_h(end),cap,dp);

     % High-pressure inlet of RO module
    cap = @(p) capAccum(p,par.pc_ro,par.Vc_ro,par.f,par);
    deltaE_ro = deltaE_NI(out.p_ro(1),out.p_ro(end),cap,dp);


    % Total change in stored energy in the system
    deltaE_sys = deltaE_a + deltaE_b ...
        + deltaE_l + deltaE_h + deltaE_ro ...
        + out.power.deltaE_battery;

    % Power flow at boundaries
    P_in = out.power.P_WEC;
    P_out = out.q_perm.*(out.p_ro - out.par.p_o);
    P_loss = out.power.P_wpLoss + out.power.P_sv...
            + out.power.P_ain + out.power.P_aout ...
            + out.power.P_bin + out.power.P_bout ...
            + out.power.P_pmLoss + out.power.P_genLoss ...
            + out.power.P_rv ...
            + out.power.P_cLoss ...
            + out.power.P_ERULoss + out.power.P_ERUelecLoss...
            + out.power.P_lPRV + out.power.P_hPRV + out.power.P_roPRV;
    P_bnds = P_in - P_out - P_loss;

    % Total balance of energy: energy added minus change in energy stored
    out.Ebal = trapz(out.t,P_bnds) - deltaE_sys;
    out.Ebal_error = out.Ebal/trapz(out.t,P_in);

    %% Mass Balance
    % Equivalent change in fluid volume in WEC-driven pump chambers
    dp = 1e2;
    cap = @(p,V) deadVCap(p,V,par);
    Va = @(theta) par.V_wecDead + par.D_WEC*(par.theta_max - theta);
    Vb = @(theta) par.V_wecDead + par.D_WEC*(par.theta_max - theta);

     % Chamber 'a'
    cap1 = @(p) cap(p,Va(out.theta(1)));
    capend = @(p) cap(p,Va(out.theta(end)));
    deltaV_a = deltaV_NI(out.par.p_o,out.p_a(end),capend,dp) ...
                - deltaV_NI(out.par.p_o,out.p_a(1),cap1,dp);

     % Chamber 'b'
    cap1 = @(p) cap(p,Vb(out.theta(1)));
    capend = @(p) cap(p,Vb(out.theta(end)));
    deltaV_b = deltaV_NI(out.par.p_o,out.p_a(end),capend,dp) ...
                - deltaV_NI(out.par.p_o,out.p_a(1),cap1,dp);

    % Change in fluid volume in accumulators
    dp = 1e2;
     % Low-pressure
    cap = @(p) capAccum(p,par.pc_l,par.Vc_l,par.f,par);
    deltaV_l = deltaV_NI(out.p_l(1),out.p_l(end),cap,dp);

     % High-pressure outlet of WEC-driven pump
    cap = @(p) capAccum(p,par.pc_h,par.Vc_h,par.f,par);
    deltaV_h = deltaV_NI(out.p_h(1),out.p_h(end),cap,dp);

     % High-pressure inlet of RO module
    cap = @(p) capAccum(p,par.pc_ro,par.Vc_ro,par.f,par);
    deltaV_ro = deltaV_NI(out.p_ro(1),out.p_ro(end),cap,dp);

    % Total change in stored volume in the system
    out.deltaV_total = deltaV_a + deltaV_b ...
                     + deltaV_l + deltaV_h + deltaV_ro;

    % Flow at boundaries
    qbnds = out.q_c ...
            - (out.q_perm + out.q_brine + out.q_lPRV);

    % Total balance of flow: flow in minus change in volume stored
    out.Vbal = trapz(out.t,qbnds) - out.deltaV_total;
    out.Vbal_error = out.Vbal/trapz(out.t,out.q_perm);

%% %%%%%%%%%%%%   FUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function dydt = sys(t,y,par)
    	[dydt, ~ , ~ ] = sys_prototype(t,y,par);
    end

    function [dydt, nonState, control] = sysPost(t,y,par)
    	[dydt, nonState, control] = sys_prototype(t,y,par);
    end

    function [dydt, nonState, control] = sys_prototype(t,y,par)
        [dydt, nonState, control] = sys_refPTO(t,y,par);
    end
end
