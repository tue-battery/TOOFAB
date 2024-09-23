%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example showing how the toolbox can be used to simulate different ageing
% phenomena. The Cycling function exhibits CC discharge and CC-CV charge.
% The ageing modeling implementations are a work in progress and are
% subject to change. 
% This file is a part of the TOOlbox for FAst Battery simulation (TOOFAB)
% Github: https://github.com/tue-battery/TOOFAB
%
% Author: Francis le Roux (f.a.l.le.roux@tue.nl)
%
% TOOFAB is licensed under the BSD 3-Clause License
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear all; close all; 
load params
%%
%This implementation uses the default set of model parameters for the DFN
%model.
% Ageing On / Off
p_25.ageing = 1;

% Full Verbosity
p_25.verbose = 2; 

% Set Simplifications (s1 - s5 ON)
p_25.set_simp = [1,1,1,1,1,0]; 

% Init Some Variables
p_25.gamma_neg_prev = 0; 
p_25.gamma_pos_prev = 0;
p_25.epss_neg_prev = 0;
p_25.epss_pos_prev = 0;

% Time parameters
p_25.dt = 10;
p_25.SaveDataTskip = 1; % How many data points to skip when saving data
p_25.TotCycle = 10; % Total Cycles

% Conditional Parameters (Used for UCV Experiments)
p_25.VOffset = 0;

% Default Tuning Parameters
p_25.i02 =3.8703e-07;% SEI exchange current density
p_25.i0lpl = 0.0063; % Lithium Plating exchange current density
p_25.i0cath = 8e-05; % Cathode Dissolution exchange current density
p_25.beta_pc_sei = -0.0001; % PC caused from SEI Rate
p_25.beta_pc_lpl = -0.0001; % PC caused from LPL Rate
p_25.beta_cath = -0.001; % Rate of Change in volume fraction positive AM due to dissolution reactoin
p_25.lam_1 = 2.58e-5; % LAM
p_25.lam_2 = -4.7525e-07; %LAM 
p_25.beta_sc2= 0.1; % Surface Crack
p_25.beta_sc1 = 0.0004; % Surface Crack
p_25.A_max = 0.05; % Surface Crack
p_25.alpha_SEI=0.7;
p_25.alpha_lpl = 0.7; 
% Select ageing mechanisms (0 - Off / 1 - On)
p_25.ageMech.LiPlating = 0;
p_25.ageMech.SurfaceCrack = 0;
p_25.ageMech.PoreClog = 0;
p_25.ageMech.LAM = 0;
p_25.ageMech.Cathode = 0;
p_25.ageMech.SEI = 1;

%Choose C-rate for charging and discharging
p_25.Crate_charge = 1; 
p_25.Crate_discharge = 0.9;

%Choose step sizes for the various phases in the cycle. Note that the
%parameter values have no significant meaning, and is only meant to show
%how the additional dt output of the input_current function can be used
p_25.dt_charge = 1;
p_25.dt_discharge = 1;
p_25.dt_CV = 1;
p_25.T_amb=273.15+25;

% Run TOOFAB DFN 
out = DFN_Ageing_Cycle2(@cycle,2e7,1,p_25);


%% Cycle Function
function [i_app_next,mem,end_simulation,dt] = cycle(k,t,i_app,V,soc,cs,ce,phis,phie,mem,p)
Cbat = 1*p.Cbat/3600; %Capacity of the battery
V_max = p.OCV_max+p.VOffset;   %Voltage at which the battery switches to CV mode
V_min = p.OCV_min;   %Voltage after which the battery stops discharging
i_min = 0.05*Cbat;  %Current level to terminate the CV mode with
Kp1 = 30;   %Proportional gain of the PI controller for the CV stage
Ki1 = 10;   %Integral gain for the PI controller
Kd1 = 0;
end_simulation=0; %flag to indicate whether the simulation should end after the iteration. Has to be set at every iteration

% Initialize the various parameters. Note how the mem struct is used to
% store the states of the PI controller 
if k==1
    mem.status = 0;
    mem.e = 0;
    mem.e2 = 0;
end

% Condition to terminate the initial CC discharge stage and start the CC
% charge stage
if V(k)<=V_min && mem.status==0
    mem.status = 1;  
end

% Condition to terminate the CC charge stage and start the CV stage
if V(k)>=V_max &&mem.status==1
    mem.status = 2; 
end

%Condition to terminate the simulation
if i_app(k)<i_min && mem.status==2 && mem.CycleNo == mem.TotCycle
    end_simulation=1; 
end

%Condition to reset CC-CV Charge Cycle 
if i_app(k)<i_min && mem.status==2
    mem.Cbat_aged(mem.CycleNo) = mem.Cbat_aged_current;
    mem.CycleTimeStamp(mem.CycleNo) = k;
    mem.CycleNo = mem.CycleNo + 1;
    mem.status = 0;
    mem.e = 0;
    mem.e2 = 0;
end


%Switch between the different stages
switch mem.status
    case 0 %CC discharge
        i_app_next = -p.Crate_discharge*Cbat; %value of i_app for the next iteration
        % i_app_next = -7.5;
        dt = p.dt_discharge; %assign the chosen step size
    case 1 %CC charge
        i_app_next = p.Crate_charge*Cbat; %value of i_app for the next iteration
        % i_app_next = 5;
        dt = p.dt_charge; %assign the chosen step size
    case 2 %CV
        %Solve the PI controller system (in this case with a forward Euler
        %discretization
        mem.e(k+1) = V_max-V(k); 
        mem.e2(k+1) = (mem.e(k+1)-mem.e(k))/p.dt;
        mem.e3(k+1) = (mem.e2(k+1) - mem.e2(k))/p.dt;
        i_app_next = i_app(k)+p.dt*(Kp1*mem.e2(k+1)+Ki1*mem.e(k+1) + Kd1*mem.e3(k+1));
        dt = p.dt_CV; %assign the chosen step size
end
end