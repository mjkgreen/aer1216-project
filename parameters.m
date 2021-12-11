% AER1216 Fall 2021 
% Fixed Wing Project Code
%
% parameters.m
%
% Initialization file which generates and stores all required data into the 
% structure P, which is then stored in the workspace. Simulink model calls 
% on this function at the start of every simulation. Code structure adapted
% from Small Unmanned Aircraft: Theory and Practice by R.W. Beard and T. W. 
% McLain. 
% 
% Inputs: 
% N/A
%
% Outputs:
% P                 structure that contains all aerodynamic, geometric, and
%                   initial condition data for the aircraft and simulation.
%
% Last updated: Pravin Wedage 2021-11-09

%% TA NOTE
% An easy way to store parameters for use in simulink is through the use of
% a structure. For example, P.g = 9.81 stores the value of gravitational
% acceleration in the field g that is contained within the structure P.
% Anytime P is called anywhere in the simulation code, the value of P.g is
% accessible. 

%% Parameter Computation
% Initial Conditions
clear all
% compute trim conditions            
P.Va0 = 15;         % initial airspeed (also used as trim airspeed)
P.Va_trim = 15; 
P.Va = P.Va_trim;


P.gravity = 9.81;
P.g = 9.81; 

% Aerosonde UAV Data
% physical parameters of airframe
P.m = 13.5;
P.ixx = 0.8244;
P.iyy = 1.135;
P.izz = 1.759;
P.ixz = 0.1204;
P.i = [P.ixx 0 P.ixz;
        0 P.iyy 0;
        P.ixz 0 P.izz];

P.s = 0.55;
P.c = 0.18994;

% aerodynamic coefficients
P.cl0 = 0.28;
P.cd0 = 0.03;
P.cm0 = -0.02338;
P.clalpha = 3.45;
P.cdalpha = 0.30;
P.cmalpha = -0.38;
P.clq = 0;
P.cdq = 0;
P.cmq = -3.6;
P.cldeltae = -0.36;
P.ddeltae = 0;
P.cmdeltae = -0.5;
P.epsilon = 0.1592;

P.cy0 = 0;
P.cl0 = 0;
P.cn0 = 0;
P.cybeta = -0.98;
P.clbeta = -0.12;
P.cnbeta = 0.25;
P.cyp = 0;
P.clp = -0.26;
P.cnp = 0.022;
P.cyr = 0;
P.clr = 0.14;
P.cnr = -0.35;
P.cydeltaa = 0;
P.cldeltaa = 0.08;
P.cndeltaa = 0.06;
P.cydeltar = -0.17;
P.cldeltar = 0.105;
P.cndeltar = -0.032;

% Control Input limits 
P.delta_e_max = deg2rad(45); % assumed symmetric
P.delta_a_max = deg2rad(45); 
P.delta_r_max = deg2rad(25);

% Initial Conditions % connects with aircraft_dynamics.m, do not modify
% structure field names
P.pn0    = 0;  % initial North position
P.pe0    = 0;  % initial East position
P.pd0    = -1000;  % initial Down position (negative altitude)
P.u0     = P.Va0; % initial velocity along body x-axis
P.v0     = 0;  % initial velocity along body y-axis
P.w0     = 0;  % initial velocity along body z-axisu_trim
P.phi0   = 0;  % initial roll angle
P.theta0 = 0;  % initial pitch angle
P.psi0   = 0;  % initial yaw angle
P.p0     = 0;  % initial body frame roll rate
P.q0     = 0;  % initial body frame pitch rate
P.r0     = 0;  % initial body frame yaw rate
P.delta_e0 =0;
P.delta_a0 =0;
P.delta_r0 =0;
P.delta_t0 =0;
                         
