function aircraft_dynamics(block)
%MSFUNTMPL_BASIC A Template for a Level-2 MATLAB S-Function
%   The MATLAB S-function is written as a MATLAB function with the
%   same name as the S-function. Replace 'msfuntmpl_basic' with the 
%   name of your S-function.
%
%   It should be noted that the MATLAB S-function is very similar
%   to Level-2 C-Mex S-functions. You should be able to get more
%   information for each of the block methods by referring to the
%   documentation for C-Mex S-functions.
%
%   Copyright 2003-2010 The MathWorks, Inc.

% AER1216 Fall 2021 
% Fixed Wing Project Code
%
% aircraft_dynamics.m
%
% Fixed wing simulation model file, based on the Aerosonde UAV, with code
% structure adapted from Small Unmanned Aircraft: Theory and Practice by 
% R.W. Beard and T. W. McLain. 
% 
% Inputs: 
% delta_e           elevator deflection [deg]
% delta_a           aileron deflection [deg]
% delta_r           rudder deflection [deg]
% delta_t           normalized thrust []
%
% Outputs:
% pn                inertial frame x (north) position [m]
% pe                inertial frame y (east) position [m]
% pd                inertial frame z (down) position [m]
% u                 body frame x velocity [m/s]
% v                 body frame y velocity [m/s]
% w                 body frame z velocity [m/s]
% phi               roll angle [rad]
% theta             pitch angle [rad]
% psi               yaw angle [rad]
% p                 roll rate [rad/s]
% q                 pitch rate [rad/s]
% r                 yaw rate [rad/s]
%
% Last updated: Pravin Wedage 2021-11-09

%% TA NOTE
% The code segements you must modify are located in the derivatives
% function in this .m file. Modify other sections at your own risk. 


%
% The setup method is used to set up the basic attributes of the
% S-function such as ports, parameters, etc. Do not add any other
% calls to the main body of the function.
%
setup(block);

end 


%% Function: setup ===================================================
% Abstract:
%   Set up the basic characteristics of the S-function block such as:
%   - Input ports
%   - Output ports
%   - Dialog parameters
%   - Options
%
%   Required         : Yes
%   C-Mex counterpart: mdlInitializeSizes
%
function setup(block)

% Register number of ports
block.NumInputPorts  = 1;
block.NumOutputPorts = 1;

% Setup port properties to be inherited or dynamic
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% Override input port properties
for i = 1:block.NumInputPorts
    block.InputPort(i).Dimensions        = 4;
    block.InputPort(i).DatatypeID  = 0;  % double
    block.InputPort(i).Complexity  = 'Real';
    block.InputPort(i).DirectFeedthrough = false; % important to be false 
end

% Override output port properties
for i = 1:block.NumOutputPorts
    block.OutputPort(i).Dimensions       = 12;
    block.OutputPort(i).DatatypeID  = 0; % double
    block.OutputPort(i).Complexity  = 'Real';
%     block.OutputPort(i).SamplingMode = 'Sample';
end

% Register parameters
block.NumDialogPrms     = 1;
P = block.DialogPrm(1).Data; % must duplicate this line in each function

% Register sample times
%  [0 offset]            : Continuous sample time
%  [positive_num offset] : Discrete sample time
%
%  [-1, 0]               : Inherited sample time
%  [-2, 0]               : Variable sample time
block.SampleTimes = [0 0];

% Register multiple instances allowable
% block.SupportMultipleExecInstances = true;

% Register number of continuous states
block.NumContStates = 12;

% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'CustomSimState',  < Has GetSimState and SetSimState methods
%    'DisallowSimState' < Error out when saving or restoring the model sim state
block.SimStateCompliance = 'DefaultSimState';

% -----------------------------------------------------------------
% The MATLAB S-function uses an internal registry for all
% block methods. You should register all relevant methods
% (optional and required) as illustrated below. You may choose
% any suitable name for the methods and implement these methods
% as local functions within the same file. See comments
% provided for each function for more information.
% -----------------------------------------------------------------

% block.RegBlockMethod('PostPropagationSetup',    @DoPostPropSetup); % discrete states only
block.RegBlockMethod('SetInputPortSamplingMode', @SetInpPortFrameData);
block.RegBlockMethod('InitializeConditions',    @InitializeConditions);
% block.RegBlockMethod('Start',                   @Start); % Initialize Conditions is used
block.RegBlockMethod('Outputs',                 @Outputs); % Required
% block.RegBlockMethod('Update',                  @Update); % only required for discrete states
block.RegBlockMethod('Derivatives',             @Derivatives); % Required for continuous states
block.RegBlockMethod('Terminate',               @Terminate); % Required

end 


%% PostPropagationSetup:
%   Functionality    : Setup work areas and state variables. Can
%                      also register run-time methods here
%   Required         : No
%   C-Mex counterpart: mdlSetWorkWidths
%
function DoPostPropSetup(block)
block.NumDworks = 1;
  
  block.Dwork(1).Name            = 'x1';
  block.Dwork(1).Dimensions      = 1;
  block.Dwork(1).DatatypeID      = 0;      % double
  block.Dwork(1).Complexity      = 'Real'; % real
  block.Dwork(1).UsedAsDiscState = true;

end


%% InitializeConditions:
%   Functionality    : Called at the start of simulation and if it is 
%                      present in an enabled subsystem configured to reset 
%                      states, it will be called when the enabled subsystem
%                      restarts execution to reset the states.
%   Required         : No
%   C-MEX counterpart: mdlInitializeConditions
%
function InitializeConditions(block)

% Rename parameters
P = block.DialogPrm(1).Data; % must duplicate this line in each function

% Initialize continuous states
block.ContStates.Data(1) = P.pn0; 
block.ContStates.Data(2) = P.pe0;
block.ContStates.Data(3) = P.pd0;
block.ContStates.Data(4) = P.u0;
block.ContStates.Data(5) = P.v0;
block.ContStates.Data(6) = P.w0;
block.ContStates.Data(7) = P.phi0;
block.ContStates.Data(8) = P.theta0;
block.ContStates.Data(9) = P.psi0;
block.ContStates.Data(10) = P.p0;
block.ContStates.Data(11) = P.q0;
block.ContStates.Data(12) = P.r0;

end 

%% Start:
%   Functionality    : Called once at start of model execution. If you
%                      have states that should be initialized once, this 
%                      is the place to do it.
%   Required         : No
%   C-MEX counterpart: mdlStart
%
function Start(block)

block.Dwork(1).Data = 0;

end 

%% Input Port Sampling Method:
function SetInpPortFrameData(block, idx, fd)
  
  block.InputPort(idx).SamplingMode = 'Sample';
  for i = 1:block.NumOutputPorts
    block.OutputPort(i).SamplingMode  = 'Sample';   
  end
end

%% Outputs:
%   Functionality    : Called to generate block outputs in
%                      simulation step
%   Required         : Yes
%   C-MEX counterpart: mdlOutputs
%
function Outputs(block)

temp_mat = zeros(block.NumContStates,1); % thirteen states
for i = 1:block.NumContStates
     temp_mat(i) = block.ContStates.Data(i);
end

block.OutputPort(1).Data = temp_mat; % states

% for i = 1:block.NumOutputPorts
%     block.OutputPort(1).Data(i) = block.ContStates.Data(i);
% end

end 


%% Update:
%   Functionality    : Called to update discrete states
%                      during simulation step
%   Required         : No
%   C-MEX counterpart: mdlUpdate
%
function Update(block)

block.Dwork(1).Data = block.InputPort(1).Data;

end 


%% Derivatives:
%   Functionality    : Called to update derivatives of
%                      continuous states during simulation step
%   Required         : No
%   C-MEX counterpart: mdlDerivatives
%
function Derivatives(block)

% Rename parameters
P = block.DialogPrm(1).Data; % must duplicate this line in each function

% compute inertial constants
% K = ;
% k1 = ;
% k2 = ;
% k3 = ;
% k4 = ;
% k5 = ;
% k6 = ;
% k7 = ;
% k8 = ;

% map states and inputs
pn    = block.ContStates.Data(1);
pe    = block.ContStates.Data(2);
pd    = block.ContStates.Data(3);
u     = block.ContStates.Data(4);
v     = block.ContStates.Data(5);
w     = block.ContStates.Data(6);
phi   = block.ContStates.Data(7);
theta = block.ContStates.Data(8);
psi   = block.ContStates.Data(9);
p     = block.ContStates.Data(10);
q     = block.ContStates.Data(11);
r     = block.ContStates.Data(12);
delta_e = block.InputPort(1).Data(1)*pi/180 ; % converted inputs to radians
delta_a = block.InputPort(1).Data(2)*pi/180 ; % converted inputs to radians
delta_r = block.InputPort(1).Data(3)*pi/180 ; % converted inputs to radians
delta_t = block.InputPort(1).Data(4);

% Air Data 
Va = sqrt(u^2+v^2+w^2);
alpha = atan(w/u);
beta = asin(v/Va);

% rotation matrix
c_eb = [cos(theta)*cos(psi), sin(phi)*sin(theta)*cos(psi) - cos(phi) * sin(psi), cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi);...
        cos(theta)*sin(psi), sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi), cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi);...
        -sin(theta), sin(phi)*cos(theta), cos(phi)*cos(theta)];
% Aerodynamic Coefficients 
% compute the nondimensional aerodynamic coefficients here
X.u = (P.cxu + 2 * P.cxe) * .5 * rho * U_e * P.s;
X.w = c_x.alpha * .5 * rho * U_e * P.s;
X.q = c_x.q * .25 * rho * U_e * P.s * P.c;
X.u_dot = c_x.u_dot * .25 * rho * P.s* P.c;
X.w_dot = c_x.alpha_dot * .25 * rho * P.s * P.c;
X.q_dot = c_x.q_dot * .125 * rho * P.s * P.c^2;
X.delta_e = c_x.delta_e * .5 * rho * U_e^2 * P.s;
X.delta_p = c_x.delta_p * .5 * rho * U_e^2 * P.s;

Z.u = (c_z.u + 2 * c_z_e) * .5 * rho * U_e * P.s;
Z.w = c_z.alpha * .5 * rho * U_e * P.s;
Z.q = c_z.q * .25 * rho * U_e * P.s * P.c;
Z.u_dot = c_z.u_dot * .25 * rho * P.s * P.c;
Z.w_dot = c_z.alpha_dot * .25 * rho * P.s * P.c;
Z.q_dot = c_z.q_dot * .125 * rho * P.s * P.c^2;
Z.delta_e = c_z.delta_e * .5 * rho * U_e^2 * P.s;
Z.delta_p = c_z.delta_p * .5 * rho * U_e^2 * P.s;

M.u = (c_m.u + 2 * c_m_e) * .5 * rho * U_e * P.s * P.c;
M.w = c_m.alpha * .5 * rho * U_e * P.s * P.c;
M.q = c_m.q * .25 * rho * U_e * P.s * P.c^2;
M.u_dot = c_m.u_dot * .25 * rho * P.s * P.c^2;
M.w_dot = c_m.alpha_dot * .25 * rho * P.s * P.c^2;
M.q_dot = c_m.q_dot * .125 * rho * P.s * P.c^3;
M.delta_e = c_m.delta_e * .5 * rho * U_e^2 * P.s * P.c;
M.delta_p = c_m.delta_p * .5 * rho * U_e^2 * P.s * P.c;

Y.v = c_y.beta * .5 * rho * U_0 * P.s;
Y.p = c_y.p * .25 * rho * U_0 * P.s * b;
Y.r = c_y.r * .25 * rho * U_0 * P.s * b;
Y.v_dot = c_y.beta_dot * .25 * rho * P.s * b;
Y.p_dot = c_y.p_dot * .125 * rho * P.s * b^2;
Y.r_dot = c_y.r_dot * .125 * rho * P.s * b^2;
Y.delta_a = c_y.delta_a * .5 * rho * U_e^2 * P.s;
Y.delta_r = c_y.delta_r * .5 * rho * U_e^2 * P.s;

L.v = c_l.beta * .5 * rho * U_0 * P.s * b;
L.p = c_l.p * .25 * rho * U_0 * P.s * b^2;
L.r = c_l.r * .25 * rho * U_0 * P.s * b^2;
L.v_dot = c_l.beta_dot * .25 * rho * P.s * b^2;
L.p_dot = c_l.p_dot * .125 * P.s * b^3;
L.r_dot = c_l.r_dot * .125 * rho * P.s * b^3;
L.delta_a = c_l.delta_a * .5 * rho * U_e^2 * P.s * b;
L.delta_r = c_l.delta_r * .5 * rho * U_e^2 * P.s * b;

N.v = c_n.beta * .5 * rho * U_0 * P.s * b;
N.p = c_n.p * .25 * rho * U_0 * P.s * b^2;
N.r = c_n.r * .25 * rho * U_0 * P.s * b^2;
N.v_dot = c_n.beta_dot * .25 * rho * P.s * b^2;
N.p_dot = c_n.p_dot * .125 * rho * P.s * b^3;
N.r_dot = c_n.r_dot * .125 * rho * P.s * b^3;
N.delta_a = c_n.delta_a * .5 * rho * U_e^2 * P.s * b;
N.delta_r = c_n.delta_r * .5 * rho * U_e^2 * P.s * b;
% aerodynamic forces and moments
% compute the aerodynamic forces and moments here

% propulsion forces and moments
% compute the propulsion forces and moments here

% gravity
% compute the gravitational forces here

% total forces and moments (body frame)

% state derivatives
% the full aircraft dynamics model is computed here
pdot = [0 0 0];
pndot = pdot(1);
pedot = pdot(2);
pddot = pdot(3);

udot = X/P.m + r*v - q*w;
vdot = Y/P.m + p*w - r*u;
wdot = Z/P.m + q*u - p*v;

phidot = 0;
thetadot = 0;
psidot = 0;

omegadot = inv(P.i)*([L;M;N] - cross([p;q;r],P.i) * [p;q;r]);
pdot = omegadot(1);
qdot = omegadot(2);
rdot = omegadot(3);

% map derivatives
block.Derivatives.Data(1) = pndot;
block.Derivatives.Data(2) = pedot;
block.Derivatives.Data(3) = pddot;
block.Derivatives.Data(4) = udot;
block.Derivatives.Data(5) = vdot;
block.Derivatives.Data(6) = wdot;
block.Derivatives.Data(7) = phidot;
block.Derivatives.Data(8) = thetadot;
block.Derivatives.Data(9) = psidot;
block.Derivatives.Data(10)= pdot;
block.Derivatives.Data(11)= qdot;
block.Derivatives.Data(12)= rdot;

end 


%% Terminate:
%   Functionality    : Called at the end of simulation for cleanup
%   Required         : Yes
%   C-MEX counterpart: mdlTerminate
%
function Terminate(block)

end 

