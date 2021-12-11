function plotmavstatevariables(uu)
% AER1216 Fall 2021 
% Fixed Wing Project Code
%
% plotmavstatevariables.m
%
% File that plots the state variables and control inputs as functions of
% time in figure 2. Model updates at each timestep. Code structure adapted 
% from Small Unmanned Aircraft: Theory and Practice by R.W. Beard and T. W.
% McLain. 
% 
% Inputs:
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
% delta_e           elevator deflection [deg]
% delta_a           aileron deflection [deg]
% delta_r           rudder deflection [deg]
% delta_t           normalized thrust []
% t                 time [s]
%
% Outputs:
% state and control variables as functions of time in figure 2. 
%
% Last updated: Pravin Wedage 2021-11-09

    % process inputs to function
    pn          = uu(1);             % North position (meters)
    pe          = uu(2);             % East position (meters)
    pd          = -uu(3);            % altitude (meters)
    u           = uu(4);             % body velocity along x-axis (meters/s)
    v           = uu(5);             % body velocity along y-axis (meters/s)
    w           = uu(6);             % body velocity along z-axis (meters/s)
    phi         = 180/pi*uu(7);      % roll angle (degrees)   
    theta       = 180/pi*uu(8);      % pitch angle (degrees)
    psi         = 180/pi*uu(9);      % yaw angle (degrees)
    p           = 180/pi*uu(10);     % body angular rate along x-axis (degrees/s)
    q           = 180/pi*uu(11);     % body angular rate along y-axis (degrees/s)
    r           = 180/pi*uu(12);     % body angular rate along z-axis (degrees/s)
    delta_e     = uu(13);     % elevator angle (degrees)
    delta_a     = uu(14);     % aileron angle (degrees)
    delta_r     = uu(15);     % rudder angle (degrees)
    delta_t     = uu(16);            % throttle setting (unitless)
    t           = uu(17);            % simulation time
    
    
    Va = sqrt(u^2 + v^2 + w^2);
    alpha = atan2(w,u);
    beta = asin(v/Va); 
    
    % compute course angle
%     chi = 180/pi*atan2(Va*sin(psi)+we, Va*cos(psi)+wn);
    
    % wrap angles
    phi = wrapTo180(phi);
    theta = wrapTo180(theta);
    psi = wrapTo180(psi);

    % define persistent variables 
    persistent pn_handle
    persistent pe_handle
    persistent pd_handle
    persistent u_handle
    persistent v_handle
    persistent w_handle
    persistent phi_handle
    persistent theta_handle
    persistent psi_handle
    persistent p_handle
    persistent q_handle
    persistent r_handle
    persistent delta_e_handle
    persistent delta_a_handle
    persistent delta_r_handle
    persistent delta_t_handle
    

  % first time function is called, initialize plot and persistent vars
    if t==0
        figure(2), clf

        subplot(8,2,1)
        hold on
        pn_handle = graph_y(t, pn, 'p_n', []);
        
        subplot(8,2,2)
        hold on
        u_handle = graph_y(t, u, 'u', []);
        
        subplot(8,2,3)
        hold on
        pe_handle = graph_y(t, pe, 'p_e', []);

        subplot(8,2,4)
        hold on
        v_handle = graph_y(t, v, 'v', []);
        
        subplot(8,2,5)
        hold on
        pd_handle = graph_y(t, pd, 'h', []);
        % altitude plot is reversed, pd here is altitude

        subplot(8,2,6)
        hold on
        w_handle = graph_y(t, w, 'w', []);

        subplot(8,2,7)
        hold on
        phi_handle = graph_y(t, phi, '\phi', []);
        
        subplot(8,2,8)
        hold on
        p_handle = graph_y(t, p, 'p', []);
        
        subplot(8,2,9)
        hold on
        theta_handle = graph_y(t, theta, '\theta', []);
        
        subplot(8,2,10)
        hold on
        q_handle = graph_y(t, q, 'q', []);
        
        subplot(8,2,11)
        hold on
        psi_handle = graph_y(t, psi, '\psi', []);
        
        subplot(8,2,12)
        hold on
        r_handle = graph_y(t, r, 'r', []);
        
        subplot(8,2,13)
        hold on
        delta_e_handle = graph_y(t, delta_e, '\delta_e', []);
        
        subplot(8,2,14)
        hold on
        delta_a_handle = graph_y(t, delta_a, '\delta_a', []);

        subplot(8,2,15)
        hold on
        delta_r_handle = graph_y(t, delta_r, '\delta_r', []);
        
        subplot(8,2,16)
        hold on
        delta_t_handle = graph_y(t, delta_t, '\delta_t', []);
        
    % at every other time step, redraw state variables
    else 
       graph_y(t, pn, 'pn', pn_handle);
       graph_y(t, pe, 'pe', pe_handle);
       graph_y(t, pd, 'h', pd_handle);
       graph_y(t, u, 'u', u_handle);
       graph_y(t, v, 'v', v_handle);
       graph_y(t, w, 'w', w_handle);
       graph_y(t, phi, '\phi', phi_handle);
       graph_y(t, theta, '\theta', theta_handle);
       graph_y(t, psi, '\psi', psi_handle);
       graph_y(t, p, 'p', p_handle);
       graph_y(t, q, 'q', q_handle);
       graph_y(t, r, 'r', r_handle);
       graph_y(t, delta_e, '\delta_e', delta_e_handle);
       graph_y(t, delta_a, '\delta_a', delta_a_handle);
       graph_y(t, delta_r, '\delta_r', delta_r_handle);
       graph_y(t, delta_t, '\delta_t', delta_t_handle);
    end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graph y with lable mylabel
function handle = graph_y(t, y, lab, handle)
  
  if isempty(handle)
    handle    = plot(t,y,'b','LineWidth',3);
    ylabel(lab)
    set(get(gca, 'YLabel'),'Rotation',0.0);
  else
    set(handle,'Xdata',[get(handle,'Xdata'),t]);
    set(handle,'Ydata',[get(handle,'Ydata'),y]);
    %drawnow
  end

%
%=============================================================================
% sat
% saturates the input between high and low
%=============================================================================
%
function out=sat(in, low, high)

  if in < low
      out = low;
  elseif in > high
      out = high;
  else
      out = in;
  end

% end sat  


