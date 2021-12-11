function drawAircraft(uu)
% AER1216 Fall 2021 
% Fixed Wing Project Code
%
% drawAircraft.m
%
% Aircraft drawing file, which draws the aircraft and its inertial frame 
% position in figure 1. Model updates at each timestep. Code structure 
% adapted from Small Unmanned Aircraft: Theory and Practice by R.W. Beard 
% and T. W. McLain. 
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
% t                 time [s]
%
% Outputs:
% aircraft position plot in figure 1. 
%
% Last updated: Pravin Wedage 2021-11-09

% Function Inputs
px = uu(1); % x position   
py = uu(2); % y position
pz = uu(3); % z position
u = uu(4); % x velocity
v = uu(5); % y velocity
w = uu(6); % z velocity
phi = uu(7); % roll angle         
theta = uu(8); % pitch angle     
psi = uu(9)- pi/2; % yaw angle, adjusted for plotting purposes     
p = uu(10); % roll angular velocity
q = uu(11); % pitch angular velocity
r = uu(12); % yaw angular velocity
t = uu(13); % time

% Air Vehicle Parameters
fuse_l1 = 0.2;      % COG to nose tip
fuse_l2 = 0.15;     % COG to nose base
fuse_l3 = 0.6;      % COG to tail [m]
fuse_h  = 0.15;     % fuselage height
fuse_w  = 0.15;     % fuselage max width [m]
wing_l  = 0.18994;  % chord length [m]
wing_w  = 2.8956;   % wingspan [m]
tail_l  = 0.12;     % tailwing chord length [m]
tail_h  = 0.5;      % vertical stab. height [m]
tail_w  = 1;        % tailwing wingspan
geo = [fuse_l1,fuse_l2,fuse_l3,fuse_h,fuse_w,...
    wing_l,wing_w,tail_l,tail_h,tail_w]; % matrix of geometric parameters

% Defining Persistent Variables
persistent vehicle_pts;
persistent vehicle_handle;
persistent fig1 ax1;

% Plot variables
fig1_as = 2; % 1/2 of the axis length of the aircraft plot

% Initialize plot
if t <= 0.1  
    % Air Vehicle Plot
    fig1 = figure(1);
    clf
    ax1 = axes('Parent',fig1);
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    view(-58,47)
%     axis([-10,10,-10,10,-10,10]);
    axis(ax1,[px-fig1_as,px+fig1_as,...
              py-fig1_as,py+fig1_as,...
              -pz-fig1_as,-pz+fig1_as]);
%     axis([-2,2,-2,2,-2,2]);
    ax1.Box = 'on';
    ax1.XGrid = 'on';
    ax1.YGrid = 'on';
    ax1.ZGrid = 'on';
    set(ax1,'Ydir','reverse');
    vehicle_pts = defineVehiclePoints(geo);
    vehicle_handle = drawVehicleBody(vehicle_pts,px,py,pz,phi,theta,psi,[]);
    hold on

else % redraw system at every other timestep
    %redraw vehicle and update view to follow vehicle
    drawVehicleBody(vehicle_pts,px,py,pz,phi,theta,psi, vehicle_handle);
    axis(ax1,[px-fig1_as,px+fig1_as,...
              py-fig1_as,py+fig1_as,...
              -pz-fig1_as,-pz+fig1_as]);
    set(ax1,'Ydir','reverse');

end

end


%% Vehicle Body Drawing Function
% Modified slightly from textbook code
function handle = drawVehicleBody(NED,px,py,pz,phi,theta,psi,handle)
NED = rotate(NED,phi,theta,psi); % rotate air vehicle by phi, theta, psi
NED = translate(NED,px,py,pz); % translate air vehicle by px, py, pz 

% Plotting of air vehicle
if isempty(handle)
    % on startup, initial points will be defined. Each face of the aircraft
    % will be defined separately.
    figure(1); hold on
    nose_top = fill3(NED(1,[1,2,3]),NED(2,[1,2,3]),-NED(3,[1,2,3]),'k');
    nose_bot = fill3(NED(1,[1,4,5]),NED(2,[1,4,5]),-NED(3,[1,4,5]),'k');
    nose_lt = fill3(NED(1,[1,3,4]),NED(2,[1,3,4]),-NED(3,[1,3,4]),'k');
    nose_rt = fill3(NED(1,[1,2,5]),NED(2,[1,2,5]),-NED(3,[1,2,5]),'k');
    fuse_top = fill3(NED(1,[2,3,6]),NED(2,[2,3,6]),-NED(3,[2,3,6]),'k');
    fuse_bot = fill3(NED(1,[4,5,6]),NED(2,[4,5,6]),-NED(3,[4,5,6]),'k');
    fuse_lt = fill3(NED(1,[3,4,6]),NED(2,[3,4,6]),-NED(3,[3,4,6]),'k');
    fuse_rt = fill3(NED(1,[2,5,6]),NED(2,[2,5,6]),-NED(3,[2,5,6]),'k');
    wing = fill3(NED(1,[7,8,9,10]),NED(2,[7,8,9,10]),-NED(3,[7,8,9,10]),'w');
    tail = fill3(NED(1,[11,12,13,14]),NED(2,[7,8,9,10]),-NED(3,[7,8,9,10]),'w');
    fin = fill3(NED(1,[6,15,16]),NED(2,[6,15,16]),-NED(3,[6,15,16]),'w');
    handle=[nose_top,nose_bot,nose_lt,nose_rt,fuse_top,fuse_bot,...
        fuse_lt,fuse_rt,wing,tail,fin]; 
    % the previous two lines group all shapes into one patch (polygon) file
    % that will be subsequently edited. 
else
    % after startup, each polygon in the aircraft shape will be updated
    % individually. The position in the handle matrix corresponds to the
    % order shown above. 
    set(handle(1),'XData',NED(1,[1,2,3]),'YData',NED(2,[1,2,3]),'ZData',-NED(3,[1,2,3]));
    set(handle(2),'XData',NED(1,[1,4,5]),'YData',NED(2,[1,4,5]),'ZData',-NED(3,[1,4,5]));
    set(handle(3),'XData',NED(1,[1,3,4]),'YData',NED(2,[1,3,4]),'ZData',-NED(3,[1,3,4]));
    set(handle(4),'XData',NED(1,[1,2,5]),'YData',NED(2,[1,2,5]),'ZData',-NED(3,[1,2,5]));
    set(handle(5),'XData',NED(1,[2,3,6]),'YData',NED(2,[2,3,6]),'ZData',-NED(3,[2,3,6]));
    set(handle(6),'XData',NED(1,[4,5,6]),'YData',NED(2,[4,5,6]),'ZData',-NED(3,[4,5,6]));
    set(handle(7),'XData',NED(1,[3,4,6]),'YData',NED(2,[3,4,6]),'ZData',-NED(3,[3,4,6]));
    set(handle(8),'XData',NED(1,[2,5,6]),'YData',NED(2,[2,5,6]),'ZData',-NED(3,[2,5,6]));
    set(handle(9),'XData',NED(1,[7,8,9,10]),'YData',NED(2,[7,8,9,10]),'ZData',-NED(3,[7,8,9,10]));
    set(handle(10),'XData',NED(1,[11,12,13,14]),'YData',NED(2,[11,12,13,14]),'ZData',-NED(3,[11,12,13,14]));
    set(handle(11),'XData',NED(1,[6,15,16]),'YData',NED(2,[6,15,16]),'ZData',-NED(3,[6,15,16]));
    drawnow
end
end

%% Rotation Function
% Adopted entirely from textbook code
function pts=rotate(pts,phi,theta,psi);
  % define rotation matrix (right handed)
  psi = psi + deg2rad(90);
  R_roll = [...
          1, 0, 0;...
          0, cos(phi), sin(phi);...
          0, -sin(phi), cos(phi)];
  R_pitch = [...
          cos(theta), 0, -sin(theta);...
          0, 1, 0;...
          sin(theta), 0, cos(theta)];
  R_yaw = [...
          cos(psi), sin(psi), 0;...
          -sin(psi), cos(psi), 0;...
          0, 0, 1];
  R = R_roll*R_pitch*R_yaw;  
  
    % note that R above either leaves the vector alone or rotates
    % a vector in a left handed rotation.  We want to rotate all
    % points in a right handed rotation, so we must transpose
  R = R';
  % rotate vertices
  pts = R*pts;
end

%% Translation Function
function pts = translate(pts,px,py,pz)
% translates the air vehicle according to the inputs. 
pts = pts + repmat([px;py;pz],1,size(pts,2));
end

%% Vehicle Vertices Defining Function
function pts = defineVehiclePoints(geo)
% Defines the vertices of the air vehicle based on the provided geometric
% parameters.
pts = [...
    geo(1), 0, geo(4)*0.1;...            % pt 1 - nose tip
    geo(2), geo(5)*0.5, -geo(4)*0.5;...  % pt 2 - nose top right
    geo(2), -geo(5)*0.5, -geo(4)*0.5;... % pt 3 - nose top left
    geo(2), -geo(5)*0.5, geo(4)*0.5;...  % pt 4 - nose bottom left
    geo(2), geo(5)*0.5, geo(4)*0.5;...   % pt 5 - nose bottom right
    -geo(3), 0, 0;...                    % pt 6 - fuselage tail
    0, geo(7)*0.5, 0;...                 % pt 7 - wing front right
    -geo(6), geo(7)*0.5, 0;...           % pt 8 - wing back right
    -geo(6), -geo(7)*0.5, 0;...          % pt 9 - wing back left
    0, -geo(7)*0.5, 0;...                % pt 10 - wing front left
    -geo(3)+geo(8), geo(10)*0.5, 0;...   % pt 11 - tail front right
    -geo(3), geo(10)*0.5, 0;...          % pt 12 - tail back right
    -geo(3), -geo(10)*0.5, 0;...         % pt 13 - tail back left
    -geo(3)+geo(8), -geo(10)*0.5, 0;...  % pt 14 - tail front left
    -geo(3)+geo(8), 0, 0;...             % pt 15 - vert. front
    -geo(3), 0, -geo(9);...              % pt 16 - tail front left
    ]';
end