function [obj] = objectivefun(in)
ro=1.29;   % [kg*m^-3];
g=9.81;    % [kg*m]
d=.063;    % [m];
m=.05;     % [kg];
eta=1;

alfa=ro*pi*d^2/(8*m);

% max_p0x = 100;
% max_p0y = 100;
% max_p0z = 50;
% 
% max_v0 = 100;
% max_elev_angle = pi;
% max_plane_angle = 2*pi;
% max_w_x = 1000;
% max_w_y = 1000;
% max_w_z = 1000;

max_p0x = 1;
max_p0y = 1;
max_p0z = 1;

max_v0 = 1;
max_elev_angle = pi;
max_plane_angle = 2*pi;
max_w_x = 1;
max_w_y = 1;
max_w_z = 1;



p0x = in(1);
p0y = in(2);
p0z = in(3);
v0 = in(4);
elev_angle = in(5);
plane_angle = in(6);
wx = in(7);
wy = in(8);
wz = in(9);
cor = in(10);
cof = in(11);
% cor = 0.8;
% cof = 0.25;
% Load the track, var name: track:
track = [];
load('../ball_tracking/image_track2.mat');

total_sim_time = size(track,1) / 30; % in secs.
time_interval = 0.0333; % in secs (30 frames/s)
time=0:time_interval:total_sim_time;

if size(time , 2) > size(track,1)
    time = time(1,1:size(track,1));
end

% This is for this specific track
%time=0:.006:0.5280;

v0x = v0 * cos(plane_angle) * cos(elev_angle);
v0y = v0 * sin(plane_angle) * cos(elev_angle);
v0z = v0 * sin(elev_angle);

w = sqrt(wx^2+wy^2+wz^2);

ICs=[p0x , p0y , p0z , v0x, v0y , v0z]; % this contains the inital condition, i.e.,

u  = @(t,x)sqrt(x(4).^2+x(5).^2+x(6).^2);

CD = @(t,x)(0.508+(1./(22.053+4.196*(u(t,x)./w).^(5/2))).^(2/5));
CM = @(t,x)(1/(2.022+.981*(u(t,x)./w)));

topspin= @(t,x)([x(4); x(5); x(6);
    (-1)*CD(t,x)*alfa*x(4)*x(4) - (eta*CM(t,x)*alfa*(d/2)*(x(5)*wz - x(6)*wy));
    (-1)*CD(t,x)*alfa*x(5)*x(5) - (eta*CM(t,x)*alfa*(d/2)*(x(6)*wx - x(4)*wz));
    (-1)*g-CD(t,x)*alfa*x(6)*x(6)-(eta*CM(t,x)*alfa*(d/2)*(x(4)*wy - x(5)*wx))]);

[tts, XYZts]  = ode113(topspin, time, ICs, []);
%[tts, XYZts]  = ode23tb(topspin, time, ICs, []);


%% Time to bounce the ball.
% If the z coordinate of the ball is zero, bounce the ball.
nonzerosZ = XYZts(XYZts(:,3) >= 0);
nb_non_zeros = size(nonzerosZ , 1);
if nb_non_zeros < size(XYZts(:,3) , 1) && nb_non_zeros > 0
    % wee need to bounce of the ball in the remaining time.
    time_elapsed = time_interval *  nb_non_zeros;
    time_remain = total_sim_time - time_elapsed;
    
    % Find the new V0, and w0 using the bounce court properties.
    curr_para = XYZts(nb_non_zeros,:);
    curr_v0 = sqrt(curr_para(4)^2 + curr_para(5)^2 + curr_para(6)^2);
    curr_elev_angle = atan(curr_para(6) / sqrt(curr_para(4)^2 + curr_para(5)^2));
    curr_plane_angle = atan(curr_para(5) / curr_para(4));
    
    %simulate the bounce using tnsbouce:
    % convert angular velocity (rad/sec) to (rev/sec)
    w_rev = w/ (2*pi); 
    cor_trans = sqrt(cor*cor);
    cof_trans = sqrt(cof*cof);
    [Vbo,Wbo]=tnsbounce([curr_v0 curr_elev_angle curr_plane_angle],w_rev,[0 0 0],[pi/2 0],cor,cof,0);
    
    v0_new = Vbo(1);
    elev_angle_new = Vbo(2);
    plane_angle_new = Vbo(3);
    % now solve for the trajectory using the new parameters.
    v0x_new = v0_new * cos(plane_angle_new) * cos(elev_angle_new);
    v0y_new = v0_new * sin(plane_angle_new) * cos(elev_angle_new);
    v0z_new = v0_new * sin(elev_angle_new);
    
    w0_new = Wbo(1);
    elev_w_angle = Wbo(2);
    plane_w_angle = Wbo(3);
    wx_new = w0_new * cos(plane_w_angle) * cos(elev_w_angle);
    wy_new = w0_new * sin(plane_w_angle) * cos(elev_w_angle);
    wz_new = w0_new * sin(elev_w_angle);
    
    p0x_new = curr_para(1);
    p0y_new = curr_para(2);
    p0z_new = 0;
    
    ICs_new = [p0x_new p0y_new p0z_new v0x_new v0y_new v0z_new];
    
    time_new = 0:time_interval:time_remain;
    
    topspin_new= @(t,x)([x(4); x(5); x(6);
    (-1)*CD(t,x)*alfa*x(4)*x(4) - (eta*CM(t,x)*alfa*(d/2)*(x(5)*wz_new - x(6)*wy_new));
    (-1)*CD(t,x)*alfa*x(5)*x(5) - (eta*CM(t,x)*alfa*(d/2)*(x(6)*wx_new - x(4)*wz_new));
    (-1)*g-CD(t,x)*alfa*x(6)*x(6)-(eta*CM(t,x)*alfa*(d/2)*(x(4)*wy_new - x(5)*wx_new))]);
    
    [tts_after, XYZafterbounce]  = ode113(topspin_new, time_new, ICs_new, []);
    
    %now join XYZafterbounce with XYZts
    before = XYZts(1:nb_non_zeros,:);
    new = [before ; XYZafterbounce];
    
    if size(new , 1) > size(time,2)
        new = new(1:size(track,1),:);
    end
    
end

% From the simulated XYZ coordinates, get the image coordinates.
% TODO: In the function, will have to use the caliberated matrix to
% project to Image coordinates

%% Plot the current trajectory:

%simulateTrajectory(in);

%% Project the 3D world coordinates to 2D image coordinates.

% Projects the 3D coordinates to 2Dcoordinates.
if nb_non_zeros < size(XYZts(:,3) , 1) && nb_non_zeros > 0
    XYK_image_sim = getImagePoints(new);
else
    XYK_image_sim = getImagePoints(XYZts);
end


% TODO: Get the ground truth (tracked ball coordinates)

%XYK_image_ground = getBallCoordinates();
%XYK_image_ground = track;

% XYK, are in N x 3 shape (x,y,k)
obj = mean( (track(:,1) - XYK_image_sim(:,1)).^2 ...
         + (track(:,2) - XYK_image_sim(:,2)).^2);

disp(obj);
end