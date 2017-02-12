%% Solve for 3D trajectory
% Trajectory equations:
% Given:    x" = -Cd*alfa*x'*x' + eta*Cm*alfa*(y'*)

%%   Solved by ODEx
% Trajectory of a Tennis Ball 
% Given:   x" = -Cd*alfa*v*x'+eta*Cm*alfa*v*z'
%          y" = -g-Cd*alfa*v*z'-eta*Cm*alfa*v*x'
%{
 x(0)=0; z(0)=0; x'(0)=v0*cos(eta);
 z'(0)=v0*sin(eta); eta=1; v0=25 [m/*sec^-1];
 eta=15 [deg]; ro=1.29 [kg*m^-3]; g=9.81 [kg*m];
 v=sqrt(x'^2+z'^2); d=.063 [m]; m=.05 [kg];
%}
clearvars
%% constants for air density, ball mass, etc
ro=1.29;   % [kg*m^-3];
g=9.81;    % [kg*m]
d=.063;    % [m];
m=.05;     % [kg];
eta=-1;     % describes direction(+-) of rotation; 
           % eta=1 is for topspin.
%%
% Different components of the anuglar velocity.
% Taking the assumption that w doesn't change.
wx = 0;
wy = -100;
wz = 100;
w = sqrt(wx^2+wy^2+wz^2);

% specify the speed, and angle.
v0 = 25.034;
elev_angle = 9*(pi/180);    % This is the angle between plane component, [vx vy] and vx from the horizontal.
plane_angle = 10 * (pi/180); % this is from the vertical in the plane. +ve counter clockwise.
v0x = v0 * cos(plane_angle) * cos(elev_angle);
v0y = v0 * sin(plane_angle) * cos(elev_angle);
v0z = v0 * sin(elev_angle);

% initial ball location.
p0x = 0;
p0y = 5;
p0z = 1;

alfa=ro*pi*d^2/(8*m);
total_sim_time = 1.5; % in secs.
time_interval = 0.0333; % in secs (30 frames/s)
time=0:time_interval:total_sim_time;

% ICs = [p0x , p0y , p0z , v0x , v0y , v0z];
ICs=[p0x , p0y , p0z , v0x, v0y , v0z]; % this contains the inital condition, i.e.,

u  = @(t,x)sqrt(x(4).^2+x(5).^2+x(6).^2);

CD = @(t,x)(0.508+(1./(22.053+4.196*(u(t,x)./w).^(5/2))).^(2/5));
CM = @(t,x)(1/(2.022+.981*(u(t,x)./w)));

%% Ball trajectory in vaccum.
vacuum = @(t,x)([x(4); x(5); x(6) ; 0 ; 0; -g]);

%% Ball trajectory in air, with no anuglar velocity
nospin = @(t,x)([x(4); x(5); x(6);
    (-1)*CD(t,x)*alfa*x(4)*x(4);
    (-1)*CD(t,x)*alfa*x(5)*x(5);
    (-1)*g-CD(t,x)*alfa*x(6)*x(6)]);

%% Ball trajectoy with angular velocity.
topspin= @(t,x)([x(4); x(5); x(6);
    (-1)*CD(t,x)*alfa*x(4)*x(4) - (eta*CM(t,x)*alfa*(d/2)*(x(5)*wz - x(6)*wy));
    (-1)*CD(t,x)*alfa*x(5)*x(5) - (eta*CM(t,x)*alfa*(d/2)*(x(6)*wx - x(4)*wz));
    (-1)*g-CD(t,x)*alfa*x(6)*x(6)-(eta*CM(t,x)*alfa*(d/2)*(x(4)*wy - x(5)*wx))]);


[tvac, XYZvac]= ode23(vacuum, time, ICs, []);
[tns, XYZns]  = ode45(nospin, time, ICs, []);
[tts, XYZts]  = ode113(topspin, time, ICs, []);

%% Time to bounce the ball.
% If the z coordinate of the ball is zero, bounce the ball.
nonzerosZ = XYZts(XYZts(:,3) >= 0);
nb_non_zeros = size(nonzerosZ , 1);
if nb_non_zeros < size(XYZts(:,3) , 1)
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
    [Vbo,Wbo]=tnsbounce([curr_v0 curr_elev_angle curr_plane_angle],w_rev,[0 0 0],[pi/2 0],0.8,0.25,0);
    
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
    
end

%% Plot the court lines.
% court dimensions in m, length, width, net height, service line :
Dx=23.7744;Dy=9.6;Dn=0.914;Ds=5.4864;

% lines
hold on
plot3(Dx*[0 1 1 0 0],Dy*[0 0 1 1 0],Dn*[0 0 0 0 0],'r'); % court
plot3(0.5*Dx*[1 1 1 1 1],Dy*[0 0 1 1 0],Dn*[0 1 1 0 0],'r'); % net
plot3(Ds*[1 1],Dy*[0 1],Dn*[0 0],'r'); % service line 1
plot3((Dx-Ds)*[1 1],Dy*[0 1],Dn*[0 0],'r'); % service line 2
plot3([Ds Dx-Ds],0.5*Dy*[1 1],Dn*[0 0],'r'); % half line
hold off

xlabel('x (m), ');
ylabel('y (m), ');
zlabel('h (m), ');
title('tennis ball trajectory vacuum , air(no ball spin), air(top spin) ');

% figure;
hold on
plot3(XYZvac(:,1), XYZvac(:,2), XYZvac(:,3), 'bo', 'linewidth', 0.5), grid
plot3(XYZns(:,1) , XYZns(:,2), XYZns(:,3), 'm', 'linewidth', 0.5)
%plot3(XYZts(:,1) , XYZts(:,2), XYZts(:,3), 'ko-', 'linewidth', 0.5)
plot3(new(:,1) , new(:,2), new(:,3), 'ko-', 'linewidth', 0.5)

rotate3d on;

