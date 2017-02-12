
%% constants for air density, ball mass, etc
ro=1.29;   % [kg*m^-3];
g=9.81;    % [kg*m]
d=.063;    % [m];
m=.05;     % [kg];

wx = 0; %learn this parameter
wy = -100;  %learn this
wz = 100;  %learn this
w = sqrt(wx^2+wy^2+wz^2); 

v0 = 25.034;  %learn
elev_angle = 9*(pi/180); %learn
plane_angle = 10 * (pi/180); %learn
v0x = v0 * cos(plane_angle) * cos(elev_angle);
v0y = v0 * sin(plane_angle) * cos(elev_angle);
v0z = v0 * sin(elev_angle);

p0x = 0;    %learn
p0y = 5;    %learn
p0z = 1;    %learn

simulated_fun = @(p0x,p0y,p0z,v0,elev_angle,plane_angle,wx,wy,wx)