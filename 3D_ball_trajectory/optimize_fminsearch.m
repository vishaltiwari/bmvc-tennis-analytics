wx_start = 0;
wy_start = -100;
wz_start = 100;

% specify the speed, and angle.
v0_start = 25.034;
elev_angle_start = -4*(pi/180);
plane_angle_start = -4 * (pi/180);

% p0x_start = 0;
% p0y_start = 5;
% p0z_start = 1;
p0x_start = 5.485;
p0y_start = 5.485;
p0z_start = 1;

cor = 0.8;
cof = 0.25;

ICs = [p0x_start,p0y_start,p0z_start,v0_start,elev_angle_start,plane_angle_start,wx_start,wy_start,wz_start,cor,cof];
%ICs = [p0x_start,p0y_start,p0z_start,v0_start,elev_angle_start,plane_angle_start,wx_start,wy_start,wz_start];
options = optimset('MaxFunEvals',Inf,'MaxIter',5000);

%[optimized_params,fval,exitflag,output] = fminsearch('objectivefun', ICs,options);
%options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective');
%[optimized_params,fval,exitflag,output] = fminsearch('objectivefun', ICs, options);

% Normalize the value parameter values:
max_p0x = 100;
max_p0y = 100;
max_p0z = 50;

max_v0 = 100;
max_elev_angle = pi;
max_plane_angle = 2*pi;
max_w_x = 1000;
max_w_y = 1000;
max_w_z = 1000;

% normalized_ICs = [p0x_start / max_p0x,p0y_start/max_p0y,p0z_start/max_p0z ...
%                     ,v0_start / max_v0 ,elev_angle_start / max_elev_angle, plane_angle_start/max_plane_angle ...
%                     ,wx_start/max_w_x ,wy_start/max_w_x ,wz_start/max_w_x,...
%                     cor,cof];

lb = [-Inf ,-Inf , -Inf , -Inf, -Inf , -Inf , -Inf , -Inf , -Inf , 0 ,0];
ub = [Inf , Inf , Inf , Inf , Inf , Inf , Inf , Inf , Inf , 1 , 1];
A=[];
b=[];
Aeq=[];
beq=[];
nonlcon = @objectivefun;
[optimized_params,fval,exitflag,output] = fmincon(nonlcon, ICs ,A,b,Aeq,beq,lb,ub,[],options);

%% visualize the trajectory.
% TODO:
simulateTrajectory(optimized_params);
