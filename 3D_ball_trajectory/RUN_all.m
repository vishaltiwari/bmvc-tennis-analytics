% Simulation models of a TENNIS ball in VACUUM and AIR with and without 
% spin using Symbolic MATH, MuPAD, ODEx solvers and SIMULINK
% This script executes all Simulink models and M-files except for MuPAD
% commands. 
% NB: place all scripts and SIMULINK models in the same directory, and then
% execute this script RUN_all.m
% NB: Plot figures require from a user to indicate graphically where to
% display final computed data in the plot area. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Written by Sulaymon L. ESHKABILOV, Ph.D
%                         October, 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
close all
%% MATLAB 2010 or earlier versions
%{
In VACUUM
x_z_vac=dsolve('D2x=0', 'D2z=-g','x(0)=0',...
'Dx(0)=v0*sin(eta*pi/180)','z(0)=h', ... 
'Dz(0)=v0*cos(eta*pi/180)', 't');
x=x_z_vac.x  %#ok
z=x_z_vac.z  %#ok
%}

%% MATLAB 2012 or later versions 
%{
In VACUUM
syms x(t) z(t) 
syms g h eta v0
Dx=diff(x); Dz=diff(z);
D2x=diff(x,2); D2z=diff(z,2);
 
x_z_vac=dsolve(D2x==0, D2z==-g,... 
    x(0)==0,Dx(0)==v0*sin(eta*pi/180),...
z(0)==h,Dz(0)==v0*cos(eta*pi/180));
x(t)=x_z_vac.x  %#ok
z(t)=x_z_vac.z  %#ok
%}
%% Analytical solution of the ball trajectory in vacuum
% Case #1
%% In MuPAD
%{
#ICs:

V0:=25; h:=1;eta:=15;g:=9.8;
ODE_eqn1:=ode({x''(t)=0, x(0)=0, x'(0)=V0*float(cos(eta*PI/180))}, x(t)); 
ODE_eqn2:=ode({z''(t)=-g, z(0)=h, z'(0)=V0*float(sin(eta*PI/180))},
z(t));
Sols1:=ode::solve(ODE_eqn1);
Sols2:=ode::solve(ODE_eqn2);
Sol_CURVE:=plot::Curve2d([24.14814566*t, 1.0 - t*(4.9*t - 6.470476128)],
t=0..2, ViewingBox = [0..35.5, 0..3.35]);
plot(Sol_CURVE, #G);
%}
%% In Symbolic MATH
% syms x(t) z(t) 
% syms g h eta v0
% Dx=diff(x); Dz=diff(z);
% D2x=diff(x,2); D2z=diff(z,2);
%
x_z_vac=dsolve('D2x=0', 'D2z=-g','x(0)=0',...
'Dx(0)=v0*cos(eta*pi/180)','z(0)=h', ... 
'Dz(0)=v0*sin(eta*pi/180)', 't');
xt=x_z_vac.x  %#ok
zt=x_z_vac.z  %#ok
g=9.8;        % gravitational accl 
h=1;          % ball hit 1 m above ground  
eta=15;       % ball hit under 15 degrees
v0=25;        % initial velocity of ball 
xt=vectorize(xt);
zt=vectorize(zt);
t=linspace(0,3.5, 200); 
xt=eval(xt);
zt=eval(zt);
z_i=find(zt==min(abs(zt)));
touch_gr=xt(z_i);
t_touch = t(z_i);
figure
plot(xt, zt,'ko', 'markersize', 4,'MarkerFaceColor','y'), grid on
xlim([0, touch_gr]);
text0=('Ball is hit h=1[m] above ground & under \eta=15^0');
gtext(text0, 'fontsize', 10);

text1=['Ball hits the ground in distance (of x) ' num2str(touch_gr) '[m]'];
gtext(text1);
text2=['Ball hits the ground after : ' num2str(t_touch)  ' [sec]'];
gtext(text2);


title('Trajectory of TENNIS ball hit in VACUUM, v_0=25 [m/s] ')
xlabel('x, (horizontal) [m]'), ylabel('z, [m]')
% 
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
h =1;      % [m];
v0=25;     % [m*sec^-1];
theta=15;  % [deg];
ro=1.29;   % [kg*m^-3];
g=9.81;    % [kg*m]
d=.063;    % [m];
m=.05;     % [kg];
w=20;      % [m*sec^-1]; angular velocity of a spinning ball
eta=1;     % describes direction(+-) of rotation; 
           % eta=1 is for topspin.
alfa=ro*pi*d^2/(8*m);
time=0:.01:1.5;
ICs=[0, 1, v0*cos(theta*pi/180), v0*sin(theta*pi/180)];

u  = @(t,x)sqrt(x(3).^2+x(4).^2);
CD = @(t,x)(0.508+(1./(22.053+4.196*(u(t,x)./w).^(5/2))).^(2/5));
CM = @(t,x)(1/(2.022+.981*(u(t,x)./w)));
vacuum = @(t,x)([x(3); x(4); 0; -g]);
nospin = @(t,x)([x(3); x(4);
    (-1)*CD(t,x)*alfa*u(t,x)*x(3);
    (-1)*g-CD(t,x)*alfa*u(t,x)*x(4)]);
topspin= @(t,x)([x(3); x(4);
    (-1)*CD(t,x)*alfa*u(t,x)*x(3)+eta*CM(t,x)*alfa*u(t,x)*x(4);
    (-1)*g-CD(t,x)*alfa*u(t,x)*x(4)-eta*CM(t,x)*alfa*u(t,x)*x(3)]);


[tvac, XZvac]= ode23(vacuum, time, ICs, []);
[tns, XZns]  = ode45(nospin, time, ICs, []);
[tts, XZts]  = ode113(topspin, time, ICs, []);

% It is also important to find when the ball hits the ground  
% (x-axis) and after how many seconds.
% For three cases: Vacuum, nospin and top-spin

% Case #1. Vacuum 
z_i=find(abs(XZvac(:,2))<=min(abs(XZvac(:,2))));
t_gr=XZvac(z_i,1);
t_t = time(z_i);

% Case #2. No-spin 
z1_i =find(abs(XZns(:,2))<=min(abs(XZns(:,2))));
t_gr1=XZns(z1_i,1);
t_t1 = time(z1_i);

% Case #3. Top-spin 
z2_i =find(abs(XZts(:,2))<=min(abs(XZts(:,2))));
t_gr2=XZts(z2_i,1);
t_t2 = time(z2_i);
figure 
plot(XZvac(:,1), XZvac(:,2), 'bo', 'linewidth', 1), grid
hold on
plot(XZns(:,1), XZns(:,2), 'm', 'linewidth', 2)
plot(XZts(:,1), XZts(:,2), 'ko-', 'linewidth', 1.0)
legend('In Vacuum','No-Spin in Air','Top-Spin in Air',0)
ylim([0, 3.35])
title 'Trajectory of a Tennis Ball hit under \eta=15^0, h=1 m'
xlabel 'Distance, [m]', ylabel 'Height, [m]'

tt1=['VACUUM: Ball hits the ground (x-axis):' num2str(t_gr) '[m]'];
gtext(tt1);
tt1a=['VACUUM: Ball hits the ground after: ' num2str(t_t) '[sec]'];
gtext(tt1a);
tt2=['Nospin: Ball hits the ground: ' num2str(t_gr1) '[m]'];
gtext(tt2);
tt2a=['Nospin: Ball hits the ground: ' num2str(t_t1)  ' [sec]'];
gtext(tt2a);
tt3=['Topspin: Ball hits the ground: ' num2str(t_gr2)  '[m]'];
gtext(tt3);
tt3a=['Topspin: Ball hits the ground: ' num2str(t_t2) '[sec]'];
gtext(tt3a);
hold off
%
%% Simulink model
eta=1;
open('TENNIS_Ball.mdl')
set_param('TENNIS_Ball','Solver','ode3','StopTime','1.5');
[tout, SIMout]=sim('TENNIS_Ball.mdl');
% plot(tout, SIMout(:,1), 'gd-', tout, SIMout(:,3),'ko', tout, SIMout(:,5),'rh'); 
% title('Simulation of the Simulink model')
% xlabel('time, [sec]'), ylabel('Height, [m]')
% ylim([0, 3.35])
% legend('Vacuum','NO-spin','TOP-spin',0)
figure
plot(SIMout(:,2), SIMout(:,1),'gd-',SIMout(:,4),SIMout(:,3),'ko',SIMout(:,6), SIMout(:,5),'rh')
ylim([0, 3.35])
title('Simulation of the Simulink model')
xlabel('ground, [m]'), ylabel('Height, [m]')
legend('Vacuum','TOP-spin','NO-spin',0)
shg
%% Spin value influence on distance covered by the ball
%{
clear all
open('TENNIS_Ball.mdl')
set_param('TENNIS_Ball','Solver','ode4','StopTime','1.5');
ETA=0:.2:1;
% Plot label settings;
Labelit = {};
Colorit = 'brgckmygrckmbgrygr';
Lineit  = '--:-:--:-:--:----:----:--';
Markit  = 'oxd+s*h+^v<p>.xsh+od+*^v';
figure 
for ii=1:length(ETA)
    eta=ETA(ii);
    [tout{ii},SIMout{ii}]=sim('TENNIS_Ball.mdl');
    Stylo = [Colorit(ii) Lineit(ii) Markit(ii)];
    plot(SIMout{:,ii}(:,4), SIMout{:,ii}(:,3), Stylo)
    Labelit{ii} = ['\eta = ' num2str(eta)];
    legend(Labelit{:},3)
    hold on
end
title(['Trajectory of a Ball with different spin values. ',...
    '\eta=1 (highest spin), \eta=0 (no spin)'])
ylim([0, 3.0])
xlabel('ground, [m]'), ylabel('Height, [m]')
hold off
shg
%}
%% Alternative version and more compact version.
clear all
open('TENNIS_Ball.mdl')
set_param('TENNIS_Ball','Solver','ode4','StopTime','1.5');
ETA=0:.2:1;
% Plot label settings;
Labelit = {};
Colorit = 'brgckmygrckmbgrygr';
Lineit  = '--:-:--:-:--:----:----:--';
Markit  = 'oxd+s*h+^v<p>.xsh+od+*^v';
figure 
for ii=1:length(ETA)
    eta=ETA(ii);
    [tout,SIMout]=sim('TENNIS_Ball.mdl');
    Stylo = [Colorit(ii) Lineit(ii) Markit(ii)];
    plot(SIMout(:,4), SIMout(:,3), Stylo)
    Labelit{ii} = ['\eta = ' num2str(eta)];
    legend(Labelit{:},3)
    hold on
end
title(['Trajectory of a Ball with different spin values. ',...
    '\eta=1 (highest spin), \eta=0 (no spin)'])
ylim([0, 3.0])
xlabel('ground, [m]'), ylabel('Height, [m]')
hold off
shg