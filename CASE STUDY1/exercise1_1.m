% EMTH171
% Case_Study1: Exercise 1: hauling a sled up incline with cable
% mzh99/jzh200 | Jiyao Zhu & Menghao Zhan

clear
clc
close all

% Motor out power(W)/ Sharf speed (RPM)
% Known values
m = 1000;  % Sled mass (Kg)
r_g = 20;  % Gear radius (m) 
r_p = 0.5; % Pulley radius (m)
ang = pi/4;% The angle of slope
f = 0.2;   % Coefficent of freaction
g = 9.81;  % Gravity (m/s^2)
b = 314.16;% rad/s
%--------------------------------Function----------------------------------
p_o = @(v) (b*(r_g*v)/r_p)-((r_g.^2*v.^2)/(r_p.^2));
p_d = @(v) m*g*sin(ang)*v + f*m*g*cos(ang)*v;
%---------------------------------Values-----------------------------------
vMin = 0;
vMax = 1000;
vstep = 0.001;
xarray = vMin : vstep : vMax;
%---------------------------------Plotting---------------------------------
figure(2)
plot(xarray,p_o(xarray), 'DisplayName', 'Motor output/W');
hold on
plot(xarray,p_d(xarray), 'DisplayName', 'Total power demand/W');
hold off
title('Exercise 1')
xlabel('Sled speed/m/s')
ylabel('Power/W')
xlim([0,12])
legend