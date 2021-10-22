% EMTH171
% Case_Study1: Exercise 2:
% Analyzing a car of petrol power supply system with several factors
% mzh99/jzh200 | Jiyao Zhu & Menghao Zhan

clear
clc
close all

% Motor out power(W)/ Sharf speed (RPM)
% Known values
m = 1500;    % Mass of car (Kg)
c_d = 0.3;   % Drag coefficient
area = 2;    % Frontal area (m^2)
c_rr = 0.010;% Rolling resistance 
a = 420;     % Engine parameters a (W.s/rad)
b = 0.44;    % Engine parameters b (W.s 2 /rad 2)
r_f = 3.50;  % Final drive ratio (m)
r_g = 0.8;   % Top gear with a gearbox ratio (m)
r_w = 0.205; % Wheel radius (m)
g = 9.81;    % Gravity coefficient (m/s^2)
air = 1.2;   % Air density (Kg/m^3)
% ========================Changeable Values in 4 Cases=====================
angle = 0;   % Gradient (degree)
acc = 1;     % Acceleration (m/s^2)
K = r_w/(r_f*r_g); % K 
%--------------------------------Function----------------------------------
% Power demand from drag resistance
p_d = @(v) c_d.*area.*(1/2).*v.^3*air;
% Power demand from rolling resistance
p_rr = @(v) c_rr.*m*g.*cos(angle).*v;
% Power demand from gradient
p_r = @(v) m.*g.*sin(angle).*v;
% Power demand from acceleration
p_a = @(v) m.*acc.*v;
% Engine power
p_e = @(v) ((a*v)/ K) - (b*v.^2)/K.^2;
% Total power demand
p_t = @(v) p_d(v) + p_r(v) + p_rr(v) + p_a(v)
%---------------------------------Values-----------------------------------
vMin = 0;
vMax = 100;
vstep = 0.01;
xarray = vMin : vstep : vMax;
%---------------------------------Plotting---------------------------------
figure(1)
plot(xarray, p_d(xarray), 'DisplayName', 'Drag')
title('Road car power')
ylim([-8e4 150000]);
xlabel('Petrol powered car speed m/s')
ylabel('Power W')
hold on
plot(xarray, p_t(xarray), 'DisplayName', 'Total power demand')
plot(xarray, p_e(xarray), 'DisplayName', 'Engine power')
plot(xarray, p_rr(xarray), 'DisplayName', 'Rolling resistance')
hold off
legend

