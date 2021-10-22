% EMTH171
% Case_Study1: Exercise 2:
% Analyzing a car of petrol power supply system with several factors
% mzh99/jzh200 | Jiyao Zhu & Menghao Zhan


clear
clc
close all

% Motor out power(W)/ Sharf speed (RPM)
% Known values
a = 420;     % Engine parameters a (W.s/rad)
b = 0.44;    % Engine parameters b (W.s 2 /rad 2)
r_f = 3.50;  % Final drive ratio (m)
r_g = 0.8;   % Top gear with a gearbox ratio (m)
r_w = 0.205; % Wheel radius (m)
%--------------------------------Function----------------------------------
p_o = @(w) (a*(w/(30/pi)))-(b*(w/(30/pi)).^2);
%---------------------------------Values-----------------------------------
vMin = 0; % x_started point
vMax = 7000; % x_ended point
vstep = 1; % The value of every interval
xarray = vMin : vstep : vMax; 
%---------------------------------Plotting---------------------------------
figure(1)
yarray = p_o(xarray);
plot(xarray,yarray/1000);
title('Mathrmatical model of engine power supply')
xlabel('Engine speed/rpm')
ylabel('Power output/kW')
xlim([0,8000])

