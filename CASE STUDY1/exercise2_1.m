%EMTH171
%Case_Study1: Exercise 1: hauling a sled up incline with cable
%mzh99

clear
clc
close all

% Motor out power(W)/ Sharf speed (RPM)
% Known values
a = 420;
b = 0.44;
r_f = 3.50;
r_g = 0.8;
r_w = 0.205;
%--------------------------------Function----------------------------------
p_o = @(w) (a*(w/(30/pi)))-(b*(w/(30/pi)).^2);
%---------------------------------Values-----------------------------------
vMin = 0;
vMax = 7000;
vstep = 1;
xarray = vMin : vstep : vMax;
%---------------------------------Plotting---------------------------------
figure(1)
yarray = p_o(xarray);
plot(xarray,yarray/1000);
title('Mathrmatical model of engine power supply')
xlabel('Engine speed/rpm')
ylabel('Power output/kW')
xlim([0,8000])
