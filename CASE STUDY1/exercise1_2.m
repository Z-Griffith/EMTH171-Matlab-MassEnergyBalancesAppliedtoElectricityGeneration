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
d_o = @(v) (b*r_g/r_p)-((2*r_g.^2*v)/(r_p.^2));
d_d = @(v) m*g*sin(ang) + f*m*g*cos(ang);
%------------------------------Interation----------------------------------
v = 2;        %Iniitial guess
N = 5;        %Guess times
store = [];   
store(1) = v ;
for ii = 1 : N
    v = v -((p_d(v)-p_o(v))/((d_d(v)-d_o(v))));
    store(ii+1) = v;
end

%---------------------------------Plotting---------------------------------
figure(2)
yarray = store;
xarray = [0:N];
plot(xarray,yarray,'-or');
title('Convergence plot: V versus number of iterations')
xlabel('Number of iterations')
ylabel('V(m/s)')
