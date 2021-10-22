% EMTH171: Case Study 1 Exercise 2 mzh99/jzh200
% Analyzing a car of petrol power supply system with several factors
clear
clc
close all
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
acc = 0;     % Acceleration (m/s^2)
K = r_w/(r_f*r_g); % K 
%--------------------------------Function----------------------------------
f = @(v) c_d*area*(1/2)*air*v.^3 + m*g*sin(angle).*v + ...
    c_rr*m*g*cos(angle).*v + m*acc.*v - (a*v./K) + ((b*v.^2)/K^2);
d = @(v) 3*c_d*area*(1/2)*air*v.^2 + m*g*sin(angle) + ...
    c_rr*m*g*cos(angle) + m*acc - a/K + ((2*b*v)./K^2);
%---------------------------------Values-----------------------------------
N = 1000;    % 100 Iterations
v = 50;     % Initial guess
tol = 1e-4; % Tolerance value
%-----------------------------Processing-----------------------------------

for ii = 1 : N
    v(ii+1) = v(ii) - f(v(ii))/d(v(ii));
    if abs(v(ii+1)-v(ii)) < tol
        break
    end
end
%---------------------------------Plotting---------------------------------
xArray = 0 : 1 : ii;   %  X-axis
yArray = v            %  Y-axis 
figure(1)
plot(xArray, yArray,'-p')
ylabel('Velocity m/s')
xlabel('Number of Iterations')
