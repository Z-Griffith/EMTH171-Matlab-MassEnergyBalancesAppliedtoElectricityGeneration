% EMTH171  29/8/2019
% Case Study 1 Exercise 2 
% Name: Hannah Tiong, Lixuan Feng
% For a petrol powered car, find the maximum speed for which power demand
% matches power output for four different cases of differing gradients and
% acceleration ot constant speed
 
clear
clc
close all

% Variables 
m = 1500;  % mass of car (kg)
Cd= 0.3;   % drag coefficient
A = 2;     % frontal area  (m^2) 
Crr = 0.01;% rolling resistance coefficient
R = 0.205; % wheel radius (m)
a = 420;   % alpha coefficient (W*s/rad)
b = 0.440; % beta coefficient (W*s^2/rad^2)
Rfd = 3.5; % final drive ratio 
Rgb = 0.8; % gearbox ratio
g = 9.81;  % gravity (m/s^2)
p = 1.2;   % air density (kg/m^3)

% Case variables, will differ depending on case
ang = 0; % gradient angle of road (rad)
acc = 0; % acceleration (m/s^2)


% Function for power demands 

% Pdrag = function for power demand from drag resistance 
Pdrag = @(V) Cd*A*(1/2)*p.*V.^3;  

% Prolling = function for power demand from rolling resistance
Prolling = @(V) Crr*m*g*cos(ang).*V;

% Pgradient = function for power demand from gradient
if ang == 0    
    Pgradient = @(V) 0; % Pgradient = 0, when gradient is 0
    PgExist = 0;        % Pgradient doesn't exist
else
    Pgradient = @(V) m*g*sin(ang).*V; % Pgradient = mgsin(x)V, when gradient != 0
    PgExist = 1;                      % Pgradient exist
end

% Paccelerate = function for power demand from acceleration

if acc == 0
    Paccelerate = @(V) 0; % Paccelerate = 0, when acceleration is 0
    PaExist = 0;          % Paccelerate doesn't exist
else
    Paccelerate = @(V) m*acc.*V;     %Pacclerate = maV, when gradient != 0
    PaExist = 1;                     % Paccelerate does exist
end

%Pdemand = total power demand
Pdemand = @(V)  Pdrag(V) + Pgradient(V) + Prolling(V) + Paccelerate(V);

% Poutput = Power output of petrol engine as a function of speed 
K = R/(Rfd*Rgb); % K is a constant, where radius of car wheel is divided by final drive and gear box gear ratio
Poutput = @(V) (a.*V./K) - ((b.*V.^2)/K^2);


% The following function is use to evaluate the speed at which the 
% power demand will match the power output of the petrol engine.
% f(V) = Pdemand - Poutput, where Pdemand = total power demand 
% and Poutput = power output
% The speed at which this function evaluates to zero is when power demand
% matches power output
% f(V) must differentiated to d(v) so that Newton's method can be use to 
% find the root of f(v)and hence the speed

f = @(V) Cd*A*(1/2)*p.*V.^3 + m*g*sin(ang).*V + Crr*m*g*cos(ang).*V ...
       + m*acc.*V - (a.*V./K) + ((b.*V.^2)./K^2);
d = @(V) 3*Cd*A*(1/2)*p.*V.^2 + m*g*sin(ang) + Crr*m*g*cos(ang) ...
       + m*acc - a/K + ((2*b.*V)./K^2);
      
% Variable for Newtons Method
N = 100;    % 100 iterations
v = 45;     % initial guess
tol = 1e-4; % tolerance value

% Newtons Method Iteration
for i = 1: N
  v(i+1) = v(i) - f(v(i))/d(v(i)); 
  diffX(i) = abs(v(i+1) -v(i));
  if diffX(i) < tol;  %, Newtons method stop iterating when difference of V approximate root < tol
      disp(v(i+1));   % display V, the road speed where power demand matches power output
      disp(i);        % display the amount of iterations needed to obtain final V solution
      break
  elseif i == N       
      disp(v(i+1));
      disp(i);
  end
end

% Values for plotting the graph

% Array of values for the x-axis
x = 0 : 0.1 : 100;
% Array of values for the y-axis, there is multiple functions of power that
% needs to be plot on the same graph
PdragArray = Pdrag(x);             % f(x) for power demand from drag resistance
PgradientArray = Pgradient(x);     % f(x) for power demand from gradient 
ProllingArray = Prolling(x);        % f(x) for power demand from rolling resistance
PaccelerateArray = Paccelerate(x);  % f(x) for power demand from acceleration
PdemandArray = Pdemand(x);          % f(x) for total power demand
PoutputArray = Poutput(x);          % f(x) for power output

% plotting the functions for power demand and power output
figure(1)
plot(x, PdragArray,'DisplayName', 'Drag');
hold on
plot(x,ProllingArray, 'DisplayName', 'Rolling Resistance');
plot(x,PdemandArray, 'DisplayName', 'Total Power demand ');
plot(x,PoutputArray, 'DisplayName', 'Power output');
hold off

   if PgExist   % plot gradient power demand if gradient !=0
       hold on 
       plot(x, PgradientArray, 'DisplayName','Gradient');
       hold off
  end
   if PaExist   % plot acceleration power demand if acceleration !=0
       hold on
       plot (x, PaccelerateArray, 'DisplayName','Acceleration');
       hold off
   end
 legend
ylim([-8e4 120000]);
xlabel('Car speed (m/s)')
ylabel('Power (W)')


% Use a plot of V verus iteration number(a "convergence plot")
% values for plotting
iteArray = 0 : 1 : i;   % iteration array for x-axis
yArray = v;   % x is a size (1XN) array of V approximation roots for N iterations 

% plotting convergence plot
figure(2)
plot(iteArray, yArray,'-o')
ylabel('V (m/s)')
xlabel('Number of iterations')
