% EMTH171
% Case Study 2?
% Name: Hannah Tiong, Lixuan Feng?
% Solution Method

clear
clc

%Variables with assosciated differential equation

Hni0 = 356.55;  % initial height of NI lake (m)
Hsi0 = 406;    % inital height of SL lake (m)
Vni0 = 0;        % inital volume of NI lake (m^3)
Vsi0 = 0;        % inital volume of SI lake (m^3)

% fixed parameters
Ani = 620e6; % surface area of NL lake (m^2)
Asi = 350e6; % surface area of SL lake (m^2)
K = 1.55;   %   m^0.5/s
L = 300;     % length of spillway weir (m)
Hwni = 357.25; % height of NI weir     (m)
Hwsi = 410;    % height of SL weir     (m)
mean = 5000;    % mean for Normal distribution for North Island Power demand (h)
sd  = 1000;     % standard deviation (h)
t = 8760;
 
p = 998;   % kg/m^3
g = 9.81;    % gravity (m/s^2)
Hgen_ni = 80;   % height of NI generator(m)
Hgen_si = 0;    % heigth of SI generator(m)
% P_hydro = power;%%%%%%
Psi_hydro = 1 ;
% Equations
%Pni_demand = @(t) 4125 + 1.4e6 * normpdf(t, mean, sd); % power demand (MW) in North Island?
%Psi_demand = 1940;    % power demand (MW) in South Isaland?
f_Fni_in = @(t)  345 + 73* sin((2*pi*(t-3624))/8760); % inlet flow (m^3/s) North Island
f_Fni_gen = f_Fni_in;
f_Fni_spill = @(Hni) K*L*(Hni-Hwni)^1.5;                          


f_Fsi_in = @(t)  593 - 183* sin((2*pi*(t-2320))/8760); % inlet flow (m^3/s) South Island
f_Fsi_gen = @(t,Hni,Hsi) (Psi_hydro*10^6/0.9*p*g*(Hsi-Hgen_si)) ;
f_Fsi_spill = @(Hsi) (K*L*(Hsi-Hwsi)^1.5);


%dhni/dt
%f1 = @(t,Hni) ((f_Fni_in - f_Fni_gen - f_Fni_spill(Hni))/Ani) * 3600 ;*3600 for (m^3/s to m^3/h)
 f1 = @(t,Hni)   (f_Fni_spill(Hni)/Ani)*3600 ; % (m^3/h)
 %f1 =   @(t,Hni) (K*L*(Hni-Hwni)^1.5)*3600;

% dhsi/dt define on for no spill as well as if spill Hsi>Hwsi
f2a = @(t,Hni, Hsi) ((f_Fsi_in(t) - f_Fsi_gen(t,Hni,Hsi) - f_Fsi_spill(Hsi))/Asi)*3600; % (m^3/h)
f2b = @(t,Hni, Hsi) ((f_Fsi_in(t) - f_Fsi_gen(t,Hni,Hsi))/Asi)*3600;

%dVni/dt
f3 = @(t, Hni) f_Fni_spill*3600;
%dVsi/dt
f4 = @(t, Hni, Hsi) f_Fsi_spill*3600;

% Eulers Method

dt = 1; % (hr)
T_end = 8760;

% Equation for power demand
f_Pni_demand = @(t) 4065 + 1.4e6 * normpdf(t, mean, sd); % power demand (MW) in North Island?
f_Psi_demand = @(t) 1940; % power demand (MW) in South Isaland?

Pni_hydro = @(t,Hni) (0.9* (345 + 73* sin((2*pi*(t-3624))/8760)) *p*g*(Hni-Hgen_ni))/10e6;

%Equation for power demand derivative
%dPni_demand = @(t)1.4e6 * (-2*(t-mean))/(2*sd^2) * normpdf(t, mean, sd);
%dPsi_demand = @(t) 0;
% Equation for flow inlet derivative
%dFni_in = (73*2*pi/8760)*cos((2*pi*(t-3624))/8760);
%dFsi_in = (-183*2*pi/8760)*cos((2*pi*(t-2320)/8760));

%Set Hni, Hsi, Vni, Vsi arrays and initialise it
Hni = zeros(1,T_end) ;
Hsi = zeros(1,T_end);
Vni = zeros(1,T_end);
Vsi = zeros(1,T_end);

Hni(1)= Hni0;
Hsi(1) = Hsi0;

% Set time, power demand, flow inlet arrays
TArray = 1 : 1 : 8760;
Pni_demand = zeros(1,T_end);
Psi_demand = zeros(1,T_end);
Fni_in = zeros(1,T_end);
Fsi_in = zeros(1,T_end);
Fni_gen = zeros(1,T_end);

% Power hydro generation of North Island
Pni_hydroArray = zeros(1, T_end);
Eni_hydro = 0;  % amount of hydro generation in NI


for i = 1: T_end -1
   
    % differential eq for height of north island
    if Hni(i)> Hwni   % if Height of lake greater than weir
         Hni(i+1) = Hni(i) + dt*f1(i,Hni(i));
    else    % Fspill = 0 fo f1 = 0
        Hni(i+1) = Hni(i) + dt*0;
    end
   
    if Hsi(i) > Hwsi % got spill
       
        Hsi(i+1) = Hsi(i) + dt*f2a(i,Hni(i),Hsi(i));  %f2a (11/10)(t,Hni,Hsi) for time being Hni has no effect
    else   % no spill, f2a omits F_spill
        Hsi(i+1) = Hsi(i) + dt*f2b(i,Hni(i),Hsi(i));
       
    end
   
   
    %Arrays for poewr demand and low array using loop
    %Calculate power demand for South Island and North Island rate
   
    Pni_demand(i) = f_Pni_demand(i);
    Psi_demand(i) = f_Psi_demand(i);
    Fni_in(i) = f_Fni_in(i);
    Fsi_in(i) = f_Fsi_in(i);
   
    % add on last value for power demand and flow inlet
    if i == T_end - 1
        Pni_demand(i+1) = f_Pni_demand(i+1);
        Psi_demand(i+1) = f_Psi_demand(i+1);
        Fni_in(i+1) = f_Fni_in(i+1);
        Fsi_in(i+1) = f_Fsi_in(i+1);
    end
   
   
    %Fsi_in(i) = f_Fsi_in(i);
    %Set generating flow for North Island, initially same as Flow in
    Fni_gen(i)= Fni_in(i);
   
   
    % rate of energy generation hydroelectric generation in NI
    %Pni_hydro = (0.9*Fni_gen(i)*p*g*(Hni(i) - Hgen_ni))/10e6 ; % divide by 10e6 : W to MW
    Pni_hydroArray(i) = Pni_hydro(i, Hni(i));
    % amount of energy generation
    Eni_hydro = Eni_hydro + Pni_hydroArray(i);
 
   
end
Eni_hydro;

geothermal_gen = 1525; % MW geothermal_generation
total_wind_capacity = 4130; % MW N:80% S:20%
N_wind_capacity = 3304; % guess value for wind north capacity
S_wind_capacity = 826; % guess value for wind south capacity
wind_gen_N = N_wind_capacity* 0.41; % north island wind generation
Total = Pni_hydroArray + geothermal_gen + wind_gen_N;

figure(6)
plot(Total, TArray);
ylabel('total of hydro, geothermal and wind genration')

% Arrays for power demand and flow inlet in one go (Array)
%Pni_demand =  f_Pni_demand(TArray);
%for n = 1 : T_end %    generate Power demand for SI with loop because value
%                         is constane Array cant work
 %Psi_demand(n) = f_Psi_demand(n);
%end
%Fni_in = f_Fni_in(TArray);
%Fsi_in = f_Fsi_in(TArray);


% Plotting

figure(1) % Power demand
plot(TArray, Pni_demand);
hold on
plot(TArray, Psi_demand);
hold off
legend("Power demand NI (MW)","Poewr Demand SI (MW)");


figure(2) % Flow rate (m^3/s)
plot(TArray, Fni_in);
hold on
plot(TArray, Fsi_in);
hold off
legend('Flow rate NI (m^3/s)','SI (m^3/s)');

figure(3) % 4 variables
plot(TArray, Hni);
%hold on
%plot(TArray, Hsi);   %%% !!! Something wrong with Hsi goes to infinity
%hold off

legend('height NI (m)');

figure(4)
plot(TArray,Pni_hydroArray); % hydro_electric geneation
ylabel('Power Hydro North Island(MW)');


% caculate the wind capacity factor
CF = @(t) 0.41 + 0.12*sin((2*pi*(t-5560))/8760);

for i = 1: T_end
    CF_Array(i) = CF(i);
end

figure(5)
plot(TArray, CF_Array)
ylabel('capacity factor for wind')