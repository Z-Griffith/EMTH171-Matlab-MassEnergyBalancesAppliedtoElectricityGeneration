% Task 3-b
close all 
clear
clc

% North Island constant factors
A_NI = 620e6; % Surface area in m^2
h_min_NI = 355.85; % Minimum height in m 
h_max_NI = 357.25; % Maximum height in m
h_gen_NI = 80; % Generator height in m 

% South constant factors
A_SI = 350e6; % Surface area in m^2
h_min_SI = 402; % Minimum height in m 
h_max_SI = 410; % Maximum height in m
h_gen_SI = 0; % Generator height in m 

% Other constants
rho = 998; % Density in kg/m^3
g = 9.81; % Earth gravity in m/s^2
K = 1.55; % Constant of spill water in m^1.5/s
L = 300;  % Length of spillway weir in m
dt = 1; % Time step for Euler method in h 
N = 8760 / dt + 1; % Number of step we should solve 

% Define Matrices
h_NI = zeros(N,1); % Defining north island lake's height matrice
h_SI = zeros(N,1); % Defining south island lake's height matrice
V_spill_NI = zeros(N,1); % Defining north island spillway flow matrice
V_spill_SI = zeros(N,1); % Defining south island spillway flow matrice


% Set initial values
h_NI(1) = (h_min_NI + h_max_NI)/2; % Set initial value of north island height
h_SI(1) = (h_min_SI + h_max_SI)/2; % Set initial value of south island height

% Assuming equal wind capacity 
wind_gen_NI = 4130/2; % Assuming wind power for north island 
wind_gen_SI = 4130/2; % Assuming wind power for north island

% Euler method to solve for 1 year (8760 hour)
for i=2:N
    t = (i-1) * dt; % calculating time based on the step
    P_demand_NI = 4065 + 1.4e6*normpdf(t,5000,1000); % Power demand of north island in MW
    P_demand_SI = 1940; % Power demand of south island in MW
    
    F_in_NI = 345 + 73*sin(2*pi*(t-3624)/8760); % Inlet flow of north island in m^3/s
    F_in_SI = 593 - 183*sin(2*pi*(t-2320)/8760); % Inlet flow of south island in m^3/s
    
    F_gen_NI = F_in_NI; % Generator flow of north island in m^3/s
    
    P_hydro_NI = 0.9*F_gen_NI*rho*g*(h_NI(i-1) - h_gen_NI)/1e6; % Hydro-electric generation of north island in MW
    
    CF = 0.41 + 0.12*sin(2*pi*(t-5660)/8760); % Calculating capacity factor
    if t>4344 && t<=5088 % Reducing wind power capacity to 50 percent in July
        CF = 0.5*CF;
    end
    P_wind_NI = wind_gen_NI * CF; % wind power of north island in MW
    P_wind_SI = wind_gen_SI * CF; % wind power of south island in MW
    
    P_geo_NI = 1605; % Geothermal power of north island in MW
    P_HVDC = P_demand_NI - P_geo_NI - P_wind_NI - P_hydro_NI; % P_HVDC in MW
    
    P_hydro_SI = P_demand_SI + P_HVDC - P_wind_SI; % Hydro-electric power for south island in MW
    
    F_gen_SI = P_hydro_SI*1e6/0.9/rho/g/(h_SI(i-1)-h_gen_SI); % Generator flow rate needed for south island in m^3/s
    
    % Checking the spill water
    if h_NI(i-1) > h_max_NI
        F_spill_NI = K*L*(h_NI(i-1) - h_max_NI)^1.5; % Spill water flow rate for north island in m^3/s
    else
        F_spill_NI = 0;
    end
    if h_SI(i-1) > h_max_SI
        F_spill_SI = K*L*(h_SI(i-1) - h_max_SI)^1.5; % Spill water flow rate for south island in m^3/s
    else
        F_spill_SI = 0;
    end
    
    % Calculating derivatives
    dh_NI_dt = (F_in_NI - F_gen_NI - F_spill_NI)*3600/A_NI; % Derivative of north island lake height to time in m/h
    dh_SI_dt = (F_in_SI - F_gen_SI - F_spill_SI)*3600/A_SI; % Derivative of south island lake height to time in m/h
    dVspill_NI_dt = F_spill_NI*3600; % Derivative of north island spill water volume to time in m^3/h
    dVspill_SI_dt = F_spill_SI*3600; % Derivative of south island spill water volume to time in m^3/h
    
    h_NI(i) = h_NI(i-1) + dt*dh_NI_dt; % Calculating north island lake level in m
    h_SI(i) = h_SI(i-1) + dt*dh_SI_dt; % Calculating south island lake level in m
    V_spill_NI(i) = V_spill_NI(i-1) + dt*dVspill_NI_dt; % Calculating north island spill water volume in m^3
    V_spill_SI(i) = V_spill_SI(i-1) + dt*dVspill_SI_dt; % Calculating south island spill water volume in m^3
end

% Plot figures
time = linspace(0,8760,N);

figure(1)
plot(time , h_SI)
xlabel('Time (h)')
ylabel('South Lake Level (m)')
grid on

figure(2)
plot(time , h_NI)
xlabel('Time (h)')
ylabel('North Lake Level (m)')
grid on

figure(3)
plot(time , V_spill_SI)
xlabel('Time (h)')
ylabel('South Spille Volume (m^3)')
grid on

figure(4)
plot(time , V_spill_NI)
xlabel('Time (h)')
ylabel('North Spille Volume (m^3)')
grid on

% Calculating 
fprintf('Minimum lake level for north island is: %f \n' , min(h_NI));
fprintf('Minimum lake level for south island is: %f \n\n' , min(h_SI));

fprintf('Spilled water for north island is: %f \n' , V_spill_NI(end));
fprintf('Spilled water for south island is: %f \n' , V_spill_SI(end));















