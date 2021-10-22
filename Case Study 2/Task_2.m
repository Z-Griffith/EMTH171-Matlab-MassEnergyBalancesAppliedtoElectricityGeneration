% Task 2
close all
clear 
clc

% Section 2-A
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

j  = 1; 
min_Money = 1e20; % Assuming money we want to minimize
for wind_gen_NI = 4130/2:5:5000/2 % Starting with wind generation of north island in MW
    % calculating wind capacity of SI
    wind_gen_SI = wind_gen_NI; % Wind generation of south island in MW

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

        money_NI(i) = 0.9*rho*(V_spill_NI(i)-V_spill_NI(i-1))*g*(h_NI(i)-h_gen_NI)/3600/1e6*100; % Money value of spilled water in north island in dollar
        money_SI(i) = 0.9*rho*(V_spill_SI(i)-V_spill_SI(i-1))*g*(h_SI(i)-h_gen_SI)/3600/1e6*100; % Money value of spilled water in south island in dollar
    end

    % Check if the level is below the minimum criteria 
    jj = 0;
    for i=2:N
        if (h_NI(i)<h_min_NI) || (h_SI(i)<h_min_SI)
            jj = 1;
        end
    end
    
    % Calculate money if the minimum criteria condition passed
    if jj == 0
        Money = sum(money_NI + money_SI); % Total money in dollar
            if Money < min_Money
                min_Money = Money;
                wind_capacity = wind_gen_NI;
                hh_NI = h_NI;
                hh_SI = h_SI;
                VV_spill_NI = V_spill_NI;
                VV_spill_SI = V_spill_SI;
            end
        j = j + 1;
    end
    
end


% Showing results of 2-A
fprintf('The minimum dollar value for splitting water is: %2.4f million dollar\n' , min_Money/1e6);
fprintf('The capacity need for wind turbines are: %d MW\n' , 2*wind_capacity);
fprintf('Wind capacity of each island: %d MW \n\n\n' , wind_capacity)

% Plotting figures
time = linspace(0,8760,N);

figure(1)
plot(time , hh_SI)
xlabel('Time (h)')
ylabel('South Lake Level (m)')
grid on

figure(2)
plot(time , hh_NI)
xlabel('Time (h)')
ylabel('North Lake Level (m)')
grid on

figure(3)
plot(time , VV_spill_SI)
xlabel('Time (h)')
ylabel('South Spille Volume (m^3)')
grid on

figure(4)
plot(time , VV_spill_NI)
xlabel('Time (h)')
ylabel('North Spille Volume (m^3)')
grid on


% Section 2-B
geo_cost = 4500; % Cost of geothermal power in dollar per installed kW 
wind_cost = 3000; % Cost of wind power in dollar per installed kW

geo_pre = 960; % Installation power of geothermal power in MW
geo_need = P_geo_NI; % Amount of geothermal power needed in MW

wind_pre = 580; % Installation power of wind power in MW
wind_need = wind_capacity; % Amount of wind power needed in MW

Cost = (geo_need - geo_pre)*1000*geo_cost + ...
       (wind_need - wind_pre)*1000*geo_pre;
fprintf('The cost of the extra generating capacity is: %2.4f million dollar. \n' , Cost/1e6)

