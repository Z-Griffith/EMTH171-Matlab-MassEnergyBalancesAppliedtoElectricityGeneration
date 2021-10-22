% EMTH171 // CaseStudy2
% Name: Menghao Zhan   //   Jiyao Zhu
clear 
clc
close all
%====================--Constant Values--======================
g = 9.81;             % Coeficience of  Gravity(m/s^2) 
tarray = 1:1:8760;    % One year shown by array (hours)
Stand_d = 1000;       % Standard deviation (h)
Mean = 5000; % Mean of distribution with North Island Power demand (h)
density = 998;        % (kg/m^3)
K = 1.55;             % (m^0.5*s^-1)
L = 300;              % length of weir (m)
%-------------------Assume NI&SI wind capacity-------------------- 
%===================--Task1--==========================
wind_geCA_NI = 4130/2;%(MW)
wind_geCA_SI = 4130/2;%(MW)
%===================--Task1--==========================
%===================--Task3(abcd)--==========================
% wind_geCA_NI = Wind_Capacity; % The value get from task2 after improving
% wind_geCA_NI = Wind_Capacity; % The value get from task2 after improving
%===================--Task3(abcd)--==========================
%===================--Known Values--==========================
Area_SI = 350e6;      % Suface area of south island(m^2)
Area_NI = 620e6;      % Suface area of north island(m^2)
Minh_SI = 402;        % Minimum height of south island(m)
Minh_NI = 355.85;     % Minimum height of north island(m)
Maxh_SI = 410;        % Maximum height of south island(m)
Maxh_NI = 357.25;     % Maximum height of north island(m)
h_ge_SI = 0;          % Height of generator of south island(m)
h_ge_NI = 80;         % Height of generator of north island(m)
Av_flow_SI = 593;     % Average flow of south island(m^3/s)
Av_flow_NI = 345;     % Average flow of north island(m^3/s) 
Maxge_SI = 3590;      % Maximum generation(MW)
Maxge_NI = 1870;      % Maximum generation(MW)
E_NI_ge = 1525;       % Energy of geothermal capacity(MW)
%===================--Initial Values--========================
hSI(1) = (Maxh_SI+Minh_SI)/2 % Initial height of SI lake (m)
hNI(1) = (Maxh_NI+Minh_NI)/2 % Initial height of NI lake (m)
vNI(1) = 0;                  % Spillage of NI lake (m^3)
vSI(1) = 0;                  % Spillage of SI lake (m^3)
P_SI_demand = 1940;          % SI power demand in MW(Eq14,15)
%=====================--For-loop--==========================
for n=2:tarray(end)
    t = n;
    % 1 ----------NI power demand in MW(Eq14,15)----------
    P_NI_demand = 4065 + 1.4e6 * normpdf(t, Mean, Stand_d);
    % 2 ----------Inlet flow rate into NI&SI lakes(Eq5,6)----------
    F_NI_inlet =  345 + 73* sin((2*pi*(t-3624))/8760);
    F_SI_inlet =  593 - 183* sin((2*pi*(t-2320))/8760);
%=====================--Task3-c--==========================
%     F_SI_inlet =  0.9*(593 - 183* sin((2*pi*(t-2320))/8760));
%=====================--Task3-c--==========================
    F_NI_ge = F_NI_inlet;
%=====================--Task3-a--==========================    
%     if t < 8760/4
%         F_NI_ge  = 0.9*F_NI_inlet; %Reducing generation flow 
%     elseif t < 8760*3/4
%         F_NI_ge = 1.04*F_NI_inlet; 
%     else
%         F_NI_ge = 0.9*F_NI_inlet;
%     end
%=====================--Task3-a--==========================
%=====================--Task3-d--==========================
    if t < 8760/6 || 8760*(5/6)<t
        % Reducing generation flow to 0.8 accln water
        F_NI_ge  = 0.8*F_NI_inlet; 
        F_SI_ge  = 0.8*F_SI_inlet; 
        CF = CF / 2;
    else
        % The generation flow dose not change /the wind power capacity changed
        F_NI_ge  = 1*F_NI_inlet;
        F_SI_ge  = 1*F_SI_inlet; 
        CF = CF * 1.3;
       % F_SI_ge = 1.2*F_SI_inlet;
    end
%=====================--Task3-d--==========================   
%=====================--Task3-b--==========================     
%     if t<=5088 && t>4344 % Reducing wind power capacity to 50%
%         CF = CF/2;
%     end
%=====================--Task3-b--==========================     
    % 4  --------The hydro-electric generation in NI(Eq7)----------
    P_NI_hydro = 0.9*F_NI_ge*density*g*(hNI(n-1)-h_ge_NI)/(1e6);
    % 6 ---- Power of wind-------
    P_NI_wind = wind_geCA_NI * CF;
    P_SI_wind = wind_geCA_SI * CF;
    %----Electrical power balances for each island -------
    P_HVDC =  P_NI_demand - E_NI_ge - P_NI_wind - P_NI_hydro;
    % 4-1--------The hydro-electric generation in SI------------
    P_SI_hydro = P_SI_demand + P_HVDC - P_SI_wind;
    % 3 ----------The generating flow for NI&SI(Eq8)----------
    F_SI_ge =  (P_SI_hydro*1e6)/((0.9*density*g)*(hSI(n-1)-h_ge_SI));
    % 10 ---- Spillway flow for each lake(Eq9) -------
    if hNI(n-1) > Maxh_NI
        F_NI_spill = K*L*(hNI(n-1)-Maxh_NI)^1.5;
        Spill_NI_dvdt = F_NI_spill*3600;% Derivative of NI spill water volume (m^3/h)
    else 
        F_NI_spill = 0;
        Spill_NI_dvdt = 0;
    end
    if hSI(n-1) > Maxh_SI
        F_SI_spill = K*L*(hSI(n-1)-Maxh_SI)^1.5;
        Spill_SI_dvdt = F_SI_spill*3600;% Derivative of SI spill water volume (m^3/h)
    else 
        F_SI_spill = 0;
        Spill_SI_dvdt = 0;
    end
    %-----------Main part----------
    dhNIt = (F_NI_inlet - F_NI_ge - F_NI_spill) * 3600/Area_NI;  % (m^3/h)
    dhSIt = (F_SI_inlet - F_SI_ge - F_SI_spill) * 3600/Area_SI;  % (m^3/h)
    %-----------Euler's method-----
    hNI(n) = hNI(n-1) + dhNIt; %Current lake level of North island(m)
    hSI(n) = hSI(n-1) + dhSIt; %Current lake level of South island(m)
    vNI(n) = vNI(n-1) + (Spill_NI_dvdt);%Current lake volume of NI lake (m^3)
    vSI(n) = vSI(n-1) + (Spill_SI_dvdt);%Current lake volume of SI lake (m^3)
    % Calculate money if the minimum criteria condition passed
%=====================--Task2---==========================
%----Find the inmproving wind genorated capacity of north island-------
    if check == 0;
        Total_cost = sum(cost_NI) + sum(cost_SI);
            if Total_cost < A_money
                A_money = Total_cost;
                Wind_Capacity = wind_geCA_NI;
                New_H_NI = hNI;
                New_H_SI = hSI;
                New_spillNI = vNI;
                New_spillSI = vSI;
            end
    end
%=====================--Task2---==========================  
end  
%=====================--Task2/3---==========================
% Calculation of total cost from North & South islands.
cost_NI = 90*density*(vNI(tarray(end)))*g*(hNI(n)-h_ge_NI)/(3600*1e6); 
cost_SI = 90*density*(vSI(tarray(end)))*g*(hSI(n)-h_ge_SI)/(3600*1e6);
total = cost_NI + cost_SI
%=====================--Task2/3---==========================
%==============================ploting================================
% figure(1)
% plot(tarray,hNI)
% hold on 
% plot(tarray,hSI)
% hold off
% xlabel('Time (h)')
% ylabel('Lake Level (m)')
% legend("North Island","South Island")
%==============================ploting================================
figure(1)
plot(tarray , hNI)
xlabel('Time(h)')
ylabel('North lake-level(m)')
figure(2)
plot(tarray , hSI)
xlabel('Time(h)')
ylabel('South lake-level(m)')
figure(3)
plot(tarray , vSI)
xlabel('Time(h)')
ylabel('South spill-volume(m^3)')
figure(4)
plot(tarray , vNI)
xlabel('Time(h)')
ylabel('North spill-volume(m^3)')
%=====================--Task2/b--==========================
E_NI_ge = 1605; % When calculating the cost the power from north island geothermal is original one
Wind_cost = (Wind_Capacity-Extra_wind)*Cost_wt*trans;
Geo_cost =  (E_NI_ge-Extra_geo)*Cost_Ngc*trans; 
Cost_extra_geo = (Wind_cost + Geo_cost)/1e6
%=====================--Task2/b--==========================



