%% Shift the standard load profile to the left such that it overlaps more with the irradiation. 
% This can be used to compare results and see how the system behaves.

% import the standard load profile
year = 100;
path_to_dataBase = '/Users/jeemijn/Desktop/NTNU_microgrids/matlab/matlab-microgrid-components/dataBase/';
filename = ([path_to_dataBase, 'LoadCurve_normalized_single_3percent_',num2str(year),'.mat']);      % Average hourly global radiation (beam + diffuse) incident on the PV array [kW/m2]. Due to the simulation step [1h], this is also [kWh/m2]
Load = importdata(filename); 

% plot using average function


% shift 3x to left



figure(10)
plot(Load)

figure(11)
plot(Load_av)
hold on
plot(P_pv_av)
hold off
