%% Shift the standard load profile to the left such that it overlaps more with the irradiation. 
% This can be used to compare results and see how the system behaves.

clear all

% import the standard load profile
year = 100;
path_to_dataBase = '/Users/jeemijn/Desktop/NTNU_microgrids/matlab/matlab-microgrid-components/dataBase/';
filename = ([path_to_dataBase, 'LoadCurve_normalized_single_3percent_',num2str(year),'.mat']);      % Average hourly global radiation (beam + diffuse) incident on the PV array [kW/m2]. Due to the simulation step [1h], this is also [kWh/m2]
Load = importdata(filename); 

% import the solar irradiation 
% (P_pv would be better but is the same shape up to temperature differences) 
% (only rough shape now since I don't want to copy >5 lines of code from opt_comp_sizes script)
irr = importdata([path_to_dataBase, 'solar_data_Phuntsholing_baseline.mat']);                       % Use \ for Windows and / for Mac and Linux
irr_scaled = irr .* 100;

% shift the load profile
Load_3_hrs_earlier = Load;                                  % vector to make a new Load profile that is shifted to the left by 3 hours
Load_6_hrs_earlier = Load;

first_3_hrs = Load(1:3);                                    % save first 3 hours of the year
Load_3_hrs_earlier(1:3) = [];                               % delete first 3 hours of the year from Load profile
Load_3_hrs_earlier = [Load_3_hrs_earlier first_3_hrs];      % add the first 3 hours at the end of the year in order to have a full year again (with same total energy) 

first_6_hrs = Load(1:6);                                    % repeat procedure for a shift of 6 hours to the left
Load_6_hrs_earlier(1:6) = [];                                  
Load_6_hrs_earlier = [Load_6_hrs_earlier first_6_hrs];    

% save shifted load profiles
save([path_to_dataBase,'LoadCurve_3_hrs_earlier.mat'],'Load_3_hrs_earlier')
save([path_to_dataBase,'LoadCurve_6_hrs_earlier.mat'],'Load_6_hrs_earlier')

% plot irradiation and average daily load profiles to check what happened
irr_scaled_av = DailyAverage(irr_scaled);
Load_av = DailyAverage(Load);
Load_av_3_hrs_earlier = DailyAverage(Load_3_hrs_earlier);
Load_av_6_hrs_earlier = DailyAverage(Load_6_hrs_earlier);

figure(10)
plot(irr_scaled_av,'r')
hold on
plot(Load_av,'b')
hold on
plot(Load_av_3_hrs_earlier,'b--')
hold on
plot(Load_av_6_hrs_earlier,'b:')
hold off
set(gca,'FontSize',12,'FontName','Times New Roman','fontWeight','bold')
xlabel('Time (hour)')
ylabel('kW')
title('An average day')
legend('irradiation x 100', 'Load', 'Load 3 hrs to left', 'Load 6 hrs to left','Location','northwest')
