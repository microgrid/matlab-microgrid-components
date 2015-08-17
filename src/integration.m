%% test for integration of areas under/between graphs

clear all

%% plotting test graphs
time = 0:24;
load_test = 1 + sin(time * (3 /(2 * pi)));
irr_test = 1 - cos(time / pi);
free = min(load_test, irr_test);            % energy for free, i.e. directly from PV without battery intervenience, is the area under this graph. 
                                            % N.B. Assumption: both load and irr are positive functions
figure(1)
plot(time, load_test)
hold on 
plot(time, irr_test)
hold on
plot(time, free)
hold off
legend('Load','Irradiance','Free area')
xlabel('Time over the day [hour]')
ylabel('Energy [kWh]')

%% integration

free_area = trapz(time, free);                                      % discrete integration using trapeziums. Might not be the best solution for non-linear data. See http://se.mathworks.com/help/matlab/math/integration-of-numeric-data.html
area_to_batt = trapz(time, irr_test) - free_area;                   % by definition this should be positive
area_load_needed_from_batt = trapz(time, load_test) - free_area;    % by definition this should be positive

unmet_load = area_load_needed_from_batt - area_to_batt;
unmet_load_perc = unmet_load / trapz(time, load_test) * 100         % equal to Loss of Load Probability. But rough estimate since SoC at end of the day influences next day!