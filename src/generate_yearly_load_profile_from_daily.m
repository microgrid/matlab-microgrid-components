%% Given a daily load profile (per hour) this script generates a yearly load profile with simply all days exactly the same

clear all

% paste here the hourly profile for 1 day from 0:00-1:00 to 23:00-00:00 in
% [kW]
daily = [0	0	0	0	0	0	0   0.01	0	0	0	0	0	0	0	0	0.02    0.02	0.012	0.032	0.038	0.018	0.006	0];
yearly = repmat(daily, 1, 365);         % copy matrix 'daily' 1 time in row direction and 365 times in column direction

if length(yearly) ~= 8760
   error('Something went wrong: the total size is not 365 days x 24 hours.') 
end

path_to_dataBase = '/Users/jeemijn/Desktop/NTNU_microgrids/matlab/matlab-microgrid-components/dataBase/';
save([path_to_dataBase,'LoadCurve_bangladesh.mat'],'yearly')
