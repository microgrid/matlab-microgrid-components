function [ av ] = DailyAverage( input )
% DailyAverage() computes the daily average of yearly input data
% INPUT: 'input' should be a row vector with hourly values throughout 1 year 
% i.e. a vector of length 24 hours x 365 days = 8760.
% OUTPUT: a row vector of length 24 giving the average daily values.

if mod(length(input), 24) ~= 0
    error('ERROR: this yearly input is not hourly i.e. not divisible by 24.')
end
nr_days = length(input) / 24;
average = zeros(1,24);                       % vector for average daily output per hour

for hour = 1:24                                     % iterate over all times 1:00, 2:00 etc.
hours_i = hour : 24 : (nr_days - 1) * 24 + hour;    % range to pick the i-th hour of each day throughout the yearly data, i.e. 1:00 of 1 January, 1:00 of 2 January etc.
    for k = hours_i                                 % sum all values of the i-th hour of the day over all days
        average(hour) = average(hour) + input(k);
    end
    average(hour) = average(hour) / nr_days;        % divide this sum by the number of days to compute average
end

av = average;

end

