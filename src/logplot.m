%% Plotting the investment cost (IC) vs. loss of load probability (LLP) curves in a logplot to check that they have the same shape for all percentages. 
% And to see whether the relation is exponential.

% Conclusion: For all percentages the plot indeed has the same shape.
% Moreover, since the loglogplot gives as good as straight lines, the relation
% between IC and LLP is (very close to) exponential.

%% for copying symbols
%   []
%   €
%   ''

%% clear
clear all;
close all;

%% importing data
data_small = importdata('MA_opt_phuentsholing_10_1.mat');
data_large = importdata('MA_opt_phuentsholing_100_10.mat');

% DESCRIPTION OF THE DATAFILES: 

% data_small describes 10 graphs; for supplying 1% to 10% of total load in steps of 1%
% data_large describes 10 graphs; for supplying 10% to 100% of total load in steps of 10%

% For example, data_small is sorted as:
% - in the first 6 columns (1-6) are the 6 optimal parameters for supplying 1% of total
% load. Namely the 6 parameters XXXXXXX
% - the 20 rows give the 20 optimal systems that are found for each
% fixed percentage of total load
% - in the second 6 columns (7-12) are the 6 optimal parameters for supplying 2% of total
% load, etc.


%% plotting 
% plotting from the importfile data_small
figure(1);
for i = 1:10
    column_loss_load = 1 + 6 * (i - 1);                 % Loss of Load Probability (LLP) of optimal system is given in first column of data
    column_investment = 6 + 6 * (i - 1);                % Investment Cost (IC) of optimal system is given in sixth column of data
    loss_load = data_small(:, column_loss_load);        
    investment = data_small(:, column_investment);      
    
%     blue_plot = plot(investment, loss_load, 'b-*');                  % linear plot
%     blue_plot = semilogx(investment, loss_load, 'b-*');              % plot with logarithmic x-axis
    blue_plot = loglog(investment, loss_load, 'b-*');              % plot with both axes logarithmic
    hold on;
end

% plotting from the importfile data_large
for i = 1:10
    column_loss_load = 1 + 6 * (i - 1);                 % Loss of Load Probability (LLP) of optimal system is given in first column of data
    column_investment = 6 + 6 * (i - 1);                % Investment Cost (IC) of optimal system is given in sixth column of data
    loss_load = data_large(:, column_loss_load);        
    investment = data_large(:, column_investment);      
    
%     red_plot = plot(investment, loss_load, 'r-*');                  % linear plot
%     red_plot = semilogx(investment, loss_load, 'r-*');              % plot with logarithmic x-axis
    red_plot = loglog(investment, loss_load, 'r-*');              % plot with both axes logarithmic
    hold on;
end

xlabel('Investment Cost (IC) [€]');
ylabel('Loss of Load Probability (LLP) [%]');
title('LLP versus IC for optimal system choices');
legend([blue_plot red_plot], '1-10% of total load','10-100% of total load');   % only show 1 of the blue and 1 of the red plots in legend
hold off;

% print('cost_vs_llp_linear_plot', '-dpng');                                       % save plot as .png in current directory
% print('cost_vs_llp_loglinear_plot', '-dpng');                                       
print('cost_vs_llp_loglog_plot', '-dpng');                                       