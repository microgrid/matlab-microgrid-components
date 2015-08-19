clear all

markertype(1:8)= cellstr('');         % create array of empty strings
markertype(1:5)= cellstr('filled');
markertype(6:8)= cellstr('o');

figure(1)
% x_val = 1:8;
% y_val = x_val;
% scatter(x_val, y_val,markertype)

x_val = [3 cellstr('') 5];
y_val = x_val;

scatter(1:2,1:2, 'filled')
hold on
scatter(x_val,y_val, 'o')







