function WZ_Replot_reverseXandTwoAxis(path)
clear all;clc;
load('spectrum_ctx.mat');

% load('spectrum_medulla.mat');

lower_bound = 0;  % Adjusted lower bound as per your request
upper_bound = 2.88;

figure,
% CEST spectrum
yyaxis right;
plot(xaxisAcqCorrect./400, ROIMeanAcqBaseCorrected, 'k*', xaxisCorrected./400, ROIMeanBase, 'k-', 'LineWidth', 2)
xlim([-5 5]);
xticks(-5:1:5);
ylim([0 1]);
yticks(0:0.2:1);
xlabel('offset (ppm)')
ylabel('S/S_0', 'FontSize', 18, 'FontWeight', 'bold');
set(gca, 'YColor', 'k', 'FontSize', 16, 'TickLength', [0.025 0.025], 'FontWeight', 'bold'); % Set right y-axis color to black and ticks larger

% Asymmetry plot
asym_x = (xaxisCorrected(ZeroIdx+1:end))./400;
asym_y = asym_Base(ZeroIdx+1:end);

% Find the indices of asym_x that correspond to the range
indices = asym_x >= lower_bound & asym_x <= upper_bound;

% Select the part of asym_x and asym_y within the range
asym_x_in_range = asym_x(indices);
asym_y_in_range = asym_y(indices);

yyaxis left;
plot(asym_x_in_range, asym_y_in_range, 'r-', 'LineWidth', 2); % Solid line for the specified range
ylim([0 0.2]);
ylabel('MTR_{asym}', 'FontSize', 18, 'FontWeight', 'bold');
set(gca, 'YColor', 'r', 'FontSize', 16, 'TickLength', [0.025 0.025], 'FontWeight', 'bold'); % Set left y-axis color to red and ticks larger

% Connect the ends of the red spectrum to zero with dashed lines
hold on;
plot([lower_bound, asym_x_in_range(1)], [0, asym_y_in_range(1)], 'r--', 'LineWidth', 2);
plot([asym_x_in_range(end), upper_bound], [asym_y_in_range(end), 0], 'r--', 'LineWidth', 2);
hold off;

% Set the grid color and enable the grid
set(gca, 'GridColor', 'k', 'GridAlpha', 0.5); % Set grid color to black and partially transparent
grid on;

set(gca,'XDir','reverse');% 
% saveas(gca,'selectfreq_medulla.fig');
% saveas(gca,'selectfreq_medulla.fig');

% close all;
end

