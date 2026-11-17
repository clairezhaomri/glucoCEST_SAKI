%% baseline spectrum
% load('spectrum_ctx.mat')
load('spectrum_medulla.mat')


lower_bound = 0;  % Adjusted lower bound as per your request
upper_bound = 2.88;

figure,
% Asymmetry plot
asym_x = (xaxisCorrected(ZeroIdx+1:end))./400;
asym_y = asym_Base(ZeroIdx+1:end);

% Find the indices of asym_x that correspond to the range
indices = asym_x >= lower_bound & asym_x <= upper_bound;

% Select the part of asym_x and asym_y within the range
asym_x_in_range = asym_x(indices);
asym_y_in_range = asym_y(indices);

plot(asym_x_in_range, asym_y_in_range, 'r--', 'LineWidth', 3); % Solid line for the specified range
ylim([0 0.15]);
ylabel('MTR_{asym}', 'FontSize', 20);
set(gca, 'FontSize', 18, 'TickLength', [0.025 0.025]); % Set left y-axis color to red and ticks larger
yticks(0:0.05:0.15);
xlabel('Offset (ppm)')
% title('ctx')
title('medulla')


% Connect the ends of the red spectrum to zero with dashed lines
hold on;

%% infusion spectrum
% load('spectrum_ctx.mat')
load('spectrum_medulla.mat')


asym_x = (xaxisCorrected(ZeroIdx+1:end))./400;
asym_y = asym_Base(ZeroIdx+1:end);

% Find the indices of asym_x that correspond to the range
indices = asym_x >= lower_bound & asym_x <= upper_bound;

% Select the part of asym_x and asym_y within the range
asym_x_in_range = asym_x(indices);
asym_y_in_range = asym_y(indices);

plot(asym_x_in_range, asym_y_in_range, 'r-', 'LineWidth', 3); % Solid line for the specified range
ylim([0 0.15]);

set(gca, 'GridColor', 'k', 'GridAlpha', 0.2); % Set grid color to black and partially transparent
grid on;

set(gca,'XDir','reverse');% 
%%
saveas(gca,'glucoCEST_pre_post_medulla.fig')
