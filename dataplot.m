clear all;
close all;

file_path = 'formation_error.csv';

% Read the file as text
file_content = fileread(file_path);

% Split into lines
lines = strsplit(file_content, '\n');

% Extract CH1–CH7 rows
ch1_line = strsplit(lines{2}, ',');  
ch2_line = strsplit(lines{3}, ',');  
ch3_line = strsplit(lines{4}, ',');  
ch4_line = strsplit(lines{5}, ',');  
ch5_line = strsplit(lines{6}, ',');  
ch6_line = strsplit(lines{7}, ',');  
ch7_line = strsplit(lines{8}, ',');  

% Convert from cell array of strings to numeric (skip the first entry 'CHx:')
ch1_data = str2double(ch1_line(2:end));
ch2_data = str2double(ch2_line(2:end));
ch3_data = str2double(ch3_line(2:end));
ch4_data = str2double(ch4_line(2:end));
ch5_data = str2double(ch5_line(2:end));
ch6_data = str2double(ch6_line(2:end));
ch7_data = str2double(ch7_line(2:end));

% Plot main figure
font_size = 20;
font_family = 'Times New Roman';

figure('Position', [100, 100, 1100, 600]);
plot(ch2_data, '--','Color', [0.2 0.5470 0.7410], 'LineWidth', 2.5, 'DisplayName', 'Reference trajectory');
hold on;
plot(ch1_data, '-m','LineWidth', 3.5, 'DisplayName', 'Motor 1');
plot(ch4_data, '-r','LineWidth', 2.8, 'DisplayName', 'Motor 2');
plot(ch5_data, '-g', 'LineWidth', 2.8, 'DisplayName', 'Motor 3');
plot(ch3_data, '--k', 'LineWidth', 2.8, 'DisplayName', 'Deviation 1');
plot(ch6_data, '--b','LineWidth', 2.8, 'DisplayName', 'Deviation 2');
plot(ch7_data, '--', 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 1.8, 'DisplayName', 'Deviation 3');

xlim([0 250]);
ylim([10 150]);
xlabel('Time step (k)');
ylabel('Tracking performance');
legend('Location', 'north', 'Orientation', 'horizontal', 'NumColumns', 4);
set(gca, 'FontSize', font_size, 'FontName', font_family);

% ---------------------------
% Inset 1
% ---------------------------
axes('Position', [0.22, 0.18, 0.25, 0.20]);  % [x, y, width, height]
box on; hold on;
plot(ch2_data, '--','Color', [0.2 0.5470 0.7410], 'LineWidth', 1.5);
plot(ch1_data, '-m','LineWidth', 1.5);
plot(ch4_data, '-r','LineWidth', 1.5);
plot(ch6_data, '--b','LineWidth', 1.5);
plot(ch5_data, '-g','LineWidth', 1.5);
plot(ch7_data, '--','Color', [0.4940 0.1840 0.5560],'LineWidth', 1.5);

xlim([78 100]);   % zoom range 1
ylim([65 92]);
set(gca, 'FontSize', 10, 'FontName', font_family);

% ---------------------------
% Inset 2
% ---------------------------
axes('Position', [0.60, 0.18, 0.25, 0.20]);  % [x, y, width, height]
box on; hold on;
plot(ch2_data, '--','Color', [0.2 0.5470 0.7410], 'LineWidth', 1.5);
plot(ch1_data, '-m','LineWidth', 1.5);
plot(ch3_data, '--k', 'LineWidth', 1.5);
plot(ch4_data, '-r','LineWidth', 1.5);
plot(ch5_data, '-g','LineWidth', 1.5);
xlim([180 210]);   % zoom range 2
ylim([112 130]);
set(gca, 'FontSize', 10, 'FontName', font_family);
