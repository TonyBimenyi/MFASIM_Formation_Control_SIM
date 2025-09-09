clear all;
close all;

file_path = 'formation_error.csv';

% Read the file as text
file_content = fileread(file_path);

% Split into lines
lines = strsplit(file_content, '\n');

% Extract CH1–CH5 rows
ch1_line = strsplit(lines{2}, ',');  % CH1: row
ch2_line = strsplit(lines{3}, ',');  % CH2: row
ch3_line = strsplit(lines{4}, ',');  % CH3: row
ch4_line = strsplit(lines{5}, ',');  % CH4: row
ch5_line = strsplit(lines{6}, ',');  % CH5: row
ch6_line = strsplit(lines{7}, ',');  % CH6: row
ch7_line = strsplit(lines{8}, ',');  % CH7: row

% Convert from cell array of strings to numeric (skip the first entry 'CHx:')
ch1_data = str2double(ch1_line(2:end));
ch2_data = str2double(ch2_line(2:end));
ch3_data = str2double(ch3_line(2:end));
ch4_data = str2double(ch4_line(2:end));
ch5_data = str2double(ch5_line(2:end));
ch6_data = str2double(ch6_line(2:end));
ch7_data = str2double(ch7_line(2:end));

% Example: Plot
font_size =20;
font_family = 'Times New Roman';

figure('Position', [100, 100, 1100, 600]);  % [x, y, width, height]
plot(ch2_data, '--','Color', [0.2 0.5470 0.7410], 'LineWidth', 2.5,  'DisplayName', 'Reference trajectory');
hold on;
plot(ch1_data, '-m','LineWidth', 3.5, 'DisplayName', 'Motor 1');
plot(ch4_data, '-r','LineWidth', 2.8, 'DisplayName', 'Motor 2');
plot(ch5_data, '-g', 'LineWidth', 2.8,'DisplayName', 'Motor 3');
plot(ch3_data, '--k', 'LineWidth', 2.8, 'DisplayName', 'Deviaton 1');


plot(ch6_data, '--b','LineWidth', 2.8, 'DisplayName', 'Deviation 2');
plot(ch7_data, '--', 'LineWidth', 2.8, 'Color', [0.4940 0.1840 0.5560],'LineWidth', 1.8, 'DisplayName', 'Deviation 3');
xlim([0 85]);
legend('Location', 'north', 'Orientation', 'horizontal', 'NumColumns', 4);
ylim([60 140]);
xlabel('Time step (k)');
ylabel('Tracking performance');
set(gca, 'FontSize', font_size, 'FontName', font_family);
% title('Tracking performance');


