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

% Convert from cell array of strings to numeric (skip the first entry 'CHx:')
ch1_data = str2double(ch1_line(2:end));
ch2_data = str2double(ch2_line(2:end));
ch3_data = str2double(ch3_line(2:end));
ch4_data = str2double(ch4_line(2:end));
ch5_data = str2double(ch5_line(2:end));

% Example: Plot
figure;
plot(ch1_data, '-m', 'DisplayName', 'CH1');
hold on;
plot(ch2_data, '-b', 'DisplayName', 'CH2');
plot(ch3_data, '-k', 'DisplayName', 'CH3');
plot(ch4_data, '-r', 'DisplayName', 'CH4');
plot(ch5_data, '-g', 'DisplayName', 'CH5');
xlim([0 250]);
legend('Location', 'north', 'Orientation', 'horizontal', 'NumColumns', 5);
ylim([60 140]);
xlabel('Time step (k)');
ylabel('Output');
title('Tracking performance');
