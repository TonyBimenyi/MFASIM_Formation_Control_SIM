clear all;
close all;

% Define parameters
file_path = 'formation_exp_data.txt'; % Replace with your file's path

% Load the content of the data
file_content = fileread(file_path);

% Extract the CH3, CH1, CH2, CH4 data using regular expressions
ch3_data_match = regexp(file_content, 'CH3_Data_OutPut\[\d+\]=\{(.*?)\};', 'tokens', 'dotall');
ch1_data_match = regexp(file_content, 'CH1_Data_OutPut\[\d+\]=\{(.*?)\};', 'tokens', 'dotall');
ch2_data_match = regexp(file_content, 'CH2_Data_OutPut\[\d+\]=\{(.*?)\};', 'tokens', 'dotall');
ch4_data_match = regexp(file_content, 'CH4_Data_OutPut\[\d+\]=\{(.*?)\};', 'tokens', 'dotall');
ch5_data_match = regexp(file_content, 'CH5_Data_OutPut\[\d+\]=\{(.*?)\};', 'tokens', 'dotall');
ch6_data_match = regexp(file_content, 'CH6_Data_OutPut\[\d+\]=\{(.*?)\};', 'tokens', 'dotall');
ch7_data_match = regexp(file_content, 'CH7_Data_OutPut\[\d+\]=\{(.*?)\};', 'tokens', 'dotall');

if ~isempty(ch3_data_match) && ~isempty(ch1_data_match) && ~isempty(ch2_data_match) && ~isempty(ch4_data_match)
    % Convert the extracted data to numerical arrays
    ch3_data = str2num(ch3_data_match{1}{1});
    ch1_data = str2num(ch1_data_match{1}{1});
    ch2_data = str2num(ch2_data_match{1}{1});
    ch4_data = str2num(ch4_data_match{1}{1});
    ch5_data = str2num(ch5_data_match{1}{1});
    ch6_data = str2num(ch6_data_match{1}{1});
    ch7_data = str2num(ch7_data_match{1}{1});
    
    % Calculate distributed errors (xi)
    xi1 = ch2_data - ch1_data; % Error for Motor 1
    xi2 = ch2_data - ch3_data; % Error for Motor 2
    xi3 = ch2_data - ch4_data; % Error for Motor 3
    
    % Calculate MSE for each error
    mse1 = mean(xi1.^2); % MSE for Motor 1
    mse2 = mean(xi2.^2); % MSE for Motor 2
    mse3 = mean(xi3.^2); % MSE for Motor 3
    
    % Display MSE values
    fprintf('MSE for Motor 1 (CH2 - CH1): %.10e\n', mse1);
    fprintf('MSE for Motor 2 (CH2 - CH3): %.10e\n', mse2);
    fprintf('MSE for Motor 3 (CH2 - CH4): %.10e\n', mse3);
    
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
xlim([0 250]);
legend('Location', 'north', 'Orientation', 'horizontal', 'NumColumns', 4);
ylim([60 140]);
xlabel('Time step (k)');
ylabel('Tracking performance');
set(gca, 'FontSize', font_size, 'FontName', font_family);
% title('Tracking performance');

    % zoom_x_start = 940; % Start of zoomed x-range
    % zoom_x_end = 980; % End of zoomed x-range
    % axes('Position', [0.23,0.577,0.30,0.25]);
    % box on; hold on;
    % plot(ch2_data, '--b', 'LineWidth', 2.5);
    % plot(ch3_data, '--k', 'LineWidth', 2.5);
    % plot(ch1_data, '-.m', 'LineWidth', 2.5);

    % plot(ch4_data, '-r', 'LineWidth', 2.5);
    % plot(ch4_data, '--g', 'LineWidth', 2.5);
    % xlim([zoom_x_start zoom_x_end]);

    % yticks([-0.4,0,0.2]);
    set(gca, 'FontSize', font_size);
    
    
    % Save the plot
    print('-dpng', 'motor_plot.png');
else
    disp('CH3_Data_OutPut[1707] not found in the file.');
end