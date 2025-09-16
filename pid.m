%% PID Obstacle Avoidance Simulation for Arduino Robot Car
% Simulates distance control using ultrasonic sensor + PID
% Robot stops if obstacle <= 10cm, else uses PID to maintain 25cm

clear; clc; close all;

%% Parameters
dt = 0.1;                   % Simulation time step (seconds)
t_end = 30;                 % Total simulation time (seconds)
time = 0:dt:t_end;

% PID Gains
Kp = 2.0;
Ki = 0.0;                   % Optional: add later
Kd = 0.0;                   % Optional: add later

% Robot Settings
setpoint = 25;              % Target distance (cm)
base_speed = 50;            % Base forward speed (0-255 scale)
stop_distance = 10;         % Stop if <= 10cm

% Initial Conditions
distance = 40;              % Start 40cm from obstacle
error = 0;
integral = 0;
prev_error = 0;
left_speed = base_speed;
right_speed = base_speed;

% Preallocate arrays for plotting
distances = zeros(size(time));
errors = zeros(size(time));
left_speeds = zeros(size(time));
right_speeds = zeros(size(time));
statuses = cell(size(time)); % "Moving" or "STOPPED"

%% Simulate Dynamic Obstacle (moves closer then away)
% Create "moving obstacle" — you can change this to static or user-defined
obstacle_position = 40 * ones(size(time));  % Start at 40cm
obstacle_position(time >= 5 & time < 15) = 40 - 2*(time(time >= 5 & time < 15) - 5);  % Move closer
obstacle_position(time >= 15 & time < 25) = 20 + 2*(time(time >= 15 & time < 25) - 15); % Move away

%% Simulation Loop
for i = 1:length(time)
    % Current obstacle distance (simulated sensor reading)
    distance = obstacle_position(i);
    distances(i) = distance;

    % 🛑 STOP if too close
    if distance <= stop_distance
        left_speed = 0;
        right_speed = 0;
        error = 0;
        integral = 0;
        prev_error = 0;
        statuses{i} = '🛑 STOPPED';
    else
        % ➤ PID CONTROL
        error = setpoint - distance;
        errors(i) = error;

        % P term
        P = Kp * error;

        % I term (optional)
        % integral = integral + error * dt;
        % I = Ki * integral;

        % D term (optional)
        % D = Kd * (error - prev_error) / dt;
        % prev_error = error;

        pid_output = P;  % + I + D;  % Add I and D later if needed
        pid_output = max(min(pid_output, 50), -50);  % Constrain

        % Calculate motor speeds
        left_speed = base_speed + pid_output;
        right_speed = base_speed - pid_output;

        % Constrain motor speeds
        left_speed = max(min(left_speed, 100), 0);
        right_speed = max(min(right_speed, 100), 0);

        statuses{i} = '🚗 MOVING';
    end

    left_speeds(i) = left_speed;
    right_speeds(i) = right_speed;
end

%% Plot Results
figure('Position', [100, 100, 1000, 800]);

% Distance + Setpoint
subplot(4,1,1);
plot(time, distances, 'b-', 'LineWidth', 2);
hold on;
plot(time, setpoint * ones(size(time)), 'r--', 'LineWidth', 1.5);
plot(time, stop_distance * ones(size(time)), 'k--', 'LineWidth', 1.5);
title('Obstacle Distance Over Time');
xlabel('Time (s)');
ylabel('Distance (cm)');
legend('Measured Distance', 'Setpoint (25cm)', 'Stop Line (10cm)');
grid on;

% Error
subplot(4,1,2);
plot(time, errors, 'm-', 'LineWidth', 2);
title('Error (Setpoint - Distance)');
xlabel('Time (s)');
ylabel('Error (cm)');
grid on;

% Motor Speeds
subplot(4,1,3);
plot(time, left_speeds, 'g-', 'LineWidth', 2);
hold on;
plot(time, right_speeds, 'c-', 'LineWidth', 2);
title('Motor Speeds (Left: Green, Right: Cyan)');
xlabel('Time (s)');
ylabel('Speed (0-100)');
legend('Left Motor', 'Right Motor');
grid on;

% Status (Stopped vs Moving)
subplot(4,1,4);
status_numeric = strcmp(statuses, '🛑 STOPPED');
stairs(time, status_numeric, 'LineWidth', 2, 'Color', [0.9 0.2 0.2]);
title('Robot Status: 1 = STOPPED, 0 = MOVING');
xlabel('Time (s)');
ylabel('Status');
ylim([-0.1, 1.1]);
grid on;

%% Optional: Animate Robot Movement (Simple)
figure;
for i = 1:10:length(time)
    clf;
    hold on;
    axis([-10 50 -10 10]);
    title(sprintf('Time: %.1fs | Distance: %.1fcm | %s', time(i), distances(i), statuses{i}));
    xlabel('Front View (Robot at 0, Obstacle at Distance)');

    % Draw robot (rectangle)
    rectangle('Position', [-5, -2, 10, 4], 'FaceColor', 'blue', 'EdgeColor', 'black');

    % Draw obstacle (red wall)
    if distances(i) < 50
        line([distances(i), distances(i)], [-5, 5], 'Color', 'red', 'LineWidth', 3);
        text(distances(i) + 1, 3, sprintf('Obstacle\n%.1f cm', distances(i)), 'FontSize', 10);
    end

    % Draw direction arrows if moving
    if left_speeds(i) > 0 || right_speeds(i) > 0
        quiver(0, 0, 5, 0, 'b', 'LineWidth', 2); % Forward arrow
        text(1, 1, sprintf('L:%d R:%d', left_speeds(i), right_speeds(i)), 'FontSize', 10);
    end

    if contains(statuses{i}, 'STOPPED')
        text(-3, -4, '🛑 STOPPED', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'red');
    end

    drawnow;
    pause(0.05);
end