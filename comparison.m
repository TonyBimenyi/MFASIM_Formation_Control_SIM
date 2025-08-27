

clc; clear; close all;

% Parameters tuned for formation control
rho = 115.5;        % Increase from 60 (faster MFA convergence)
lamda = 900;     % Reduce regularization (was 300) for more responsiveness
eta = 58;         % Slightly faster Phi update

mu = 1000;      % Phi update regularization (unchanged)
epsilon = 1e-5; % Stability threshold (unchanged)
alpha = 0.5;    % Integral gain (unchanged)
T = 0.01;        % Sampling time (unchanged)
gamma1 = 0.32;   % was 0.07
gamma2 = 0.12;   % was 0.07
gamma3 = 0.18;   % was 0.1
gamma4 = 0.2;   % was 0.1

beta =45;      % Stronger SMC gain
sigma = 95;     % SMC regularization (unchanged)
tau = 1e-6;     % Slightly increased for smoother control
nena = 1e-5;    % Unchanged
rT = 1024;      % Unchanged
L = 900;        % Unchanged
m = 200;        % Changed to 200 simulation steps
n = 1000;        % Input scaling
timeV = 300;

% Time vector for plotting


% Initialize arrays
phi1 = zeros(m+1, 1); phi2 = zeros(m+1, 1); phi3 = zeros(m+1, 1); phi4 = zeros(m+1, 1);
mfa1 = zeros(m+1, 1); mfa2 = zeros(m+1, 1); mfa3 = zeros(m+1, 1); mfa4 = zeros(m+1, 1);
sm1 = zeros(m, 1); sm2 = zeros(m, 1); sm3 = zeros(m, 1); sm4 = zeros(m, 1);
y1 = zeros(m+1, 1); y2 = zeros(m+1, 1); y3 = zeros(m+1, 1); y4 = zeros(m+1, 1);
mfa_y1 = zeros(m+1, 1); mfa_y2 = zeros(m+1, 1); mfa_y3 = zeros(m+1, 1); mfa_y4 = zeros(m+1, 1);
u1 = zeros(m, 1); u2 = zeros(m, 1); u3 = zeros(m, 1); u4 = zeros(m, 1);
xi1 = zeros(m, 1); xi2 = zeros(m, 1); xi3 = zeros(m, 1); xi4 = zeros(m, 1);
v1 = zeros(m, 1); v2 = zeros(m, 1); v3 = zeros(m, 1); v4 = zeros(m, 1);
s1 = zeros(m, 1); s2 = zeros(m, 1); s3 = zeros(m, 1); s4 = zeros(m, 1);
omega1 = zeros(m, 1); omega2 = zeros(m, 1); omega3 = zeros(m, 1); omega4 = zeros(m, 1);
yd = zeros(m+1,1);

for k = 1:1:m+1
    yd(k) = 2 * sin(k * pi / 50) * exp(-0.01 * k);
end

for k = 1:1:m+1
    yd(k) = 3 + 0.5 * sin(0.1 * k);  % Always > 0

end

% yd = 3*ones(m+1,1);
% block = 50;
% levels = [3.0 3.5 3.0 2.5 3.0 3.5];       % up-down variation
% for i = 1:numel(levels)
%     i1 = (i-1)*block + 1;
%     i2 = min(i*block, m+1);
%     yd(i1:i2) = levels(i);
% end



% Main Loop
for k = 1:m
    % Phi updates
    if k == 1
        phi1(k) = 0.65; phi2(k) = 0.65; phi3(k) = 0.65; phi4(k) = 0.65;
    elseif k == 2
        phi1(k) = phi1(k-1) + (eta * u1(k-1) / (mu + u1(k-1)^2)) * (y1(k) - phi1(k-1)*u1(k-1));
        phi2(k) = phi2(k-1) + (eta * u2(k-1) / (mu + u2(k-1)^2)) * (y2(k) - phi2(k-1)*u2(k-1));
        phi3(k) = phi3(k-1) + (eta * u3(k-1) / (mu + u3(k-1)^2)) * (y3(k) - phi3(k-1)*u3(k-1));
        phi4(k) = phi4(k-1) + (eta * u4(k-1) / (mu + u4(k-1)^2)) * (y4(k) - phi4(k-1)*u4(k-1));
    else
        phi1(k) = phi1(k-1) + (eta * (u1(k-1) - u1(k-2)) / (mu + (u1(k-1) - u1(k-2))^2)) * (y1(k) - y1(k-1) - phi1(k-1) * (u1(k-1) - u1(k-2)));
        phi2(k) = phi2(k-1) + (eta * (u2(k-1) - u2(k-2)) / (mu + (u2(k-1) - u2(k-2))^2)) * (y2(k) - y2(k-1) - phi2(k-1) * (u2(k-1) - u2(k-2)));
        phi3(k) = phi3(k-1) + (eta * (u3(k-1) - u3(k-2)) / (mu + (u3(k-1) - u3(k-2))^2)) * (y3(k) - y3(k-1) - phi3(k-1) * (u3(k-1) - u3(k-2)));
        phi4(k) = phi4(k-1) + (eta * (u4(k-1) - u4(k-2)) / (mu + (u4(k-1) - u4(k-2))^2)) * (y4(k) - y4(k-1) - phi4(k-1) * (u4(k-1) - u4(k-2)));
    end

    % Stability protection
    if k > 2 && (abs(phi1(k)) <= epsilon || abs(u1(k-1) - u1(k-2)) <= epsilon || sign(phi1(k)) ~= sign(phi1(1)))
        phi1(k) = phi1(1);
    end
    if k > 2 && (abs(phi2(k)) <= epsilon || abs(u2(k-1) - u2(k-2)) <= epsilon || sign(phi2(k)) ~= sign(phi2(1)))
        phi2(k) = phi2(1);
    end
    if k > 2 && (abs(phi3(k)) <= epsilon || abs(u3(k-1) - u3(k-2)) <= epsilon || sign(phi3(k)) ~= sign(phi3(1)))
        phi3(k) = phi3(1);
    end
    if k > 2 && (abs(phi4(k)) <= epsilon || abs(u4(k-1) - u4(k-2)) <= epsilon || sign(phi4(k)) ~= sign(phi4(1)))
        phi4(k) = phi4(1);
    end

    

    v1(k) = 2 + 0.2 * sin(k * pi / 210) * exp(-0.01 * k) ;   % Top agent (further above leader)
    v2(k) = 1 + 0.2 * sin(k * pi / 800) * exp(-0.01 * k) ; % Further above leader
    v3(k) = -1 + 0.2 * sin(k * pi / 210) * exp(-0.01 * k) ; % Further below leader
    v4(k) = -2 + 0.2 * sin(k * pi / 210) * exp(-0.01 * k) ;  % Bottom agent (further below leader)

    % v1(k) = 2.0;
    % v2(k) = -2.0;
    % v3(k) = 2.2;
    % v4(k) = -2.2;

    % % Error dynamics
    xi1(k) = yd(k) - 2*y1(k) + y4(k) + 6.11;
    xi2(k) = y1(k) - 2*y2(k) + y3(k) + 0.9;
    xi3(k) = y2(k) + yd(k) - 2*y3(k) - 2.855;
    xi4(k) = y1(k) + y3(k) - 2*y4(k) - 4.999;

    % xi1(k) = yd(k) - 2*y1(k) + y4(k) ;
    % xi2(k) = y1(k) - 2*y2(k) + y3(k);
    % xi3(k) = y2(k) + yd(k) - 2*y3(k);
    % xi4(k) = y1(k) + y3(k) - 2*y4(k);
    

    if k == 1
        integral_xi1 = 0;
        integral_xi2 = 0;
        integral_xi3 = 0;
        integral_xi4 = 0;
    else
        integral_xi1 = integral_xi1 + T * xi4(k);
        integral_xi2 = integral_xi2 + T * (xi1(k) + xi3(k));
        integral_xi3 = integral_xi3 + T * xi2(k);
        integral_xi4 = integral_xi4 + T * (xi1(k) + xi3(k));
    end

    % Sliding surfaces
    if k == 1
        s1(k) = 0; s2(k) = 0; s3(k) = 0; s4(k) = 0;
    else
        s1(k) = xi1(k) + alpha * integral_xi1;
        s2(k) = xi2(k) + alpha * integral_xi2;
        s3(k) = xi3(k) + alpha * integral_xi3;
        s4(k) = xi4(k) + alpha * integral_xi4;
    end

    % MFA updates
    if k == 1
        mfa1(k) = 0; mfa2(k) = 0; mfa3(k) = 0; mfa4(k) = 0;
    else
        mfa1(k) = mfa1(k-1) + (rho * phi1(k)) / (lamda + abs(phi1(k)^2)) * xi1(k);
        mfa2(k) = mfa2(k-1) + (rho * phi2(k)) / (lamda + abs(phi2(k)^2)) * xi2(k);
        mfa3(k) = mfa3(k-1) + (rho * phi3(k)) / (lamda + abs(phi3(k)^2)) * xi3(k);
        mfa4(k) = mfa4(k-1) + (rho * phi4(k)) / (lamda + abs(phi4(k)^2)) * xi4(k);
    end

    % SMC updates
    if k == 1
        sm1(k) = 0; sm2(k) = 0; sm3(k) = 0; sm4(k) = 0;
    else
        sm1(k) = sm1(k-1) + (beta * phi1(k)) / (sigma + (phi1(k))^2) * ...
            ((xi1(k) + (y4(k) - y4(k-1)) + 1*  (yd(k+1) - yd(k)) + ((v4(k)-v4(k-1))-(v1(k)-v1(k-1))) + 1* (v1(k)-v1(k-1))) / (1 + 1) ...
            - (xi1(k) - s1(k)) / ((1+alpha*T) * 2) + tau * sign(s1(k)));

        sm2(k) = sm2(k-1) + (beta * phi2(k)) / (sigma + (phi2(k))^2) * ...
            ((xi2(k) + (y1(k) - y1(k-1) + y3(k) - y3(k-1)) + 0*(yd(k+1) - yd(k)) + (v1(k)-v1(k-1)+v3(k)-v3(k-1))-((v2(k)-v2(k-1)))  + 0*(v2(k)-v2(k-1))) / 2 ...
            - (xi2(k) - s2(k)) / ((1+alpha*T) * 2) + tau * sign(s2(k)));

        sm3(k) = sm3(k-1) + (beta * phi3(k)) / (sigma + (phi3(k))^2) * ...
            ((xi3(k) + (y2(k) - y2(k-1)) + 1* (yd(k+1) - yd(k)) + ((v2(k)-v2(k-1))-(v3(k)-v3(k-1))) + 1* (v3(k)-v3(k-1))) / (1 + 1) ...
            - (xi3(k) - s3(k)) / ((1+alpha*T) * 2) + tau * sign(s3(k)));

        sm4(k) = sm4(k-1) + (beta * phi4(k)) / (sigma + (phi4(k))^2) * ...
            ((xi4(k) + (y1(k) - y1(k-1) + y3(k) - y3(k-1)) + 0*(yd(k+1) - yd(k)) + ((v1(k)-v1(k-1)+v3(k)-v3(k-1))-(v4(k)-v4(k-1))) + 0*(v1(k)-v1(k-1))) / 2 ...
            - (xi4(k) - s4(k)) / ((1+alpha*T) * 2) + tau * sign(s4(k)));

            % sm1(k) = sm1(k-1) + (beta * phi1(k)) / (sigma + (phi1(k))^2) * ...
            %     ((xi1(k) + (y4(k) - y4(k-1)) + 1*  (yd(k+1) - yd(k)) + ((v1(k)-v1(k-1))-(v4(k)-v4(k-1))) + 1* (v1(k)-v1(k-1))) / (1 + 1) ...
            %     - (xi1(k) - s1(k)) / ((1+alpha*T) * 2) + tau * sign(s1(k)));
    
            % sm2(k) = sm2(k-1) + (beta * phi2(k)) / (sigma + (phi2(k))^2) * ...
            %     ((xi2(k) + (y1(k) - y1(k-1) + y3(k) - y3(k-1)) + 0*(yd(k+1) - yd(k)) + ((v2(k)-v2(k-1))) - (v1(k)-v1(k-1)+v3(k)-v3(k-1)) + 0*(v2(k)-v2(k-1))) / 2 ...
            %     - (xi2(k) - s2(k)) / ((1+alpha*T) * 2) + tau * sign(s2(k)));
    
            % sm3(k) = sm3(k-1) + (beta * phi3(k)) / (sigma + (phi3(k))^2) * ...
            %     ((xi3(k) + (y2(k) - y2(k-1)) + 1* (yd(k+1) - yd(k)) + ((v3(k)-v3(k-1))-(v2(k)-v2(k-1))) + 1* (v3(k)-v3(k-1))) / (1 + 1) ...
            %     - (xi3(k) - s3(k)) / ((1+alpha*T) * 2) + tau * sign(s3(k)));
    
            % sm4(k) = sm4(k-1) + (beta * phi4(k)) / (sigma + (phi4(k))^2) * ...
            %     ((xi4(k) + (y1(k) - y1(k-1) + y3(k) - y3(k-1)) + 0*(yd(k+1) - yd(k)) + ((v4(k)-v4(k-1))-(v1(k)-v1(k-1)+v3(k)-v3(k-1))) + 0*(v1(k)-v1(k-1))) / 2 ...
            %     - (xi4(k) - s4(k)) / ((1+alpha*T) * 2) + tau * sign(s4(k)));
    end

    % Control signals
    if k == 1
        u1(k) = 0; u2(k) = 0; u3(k) = 0; u4(k) = 0;
    else
        u1(k) = mfa1(k) + gamma1 * sm1(k);
        u2(k) = mfa2(k) + gamma2 * sm2(k);
        u3(k) = mfa3(k) + gamma3 * sm3(k);
        u4(k) = mfa4(k) + gamma4 * sm4(k);
    end

    % System dynamics (zero disturbances assumed)
    if k == 1
        y1(k) = 4.8; y2(k) = 3.8; y3(k) = 1.8; y4(k) = 0.8;
    end
    a =0.5;
    a4 =0.08;
    a3 =0.3;
    % a2 =0.33;
    % nonlinearity1 = 3; % Coefficient for cubic nonlinearity
    nonlinearity4 = 1; % Coefficient for cubic nonlinearity
    nonlinearity1 = 2; % Coefficient for cubic nonlinearity

    
    y1(k+1) = a * y1(k) + (n / (rT * 0.7)) * u1(k) + nonlinearity1  ;  % Slightly slower than before
    y2(k+1) = a * y2(k) + (n / (rT * 0.1)) * u2(k) + nonlinearity1;  % Equal inertia to y1
    y3(k+1) = a3 * y3(k) + (n / (rT * 0.2)) * u3(k) + nonlinearity1;  % Keep above yd, but reduce overshoot
    y4(k+1) = a4 * y4(k) + (n / (rT * 0.2)) * u4(k) + nonlinearity4;  % Reduce aggressiveness


    mfa_y1(k+1) = 0.5 * y1(k) + (n / (rT * 0.7)) * mfa1(k) + nonlinearity1  ;  % Slightly slower than before
    mfa_y2(k+1) = a * y2(k) + (n / (rT * 0.1)) * mfa2(k) + nonlinearity1;  % Equal inertia to y1
    mfa_y3(k+1) = a3 * y3(k) + (n / (rT * 0.2)) * mfa3(k) + nonlinearity1;  % Keep above yd, but reduce overshoot
    mfa_y4(k+1) = a4 * y4(k) + (n / (rT * 0.2)) * mfa4(k) + nonlinearity4;  % Reduce aggressiveness


        % % Plant model update with nonlinear term and feedforward
        % a = 0.5;
        % b1 = 1.2 * n / (rT * 0.2);
        % b2 = 1.15 * n / (rT * 0.2);
        % b3 = 1.2 * n / (rT * 0.2);
        % b4 = 1.15 * n / (rT * 0.2);
        
        % nonlinearity1 = 0.02; % Coefficient for cubic nonlinearity
        % nonlinearity2 = 0.02; % Coefficient for cubic nonlinearity
        % nonlinearity3 = 0.02; % Coefficient for cubic nonlinearity
        % nonlinearity4 = 0.02; % Coefficient for cubic nonlinearity
        % ff_gain = 0.45; % Feedforward gain
        
        % % Add cubic nonlinearity and feedforward term
        % y1(k+1) = a * y1(k) + b1 * u1(k) - nonlinearity1 * y1(k)^3 + ff_gain;
        % y2(k+1) = a * y2(k) + b2 * u2(k) - nonlinearity2 * y2(k)^2 + ff_gain;
        % y3(k+1) = a * y3(k) + b3 * u3(k) - nonlinearity3 * y3(k)^3 + ff_gain;
        % y4(k+1) = a * y4(k) + b4 * u4(k) - nonlinearity4 * y4(k)^2 + ff_gain;
end

% Time vector for plotting
t_plot = (1:m);

% Compute expected outputs for formation tracking
expected1 = yd + [v1; v1(end)];
expected2 = yd + [v2; v2(end)];
expected3 = yd + [v3; v3(end)];
expected4 = yd + [v4; v4(end)];



% Plot
figure('Position', [100, 100, 1100, 600]);  % [x, y, width, height]
markersize = 7;

% Reference trajectory yd(k) with solid line + marker
plot(t_plot, yd(1:m), '-', ...
     'Color', [0 0 0], ...
     'LineWidth', 1.5, ...
     'MarkerIndices', 1:20:m, ...
     'DisplayName', 'y_d(k)'); hold on;

% Agent outputs (dashed + distinct colors + markers + marker size)
plot(t_plot, y1(1:m), '--', ...
     'Color', [0 0.4470 0.7410], ...
     'LineWidth', 2.85, ...
     'Marker', 's', ...
     'MarkerSize', markersize, ...
     'MarkerIndices', 5:9:m, ...
     'DisplayName', 'y_1(k)');


plot(t_plot, y2(1:m), '--', ...
     'Color', [0.8500 0.3250 0.0980], ...
     'LineWidth', 2.85, ...
     'Marker', '^', ...
     'MarkerSize', markersize, ...
     'MarkerIndices', 10:14:m, ...
     'DisplayName', 'y_2(k)');

plot(t_plot, y3(1:m), '--', ...
     'Color', [0.9290 0.6940 0.1250], ...
     'LineWidth', 2.85, ...
     'Marker', 'v', ...
     'MarkerSize', markersize, ...
     'MarkerIndices', 5:13:m, ...
     'DisplayName', 'y_3(k)');

plot(t_plot, y4(1:m), '--', ...
     'Color', [0.4940 0.1840 0.5560], ...
     'LineWidth', 2.85, ...
     'Marker', 's', ...
     'MarkerSize', markersize, ...  % Slightly larger for 'x' marker for visibility
     'MarkerIndices', 5:13:m, ...
     'DisplayName', 'y_4(k)');


% Expected outputs (solid lines, lighter/different hues)
plot(t_plot, expected1(1:m), '-', ...
     'Color', [0.3010 0.7450 0.9330], ...
     'LineWidth', 1.5, ...
     'DisplayName', 'v_1(k)');

plot(t_plot, expected2(1:m), '-', ...
     'Color', [0.6350 0.0780 0.1840], ...
     'LineWidth', 1.5, ...
     'DisplayName', 'v_2(k)');

plot(t_plot, expected3(1:m), '-', ...
     'Color', [0.4660 0.6740 0.1880], ...
     'LineWidth', 1.5, ...
     'DisplayName', 'v_3(k)');

plot(t_plot, expected4(1:m), '-', ...
    'Color', [1.0 0.4 0.7], ...  % pink
        'LineWidth', 1.5, ...
        'DisplayName', 'v_4(k)');

% % Plot actual outputs (thicker lines, semi-transparent)
% plot(t_plot, expected1(1:m), '-', ...
%     'Color', [0.3010 0.7450 0.9330 0.1], ... % RGBA: blue, 40% opacity
%     'LineWidth', 12, ...
%     'DisplayName', 'v1(k)');

% plot(t_plot, expected2(1:m), '-', ...
%     'Color', [0.6350 0.0780 0.1840 0.1], ... % dark red
%     'LineWidth', 12, ...
%     'DisplayName', 'v2(k)');

% plot(t_plot, expected3(1:m), '-', ...
%     'Color', [0.4660 0.6740 0.1880 0.1], ... % green
%     'LineWidth', 12, ...
%     'DisplayName', 'v3(k)');

% plot(t_plot, expected4(1:m), '-', ...
%     'Color', [1.0 0.4 0.7 0.1], ... % pink
%     'LineWidth', 12, ...
%     'DisplayName', 'v4(k)');
   

% Legend, Axis, Labels
font_size = 20;
font_family = 'Times New Roman';
legend('Location', 'north', 'Orientation', 'horizontal', 'FontSize', font_size, 'FontName', font_family);
ylim([0 7]);       
xlabel('Time step (k)','FontSize', font_size, 'FontName', font_family);
ylabel('Tracking performance','FontSize', font_size, 'FontName', font_family);

% ✅ Change tick font size and font
set(gca, 'FontSize', font_size, 'FontName', font_family);
grid off;
% title('Agent Outputs and Expected Formation Tracking vs Reference yd');



figure('Position', [100, 100, 1100, 600]);  % [x, y, width, height]
hold on; grid off;
plot(t_plot, xi1, 'r--', 'LineWidth', 1.5);
plot(t_plot, xi2, 'g--', 'LineWidth', 1.5);
plot(t_plot, xi3, 'b--', 'LineWidth', 1.5);
plot(t_plot, xi4, 'm--', 'LineWidth', 1.5);
legend('\xi_1(k)', '\xi_2(k)', '\xi_3(k)', '\xi_4(k)');
xlabel('Time step (k)');
ylabel('Formation error');
set(gca, 'FontSize', font_size, 'FontName', font_family);


% --- Small inset axes (e.g., zoom-in view) ---
axes('Position',[0.30 0.55 0.35 0.3]);  % [x y width height] (normalized units)
box on; hold on; grid off;

% Replot inside small axes (maybe focus on first 50 steps)
idx = (t_plot >= 150 & t_plot <= 200);  % Zoom region
plot(t_plot(idx), xi1(idx), 'r--', 'LineWidth', 1.2);
plot(t_plot(idx), xi2(idx), 'g--', 'LineWidth', 1.2);
plot(t_plot(idx), xi3(idx), 'b--', 'LineWidth', 1.2);
plot(t_plot(idx), xi4(idx), 'm--', 'LineWidth', 1.2);

set(gca, 'FontSize', font_size-2, 'FontName', font_family);

% figure('Position', [100, 100, 1100, 600]);  % [x, y, width, height]

% % --- Subplot 1: xi1 ---
% subplot(2,2,1); hold on; grid off;
% plot(t_plot, xi1, 'r--', 'LineWidth', 1.5);
% ylim([-3 10]);
% xlabel('Time step (k)');
% ylabel('Formation error');
% set(gca, 'FontSize', font_size, 'FontName', font_family);

% % Inset zoom for xi1
% axes('Position',[0.20 0.70 0.15 0.15]);  % adjust location
% box on; hold on; grid off;
% idx = (t_plot >= 150 & t_plot <= 200);
% plot(t_plot(idx), xi1(idx), 'r--', 'LineWidth', 1.2);
% set(gca, 'FontSize', font_size-3, 'FontName', font_family);

% % --- Subplot 2: xi2 ---
% subplot(2,2,2); hold on; grid off;
% plot(t_plot, xi2, 'g--', 'LineWidth', 1.5);
% xlabel('Time step (k)');
% ylabel('Formation error');
% set(gca, 'FontSize', font_size, 'FontName', font_family);

% % Inset zoom for xi2
% axes('Position',[0.63 0.70 0.15 0.15]);
% box on; hold on; grid off;
% plot(t_plot(idx), xi2(idx), 'g--', 'LineWidth', 1.2);
% set(gca, 'FontSize', font_size-3, 'FontName', font_family);

% % --- Subplot 3: xi3 ---
% subplot(2,2,3); hold on; grid off;
% plot(t_plot, xi3, 'b--', 'LineWidth', 1.5);
% xlabel('Time step (k)');
% ylabel('Formation error');
% set(gca, 'FontSize', font_size, 'FontName', font_family);

% % Inset zoom for xi3
% axes('Position',[0.20 0.23 0.15 0.15]);
% box on; hold on; grid off;
% plot(t_plot(idx), xi3(idx), 'b--', 'LineWidth', 1.2);
% set(gca, 'FontSize', font_size-3, 'FontName', font_family);

% % --- Subplot 4: xi4 ---
% subplot(2,2,4); hold on; grid off;
% plot(t_plot, xi4, 'm--', 'LineWidth', 1.5);
% xlabel('Time step (k)');
% ylabel('Formation error');
% set(gca, 'FontSize', font_size, 'FontName', font_family);

% % Inset zoom for xi4
% axes('Position',[0.63 0.23 0.15 0.15]);
% box on; hold on; grid off;
% plot(t_plot(idx), xi4(idx), 'm--', 'LineWidth', 1.2);
% set(gca, 'FontSize', font_size-3, 'FontName', font_family);


% sgtitle('Internal Signals of All Agents'); % global title


%% === MSI of xi for Hybrid (MFA+SMC) vs MFA-only ===
% Define MSI as mean of squared xi over k = 1..m
msi = @(x) mean(x(1:m).^2);

% xi for MFA-only (use the same formation-error definitions but with mfa_y*)
xi_mfa1 = zeros(m,1); xi_mfa2 = zeros(m,1);
xi_mfa3 = zeros(m,1); xi_mfa4 = zeros(m,1);

for k = 1:m
    % Use the same offsets you used above (+6.11, +0.9, -2.855, -4.999)
    xi_mfa1(k) = yd(k) - 2*mfa_y1(k) + mfa_y4(k) + 6.11;
    xi_mfa2(k) = mfa_y1(k) - 2*mfa_y2(k) + mfa_y3(k) + 0.9;
    xi_mfa3(k) = mfa_y2(k) + yd(k) - 2*mfa_y3(k) - 2.855;
    xi_mfa4(k) = mfa_y1(k) + mfa_y3(k) - 2*mfa_y4(k) - 4.999;
end

% MSI per agent (Hybrid)
MSI_hybrid = struct();
MSI_hybrid.xi1 = msi(xi1);
MSI_hybrid.xi2 = msi(xi2);
MSI_hybrid.xi3 = msi(xi3);
MSI_hybrid.xi4 = msi(xi4);
MSI_hybrid.avg = mean([MSI_hybrid.xi1, MSI_hybrid.xi2, MSI_hybrid.xi3, MSI_hybrid.xi4]);

% MSI per agent (MFA-only)
MSI_mfa = struct();
MSI_mfa.xi1 = msi(xi_mfa1);
MSI_mfa.xi2 = msi(xi_mfa2);
MSI_mfa.xi3 = msi(xi_mfa3);
MSI_mfa.xi4 = msi(xi_mfa4);
MSI_mfa.avg = mean([MSI_mfa.xi1, MSI_mfa.xi2, MSI_mfa.xi3, MSI_mfa.xi4]);

% Pretty print
fprintf('\n=== Mean Squared Index (MSI) of \\xi ===\n');
fprintf('Hybrid (MFA+SMC):  xi1 = %.6g, xi2 = %.6g, xi3 = %.6g, xi4 = %.6g,  AVG = %.6g\n', ...
    MSI_hybrid.xi1, MSI_hybrid.xi2, MSI_hybrid.xi3, MSI_hybrid.xi4, MSI_hybrid.avg);
fprintf('MFA-only:          xi1 = %.6g, xi2 = %.6g, xi3 = %.6g, xi4 = %.6g,  AVG = %.6g\n\n', ...
    MSI_mfa.xi1, MSI_mfa.xi2, MSI_mfa.xi3, MSI_mfa.xi4, MSI_mfa.avg);

    