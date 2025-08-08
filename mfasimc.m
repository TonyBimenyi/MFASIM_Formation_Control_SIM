clc; clear; close all;

% Parameters tuned for formation control
rho = 90;        % Increase from 60 (faster MFA convergence)
lamda = 900;     % Reduce regularization (was 300) for more responsiveness
eta = 28;         % Slightly faster Phi update

mu = 80.5;      % Phi update regularization (unchanged)
epsilon = 1e-5; % Stability threshold (unchanged)
alpha = 0.5;    % Integral gain (unchanged)
T = 0.06;        % Sampling time (unchanged)
gamma1 = 0.12;   % was 0.07
gamma2 = 0.12;   % was 0.07
gamma3 = 0.18;   % was 0.1
gamma4 = 0.18;   % was 0.1

beta =50;      % Stronger SMC gain
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
u1 = zeros(m, 1); u2 = zeros(m, 1); u3 = zeros(m, 1); u4 = zeros(m, 1);
xi1 = zeros(m, 1); xi2 = zeros(m, 1); xi3 = zeros(m, 1); xi4 = zeros(m, 1);
v1 = zeros(m, 1); v2 = zeros(m, 1); v3 = zeros(m, 1); v4 = zeros(m, 1);
s1 = zeros(m, 1); s2 = zeros(m, 1); s3 = zeros(m, 1); s4 = zeros(m, 1);
omega1 = zeros(m, 1); omega2 = zeros(m, 1); omega3 = zeros(m, 1); omega4 = zeros(m, 1);
yd = zeros(m+1,1);

for k = 1:1:m+1
    yd(k) = 2 * sin(k * pi / 50) * exp(-0.01 * k);
end


% Main Loop
for k = 1:m
    % Phi updates
    if k == 1
        phi1(k) = 0.5; phi2(k) = 0.5; phi3(k) = 0.5; phi4(k) = 0.5;
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

    

    v1(k) = 7 + 0.2 * sin(k * pi / 50) * exp(-0.01 * k) + 0.5;   % Top agent (further above leader)
    v2(k) = 5.2 + 0.2 * sin(k * pi / 50) * exp(-0.01 * k) + 0.2; % Further above leader
    v3(k) = 3.2 + 0.2 * sin(k * pi / 50) * exp(-0.01 * k) + 0.2; % Further below leader
    v4(k) = 1 + 0.2 * sin(k * pi / 50) * exp(-0.01 * k) + 0.5;  % Bottom agent (further below leader)

    % v1(k) = 2.0;
    % v2(k) = -2.0;
    % v3(k) = 2.2;
    % v4(k) = -2.2;

    % Error dynamics
    xi1(k) = yd(k) - 2*y1(k) + y4(k) + 1;
    xi2(k) = y1(k) - 2*y2(k) + y3(k) + 1;
    xi3(k) = y2(k) + yd(k) - 2*y3(k) + 1;
    xi4(k) = y1(k) + y3(k) - 2*y4(k) + 1;

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
            ((xi1(k) + (y4(k) - y4(k-1)) + (yd(k+1) - yd(k)) + ((v1(k)-v1(k-1))-(v4(k)-v4(k-1))) + (v1(k)-v1(k-1))) / (1 + 1) ...
            - (xi1(k) - s1(k)) / ((1+alpha*T) * 2) + tau * sign(s1(k)));

        sm2(k) = sm2(k-1) + (beta * phi2(k)) / (sigma + (phi2(k))^2) * ...
            ((xi2(k) + (y1(k) - y1(k-1) + y3(k) - y3(k-1)) + 0*(yd(k+1) - yd(k)) + ((v2(k)-v2(k-1))-(v1(k)-v1(k-1)+v3(k)-v3(k-1))) + 0*(v1(k)-v1(k-1))) / 2 ...
            - (xi2(k) - s2(k)) / ((1+alpha*T) * 2) + tau * sign(s2(k)));

        sm3(k) = sm3(k-1) + (beta * phi3(k)) / (sigma + (phi3(k))^2) * ...
            ((xi3(k) + (y2(k) - y2(k-1)) + (yd(k+1) - yd(k)) + ((v3(k)-v3(k-1))-(v2(k)-v2(k-1))) + (v3(k)-v3(k-1))) / (1 + 1) ...
            - (xi3(k) - s3(k)) / ((1+alpha*T) * 2) + tau * sign(s3(k)));

        sm4(k) = sm4(k-1) + (beta * phi4(k)) / (sigma + (phi4(k))^2) * ...
            ((xi4(k) + (y1(k) - y1(k-1) + y3(k) - y3(k-1)) + 0*(yd(k+1) - yd(k)) + ((v4(k)-v4(k-1))-(v1(k)-v1(k-1)+v3(k)-v3(k-1))) + 0*(v1(k)-v1(k-1))) / 2 ...
            - (xi4(k) - s4(k)) / ((1+alpha*T) * 2) + tau * sign(s4(k)));
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
        y1(k) = 0; y2(k) = 0; y3(k) = 0; y4(k) = 0;
    end
    y1(k+1) = 0.15 * y1(k) + (n / (rT * 0.2)) * u1(k);  % Slightly slower than before
    y2(k+1) = 0.15 * y2(k) + (n / (rT * 0.2)) * u2(k);  % Equal inertia to y1
    y3(k+1) = 0.18 * y3(k) + (n / (rT * 0.2)) * u3(k);  % Keep above yd, but reduce overshoot
    y4(k+1) = 0.18 * y4(k) + (n / (rT * 0.2)) * u4(k);  % Reduce aggressiveness
end

% Time vector for plotting
t_plot = (1:m);

% Compute expected outputs for formation tracking
expected1 = yd + [v1; v1(end)];
expected2 = yd + [v2; v2(end)];
expected3 = yd + [v3; v3(end)];
expected4 = yd + [v4; v4(end)];

% Plot
figure;

% Reference trajectory yd(k) with solid line + marker
plot(t_plot, yd(1:m), '-', ...
     'Color', [0 0 0], ...
     'LineWidth', 1.5, ...
     'MarkerIndices', 1:20:m, ...
     'DisplayName', 'yd(k)'); hold on;

% Agent outputs (dashed + distinct colors + markers)
plot(t_plot, y1(1:m), '--', ...
     'Color', [0 0.4470 0.7410], ...
     'LineWidth', 2.5, ...
     'Marker', 's', ...
     'MarkerIndices', 5:20:m, ...
     'DisplayName', 'y1(k)');

plot(t_plot, y2(1:m), '--', ...
     'Color', [0.8500 0.3250 0.0980], ...
     'LineWidth', 2.5, ...
     'Marker', '^', ...
     'MarkerIndices', 10:20:m, ...
     'DisplayName', 'y2(k)');

plot(t_plot, y3(1:m), '--', ...
     'Color', [0.9290 0.6940 0.1250], ...
     'LineWidth', 2.5, ...
     'Marker', 'v', ...
     'MarkerIndices', 15:20:m, ...
     'DisplayName', 'y3(k)');

plot(t_plot, y4(1:m), '--', ...
     'Color', [0.4940 0.1840 0.5560], ...
     'LineWidth', 2.5, ...
     'Marker', 'x', ...
     'MarkerIndices', 20:20:m, ...
     'DisplayName', 'y4(k)');

% Expected outputs (solid lines, lighter/different hues)
plot(t_plot, expected1(1:m), '-', ...
     'Color', [0.3010 0.7450 0.9330], ...
     'LineWidth', 1.5, ...
     'DisplayName', 'v1(k)');

plot(t_plot, expected2(1:m), '-', ...
     'Color', [0.6350 0.0780 0.1840], ...
     'LineWidth', 1.5, ...
     'DisplayName', 'v2(k)');

plot(t_plot, expected3(1:m), '-', ...
     'Color', [0.4660 0.6740 0.1880], ...
     'LineWidth', 1.5, ...
     'DisplayName', 'v3(k)');

plot(t_plot, expected4(1:m), '-', ...
    'Color', [1.0 0.4 0.7], ...  % pink
        'LineWidth', 1.5, ...
        'DisplayName', 'v4(k)');

% Plot actual outputs (thicker lines, semi-transparent)
plot(t_plot, expected1(1:m), '-', ...
    'Color', [0.3010 0.7450 0.9330 0.1], ... % RGBA: blue, 40% opacity
    'LineWidth', 12, ...
    'DisplayName', 'v1(k)');

plot(t_plot, expected2(1:m), '-', ...
    'Color', [0.6350 0.0780 0.1840 0.1], ... % dark red
    'LineWidth', 12, ...
    'DisplayName', 'v2(k)');

plot(t_plot, expected3(1:m), '-', ...
    'Color', [0.4660 0.6740 0.1880 0.1], ... % green
    'LineWidth', 12, ...
    'DisplayName', 'v3(k)');

plot(t_plot, expected4(1:m), '-', ...
    'Color', [1.0 0.4 0.7 0.1], ... % pink
    'LineWidth', 12, ...
    'DisplayName', 'v4(k)');
   

% Legend, Axis, Labels
legend('Location', 'north', 'Orientation', 'horizontal', 'NumColumns', 5);
% ylim([-6 7.5]);       
xlabel('Time step (k)');
ylabel('Tracking performance');
grid on;
% title('Agent Outputs and Expected Formation Tracking vs Reference yd');

% figure; hold on; grid on;
% plot(t_plot, xi1, 'r--', 'LineWidth', 1.5);
% plot(t_plot, xi2, 'g--', 'LineWidth', 1.5);
% plot(t_plot, xi3, 'b--', 'LineWidth', 1.5);
% plot(t_plot, xi4, 'm--', 'LineWidth', 1.5);
% legend('xi1(k)', 'xi2(k)', 'xi3(k)', 'xi4(k)');
% xlabel('Time step (k)');
% ylabel('Internal Variable (e.g., integral of error)');
% title('Internal Signals: \xi_i(k)');
