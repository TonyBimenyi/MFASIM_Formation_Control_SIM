clc; clear; close all;

% Parameters tuned for formation control
rho = 40;       % MFA adaptation speed (faster convergence)
eta = 2;        % Phi adaptation speed (faster phi estimation)
lamda = 300;    % MFA regularization (less restrictive)
mu = 80.5;      % Phi update regularization (unchanged)
epsilon = 1e-5; % Stability threshold (unchanged)
alpha = 0.5;    % Integral gain (unchanged)
T = 0.1;        % Sampling time (unchanged)
gamma1 = 0.07;  % Increased SMC weight
gamma2 = 0.07;
gamma3 = 0.1;
gamma4 = 0.1;
beta = 15;      % Stronger SMC gain
sigma = 95;     % SMC regularization (unchanged)
tau = 1e-4;     % Slightly increased for smoother control
nena = 1e-5;    % Unchanged
rT = 1024;      % Unchanged
L = 900;        % Unchanged
m = 100;        % Simulation steps
n = 700;        % Input scaling


% Time vector for plotting
k = 1:m;
t = (k-1) * T;

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
yd = zeros(m+1, 1);

% Desired signal (Reference trajectory) y0(k)
% for k = 1:m+1
%     yd(k) = 0.5 * sin(k * pi / 25) + 0.2 * cos(k * pi / 10);
% end
for k = 1:m+1
    yd(k) = 0.7 * sin(k * pi / 50);
end


% Main Loop
for k = 1:m-1 % Loop to m-1 to avoid index out of bounds
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

 

    v1(k) = 2 + 0.1 * sin(k * pi / 40);   % Top agent (above leader)
    v2(k) = 1.6 + 0.05 * sin(k * pi / 30);  % Slightly above leader
    v3(k) = -1.6 - 0.05 * sin(k * pi / 30); % Slightly below leader
    v4(k) = -2 - 0.1 * sin(k * pi / 40);  % Bottom agent (below leader)
    
    
    % Error dynamics
    xi1(k) = yd(k) - 2*y1(k) + y4(k) - v1(k);
    xi2(k) = y1(k) - 2*y2(k) + y3(k) - v2(k);
    xi3(k) = y2(k) + yd(k) - 2*y3(k) - v3(k);
    xi4(k) = y1(k) + y3(k) - 2*y4(k) - v4(k);

    if k == 1
        integral_xi1 = 0;
        integral_xi2 = 0;
        integral_xi3 = 0;
        integral_xi4 = 0;
    else
        integral_xi1 = integral_xi1 + T * xi4(k);
        integral_xi2 = integral_xi2 + T * xi1(k)+xi3(k);
        integral_xi3 = integral_xi3 + T * xi2(k);
        integral_xi4 = integral_xi4 + T * xi1(k)+xi3(k);
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
            ( (xi1(k) + (y4(k) - y4(k-1)) + (yd(k+1) - yd(k)) + ((v1(k)-v1(k-1))-(v4(k)-v4(k-1))) + (v1(k)-v1(k-1))) / (1 + 1) ...
            - (xi1(k) - s1(k)) / ((1+alpha*T) * 2) + tau * sign(s1(k)) );

        sm2(k) = sm2(k-1) + (beta * phi2(k)) / (sigma + (phi2(k))^2) * ...
            ( (xi2(k) + (y1(k) - y1(k-1) + y3(k) - y3(k-1)) + 0*(yd(k+1) - yd(k)) + ((v2(k)-v2(k-1))-(v1(k)-v1(k-1)+v3(k)-v3(k-1))) + 0*(v1(k)-v1(k-1))) / 2 ...
            - (xi2(k) - s2(k)) / ((1+alpha*T) * 2) + tau * sign(s2(k)) );

        sm3(k) = sm3(k-1) + (beta * phi3(k)) / (sigma + (phi3(k))^2) * ...
            ( (xi3(k) + (y2(k) - y2(k-1)) + (yd(k+1) - yd(k)) + ((v3(k)-v3(k-1))-(v2(k)-v2(k-1))) + (v3(k)-v3(k-1))) / (1 + 1) ...
            - (xi3(k) - s3(k)) / ((1+alpha*T) * 2) + tau * sign(s3(k)) );

        sm4(k) = sm4(k-1) + (beta * phi4(k)) / (sigma + (phi4(k))^2) * ...
            ( (xi4(k) + (y1(k) - y1(k-1) + y3(k) - y3(k-1)) + 0*(yd(k+1) - yd(k)) + ((v4(k)-v4(k-1))-(v1(k)-v1(k-1)+v3(k)-v3(k-1))) + 0*(v1(k)-v1(k-1))) / 2 ...
            - (xi4(k) - s4(k)) / ((1+alpha*T) * 2) + tau * sign(s4(k)) );
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
    y1(k+1) = 0.1 * y1(k) + (n / (rT * 0.05)) * u1(k);  % Fast response (close to yd)
    y2(k+1) = 0.1 * y2(k) + (n / (rT * 0.05)) * u2(k);  % Slightly slower than y1 (below yd)
    y3(k+1) = 0.2 * y3(k) + (n / (rT * 0.1)) * u3(k);  % Slightly above yd
    y4(k+1) = 0.2 * y4(k) + (n / (rT * 0.1)) * u4(k);  % Further below yd
    
end

% Time vector for plotting (from 0 to m*T)
t_plot = (0:m) * T;

% Plot y1, y2, y3, y4 over time
figure;
plot(t_plot, y1, 'LineWidth', 1.5); hold on;
plot(t_plot, y2, 'LineWidth', 1.5);
plot(t_plot, y3, 'LineWidth', 1.5);
plot(t_plot, y4, 'LineWidth', 1.5);
plot(t_plot, yd, 'LineWidth', 1.5, 'LineStyle', '--');
xlabel('Time (s)');
ylabel('Output y_i');
legend('y1', 'y2', 'y3', 'y4', 'yd');
title('Formation Control: Outputs y1, y2, y3, y4 vs Reference yd');
grid on;

% % Plot tracking errors
% figure;
% plot(t(1:m), xi1, 'LineWidth', 1.5); hold on;
% plot(t(1:m), xi2, 'LineWidth', 1.5);
% plot(t(1:m), xi3, 'LineWidth', 1.5);
% plot(t(1:m), xi4, 'LineWidth', 1.5);
% plot(t(1:m), v1, 'LineWidth', 1.5, 'LineStyle', '--');
% plot(t(1:m), v2, 'LineWidth', 1.5, 'LineStyle', '--');
% plot(t(1:m), v3, 'LineWidth', 1.5, 'LineStyle', '--');
% plot(t(1:m), v4, 'LineWidth', 1.5, 'LineStyle', '--');
% xlabel('Time (s)');
% ylabel('Tracking Errors xi_i and Desired Deviations v_i');
% legend('xi1', 'xi2', 'xi3', 'xi4', 'v1', 'v2', 'v3', 'v4');
% title('Tracking Errors vs Desired Deviations');
% grid on;