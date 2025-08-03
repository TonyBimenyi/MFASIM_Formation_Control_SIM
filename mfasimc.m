clc; clear; close all;

% Parameters
rho = 10.5;
eta = 5;
lamda = 200;
mu = 80.5;
epsilon = 1e-5;
alpha = 1;
T = 0.1;
gamma1 = 0.05;
gamma2 = 0.05;
gamma3 = 0.05;
gamma4 = 0.05;
beta = 10;
sigma = 95;
tau = 1e-5;
nena = 1e-5;
rT = 1024;
L = 200;
m = 500;
n = 600;

% Time vector for plotting
k = 1:m;
t = (k-1) * T;

phi1 = zeros(m+1, 1); phi2 = zeros(m+1, 1); phi3 = zeros(m+1, 1); phi4 = zeros(m+1, 1);
mfa1 = zeros(m+1, 1); mfa2 = zeros(m+1, 1); mfa3 = zeros(m+1, 1); mfa4 = zeros(m+1, 1);
sm1 = zeros(m, 1); sm2 = zeros(m, 1); sm3 = zeros(m, 1); sm4 = zeros(m, 1);
y1 = zeros(m+1, 1); y2 = zeros(m+1, 1); y3 = zeros(m+1, 1); y4 = zeros(m+1, 1);
u1 = zeros(m, 1); u2 = zeros(m, 1); u3 = zeros(m, 1); u4 = zeros(m, 1);
xi1 = zeros(m, 1); xi2 = zeros(m, 1); xi3 = zeros(m, 1); xi4 = zeros(m, 1);
v1 = zeros(m, 1); v2 = zeros(m, 1); v3 = zeros(m, 1); v4 = zeros(m, 1);
s1 = zeros(m, 1); s2 = zeros(m, 1); s3 = zeros(m, 1); s4 = zeros(m, 1);
omega1 = zeros(m, 1); omega2 = zeros(m, 1); omega3 = zeros(m, 1); omega4 = zeros(m, 1);

for k = 1:m
    % Phi updates
    if k == 1
        phi1(k) = 1; phi2(k) = 1; phi3(k) = 1; phi4(k) = 1;
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

    % Error dynamics
    xi1(k) = yd(k) - 2*y1(k) + y4(k);
    xi2(k) = y1(k) - 2*y2(k) + y3(k);
    xi3(k) = y2(k) + yd(k) - 2*y3(k);
    xi4(k) = y1(k) + y3(k) - 2*y4(k);

    if k == 1
        integral_xi1 = 0;
        integral_xi2 = 0;
        integral_xi3 = 0;
        integral_xi4 = 0;
    else
        integral_xi1 = integral_xi1 + T * xi1(k);
        integral_xi2 = integral_xi2 + T * xi2(k);
        integral_xi3 = integral_xi3 + T * xi3(k);
        integral_xi4 = integral_xi4 + T * xi4(k);
    end


     % Sliding surfaces
    if k == 1
        s1(k) = 0; s2(k) = 0; s3(k) = 0; s4(k) = 0;
    else
        s1(k) = xi1(k) + alpha * integral_xi1;
        s1(k) = xi1(k) + alpha * integral_xi1;
    end

end
