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
