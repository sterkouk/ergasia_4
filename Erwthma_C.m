clear all
close all
load('speakerA.mat');
load('speakerB.mat');
% Initialization
n = length(u);
factors = 6600; 
% Wiener
P = xcorr(d, u);
P = P/length(u);
P = P(n:n+factors-1);
r = var(u)*autocorr(u, factors - 1); 
R = toeplitz(r);
wo = R \ P; 
% calculations
y = zeros(n, 1);
for i = factors + 1:n
     y(i) = wo' * u(i:-1:(i - factors + 1));
end
e = d - y; 
sound(e, fs);
% LMS
% initialization
mu = 0.0009;
w = zeros(factors, 1);
e = zeros(n, 1);
y = zeros(n, 1);
% calculations
for i = (factors + 1):n
    y(i) = w' * u(i: - 1:(i - factors + 1));
    e(i) = d(i) - y(i);
    w = w + mu * e(i) * u(i:-1:(i - factors + 1));
end
sound(e, fs);
% normalized LMS
mu = 0.5; 
a = 100;
w = zeros(factors, 1);
e = zeros(n, 1);
y = zeros(n, 1);
for i = (factors + 1):n
    y(i) = w'*u(i: - 1:(i - factors + 1));
    e(i) = d(i) - y(i);
    w = w + mu * e(i) * u(i:-1:(i - factors + 1)) / (a + u(i:-1:(i - factors + 1))' * u(i:-1:(i - factors + 1)));
end
sound(e, fs);
% RLS
l = 1;
de = 0.005;
P = (1 / de) * eye(factors, factors);
w = zeros(factors, 1);
e = zeros(n, 1);
y = zeros(n, 1);
%factors = 500;
% calculations
for i = (factors + 1):n
    y(i) = w' * u(i:-1:(i - factors + 1));
    k = ((l^-1) * P * u(i:-1:i - factors + 1) / (1 + (l^-1) * u(i:-1:i - factors + 1)' * P * u(i:-1:(i - factors + 1))));
    e(i) = d(i) - y(i);
    w = w + k * e(i);
    P = (l^-1) *P - (l^-1) * k * u(i:-1:(i - factors + 1))' * P;
end
sound(e, fs)