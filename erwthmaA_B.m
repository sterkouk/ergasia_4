clear all
close all
% Initialiaztion
n = 5000; 
factors = 100; 
% noise v1
s2d = 0.42;
v1 = sqrt(s2d) * randn(n,1);
v1 = v1 - mean(v1);
% noise v2
s2d = 0.72;
v2 = sqrt(s2d) * randn(n,1);
v2 = v2 - mean(v2);
% signal u
u = zeros(n,1);
for i = 1:3
   u(i) = v1(i); 
end
for i=4:n
    u(i) = -0.87 * u(i - 1) - 0.22 * u(i - 2) - 0.032 * u(i - 3) + v1(i);
end
% signal s
s = zeros(n, 1);
s(1) = -0.13 * u(1);
s(2) = 0.67 * u(1) - 0.13 * u(2);
s(3) = -0.18 * u(1) + 0.67 * u(2) - 0.13 * u(3);%First manually create the 3 values
for i = 4:n
   s(i) = -0.13 * u(i) + 0.67 * u(i-1) - 0.18 * u(i-2) - 0.39 * u(i-3);
end
% signal x
x = zeros(n, 1);
x(1) = v2(1);
x(2) = -0.57 * x(1) + v2(2);
x(3) = -0.57 * x(2) - 0.16 * x(1) +  v2(3);
for i = 4:n
    x(i) = -0.57 * x(i-1) - 0.16 * x(i-2) - 0.08 * x(i-3) + v2(i);
end
% signal d
d = s + x;
% Calculate P, R, wo of the Wiener Filter 
r = var(u)*autocorr(u, factors - 1); 
R = toeplitz(r);
P = xcorr(d, u);
P = P/length(u);
P = P(n:n + factors - 1);
wo = R \ P; 
tic
% calc y
y = zeros(n, 1);
for i = factors + 1:n
     y(i) = wo' * u(i:-1:(i - factors + 1));
end
% calc error 
e = d - y;
errorWiener = (e - x).^2;
% plot
figure(1)
semilogy(errorWiener)
title('Wiener Error')
ylabel('Error');
xlabel('step (n)');
display('Elapsed times')
fprintf('Wiener: \t\t\t%f\n',toc);
% LMS
% initialization
mu = 0.0009;
wo = zeros(factors, 1);
y = zeros(n, 1);
e = zeros(n, 1);
J_lms = zeros(n, 1);
tic
% calculations
for i = (factors + 1):n
    y(i) = wo' * u(i:-1:(i - factors + 1));
    e(i) = d(i) - y(i);
    wo = wo + mu * e(i) * u(i:-1:(i - factors + 1));
    J_lms(i) = (e(i) - x(i))^2;
end
fprintf('LMS: \t\t\t\t%f\n',toc);

% plot
figure(2)
semilogy(J_lms)
title('LMS Error')
ylabel(' Error ');
xlabel('step (n)');
% normalized LMS
% init
mu = 0.06; 
a = 100;
wo = zeros(factors, 1);
y = zeros(n, 1);
e = zeros(n, 1);
J_nlms = zeros(n, 1);
tic
% calculations
for i = (factors + 1):n
    y(i) = wo' * u(i:-1:(i - factors + 1));
    e(i) = d(i) - y(i);
    wo = wo + mu * e(i) * u(i:-1:(i - factors + 1)) / (a + u(i:-1:(i - factors + 1))' * u(i:-1:(i - factors + 1)));
    J_nlms(i) = (e(i) - x(i))^2;
end
fprintf('normalized LMS: \t%f\n',toc);
% plot
figure(3)
semilogy(J_nlms)
title('Normalized LMS Error')
ylabel('Square Error ( E[e^2(n)] )');
xlabel('step (n)');
% RLS
% initialization
l = 1; % lamda
de = 0.005; % delta
y = zeros(n, 1);
e = zeros(n, 1);
J_rls = zeros(n, 1);
wo = zeros(factors, 1);
P = (1 / de) * eye(factors, factors);
tic
% calculations
for i = (factors + 1):n
    y(i) = wo' * u(i:-1:(i - factors + 1));
    k = ((l^-1) * P * u(i:-1:i - factors + 1) / (1 + (l^-1) * u(i:-1:i - factors + 1)'* P * u(i:-1:(i - factors + 1))));
    e(i) = d(i) - y(i);
    wo = wo + k * e(i);
    P = (l^-1) * P - (l^-1) * k * u(i:-1:(i - factors + 1))' * P;
    J_rls(i) = (e(i) - x(i))^2;
end
fprintf('RLS: \t\t\t\t%f\n',toc);
% plot
figure(4)
semilogy(J_rls)
title('RLS Error')
ylabel('Error');
xlabel('step (n)');

%% plot All

figure(5)
title('Square Mean Error')
semilogy(1:n, errorWiener, 1:n, J_lms, 1:n, J_nlms, 1:n, J_rls)
legend('Wiener', 'LMS', 'normalized LMS', 'RLS');
ylabel('Error');
xlabel('step (n)');