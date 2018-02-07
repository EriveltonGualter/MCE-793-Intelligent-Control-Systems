% Intelligent Control System – Interactive Session Dynamic Simulations
% 01/24/2017
% Erivelton Gualter

clear all, close all

% Parameters
m = 1;      % mass rod [kg]
l = 1;      % length rod [m]
g = -9.81;   % gravitacional acceleration [m/s2]

I = m*l^2/12;
r = l/2;

% Simulation parameters
tf = 5;     % Final time [s]
dt = 1e-3;  % Integration step size [s]
t = 0:dt:tf; % array time

% Initial Conditional
x10 = pi/4;
x20 = 0;

X = [x10; x20];

for i=1 : length(t)-1
    
    x1 = X(1,i);
    x2 = X(2,i);  
    
    
    xd1 = x2;
    xd2 = -m*g*r*sin(x1) / (m*r^2+I);    
    
    XDOT = [xd1; xd2];
    
    X(:,i+1) = X(:,i) + XDOT* dt;
    
    T(i) = x2^2 * (m*r^2 + I) / 2;
    U(i) = m*g*(r-r*cos(x1));
end

mySinglePendulumAnimation( t, X(1,:), l );

subplot(211); plot(t,X, 'LineWidth', 2); legend('Angle','Angular Velocity');
title('Simulation'); xlabel('Time (s)'); ylabel('Angle');

t = t(2:end); 
subplot(212); plot(t,T,t,U, t,T+U, 'LineWidth', 2); 
legend('Kinetic Energy','Potential Energy','Total Energy');
xlabel('Time (s)'); ylabel('Angular Velocity');
