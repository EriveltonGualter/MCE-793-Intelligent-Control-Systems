% Erivelton Gualter
% interactive 7

clear all
close all
clc

% Initial Conditions
tau     = 1;
alpha_z = 1;
beta_z  = 1;
alpha_x = 1;
g       = 1;
Y0      = -5;
c = linspace(0,1,11)';
sigma = 0.1;

W = rand(11,1);

%%

% Simulation parameters
dt = 0.01;
tf = 10;
t = 0:dt:tf;

X = 1;
Y = Y0;
Z = tau*0;https://us.pycon.org/2018/about/volunteers/
for i=1:length(t)    
%     F(i) = 
    F(i) = dmpbasis(X(i), c , sigma, W, g, Y)
    
    [YDOT,ZDOT] = dmpderivative(Y(i), Y(i), g, F(i), alpha_z, beta_z, tau);
    XDOT = -1/tau*alpha_x*X(i);
    
    Y(i+1) = Y(i) + YDOT*dt;
    Z(i+1) = Z(i) + ZDOT*dt;
    X(i+1) = X(i) + XDOT*dt;
end