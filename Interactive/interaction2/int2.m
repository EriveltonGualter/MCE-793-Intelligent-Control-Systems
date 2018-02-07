% Erivelton Gualter dos Santos

clear all, close all

load('pendulumdata.mat'); % Load data 

g= 9.81;   % gravitational acceleration

t = [qDDOT];
phi = [3*T -3*g*cos(q)/2];

Wml = inv(phi'*phi)*phi'*t;