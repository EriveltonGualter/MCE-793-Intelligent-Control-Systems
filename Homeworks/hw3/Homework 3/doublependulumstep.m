%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% doublependulumstep.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Eric Schearer
% Date created: 2/23/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [qOut,qDotOut,qDDotOut] = doublependulumstep(qIn,qDotIn,torqueIn,dt)

% constants
m1 = 2.4;
m2 = 1.3;
l1 = 0.8;
l2 = 0.6;
r1 = l1/2;
r2 = l2/2;
I1 = m1*l1^2/12;
I2 = m2*l2^2/12;
g = 9.81;


qDDotOut = compute_accel2(qIn(1),qIn(2),qDotIn(1),qDotIn(2),torqueIn(1),torqueIn(2),g,l1,l2,m1,m2);
qDotOut  = [qDotIn(1) + qDDotOut(1)*dt    qDotIn(2) + qDDotOut(2)*dt];
qOut     = [qIn(1)    + qDotIn(1)*dt      qIn(2) + qDotIn(2)*dt];
  





    