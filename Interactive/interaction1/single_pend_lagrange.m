%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% double_pend_lagrange.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: This derives the equations of motion for a double pendulum.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Eric Schearer
% Date created: 1/25/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

% define various constants
syms l r m I g real
% define sysbolic variables
syms t real
q = sym('q(t)');
syms Q Qdot Qddot T real

% position of CoMs
x = r*cos(q);
y = r*sin(q);

% velocity of CoMs
x1dot = diff(x,t);
y1dot = diff(y,t);
%%
% compute potential energy
U = m*g*y;
Usave = subs(U,{q,diff(q,t)},{Q,Qdot})

% compute kinetic energy
Ttrans = simplify(0.5*m*(x1dot^2 + y1dot^2)); 
Trot   = 0.5*I*diff(q,t)^2; 
T = Ttrans + Trot;
Tsave = subs(T,{q,diff(q,t)},{Q,Qdot})

% compute the Lagrangian
L = simplify(T - U);
% substitute q1dot and q2dot so Matlab can do partial derivatives
L = subs(L,{q,diff(q,t)},{Q,Qdot});
% compute derivatives for Euler-Lagrange equations
dLdQdot    = diff(L,Qdot);
dLdQ       = diff(L,Q);
% substitute back in so Matlab can do time derivatives
dLdQdot    = subs(dLdQdot,{Q,Qdot},{q,diff(q,t)});
dLdQ       = subs(dLdQ   ,{Q,Qdot},{q,diff(q,t)});

% now compute the Euler-Lagrange equations
eqn = diff(dLdQdot,t) - dLdQ - T;
% substitute in position, velocity, acceleration
eqn = subs(eqn,{q,diff(q,t),diff(q,t,t)},{Q,Qdot,Qddot});

% solve for the accelerations
accel = solve(eqn,Qddot);

% make these symbolic expressions into Matlab functions
matlabFunction(Usave,'file','compute_potential')
matlabFunction(Tsave,'file','compute_kinetic')
matlabFunction(accel,'file','compute_accel')