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
syms l1 l2 r1 r2 m1 m2 I1 I2 g real
% define sysbolic variables
syms t real
q1 = sym('q1(t)');
q2 = sym('q2(t)');
syms Q1 Q2 Q1dot Q2dot Q1ddot Q2ddot T1 T2 real

% position of CoMs
x1 = r1*cos(q1);
y1 = r1*sin(q1);
x2 = l1*cos(q1) + r2*cos(q1+q2);
y2 = l1*sin(q1) + r2*sin(q1+q2);

% velocity of CoMs
x1dot = diff(x1,t);
y1dot = diff(y1,t);
x2dot = diff(x2,t);
y2dot = diff(y2,t);

% compute potential energy
U = m1*g*y1 + m2*g*y2;
Usave = subs(U,{q1,q2,diff(q1,t),diff(q2,t)},{Q1,Q2,Q1dot,Q2dot});

% compute kinetic energy
T1trans = simplify(0.5*m1*(x1dot^2 + y1dot^2)); 
T1rot   = 0.5*I1*diff(q1,t)^2; 
T2trans = simplify(0.5*m2*(x2dot^2 + y2dot^2)); 
T2rot   = 0.5*I2*(diff(q1,t) + diff(q2,t))^2;
T = T1trans + T1rot + T2trans + T2rot;
Tsave = subs(T,{q1,q2,diff(q1,t),diff(q2,t)},{Q1,Q2,Q1dot,Q2dot});

% compute the Lagrangian
L = simplify(T - U);
% substitute q1dot and q2dot so Matlab can do partial derivatives
L = subs(L,{q1,q2,diff(q1,t),diff(q2,t)},{Q1,Q2,Q1dot,Q2dot});
% compute derivatives for Euler-Lagrange equations
dLdQ1dot    = diff(L,Q1dot);
dLdQ1       = diff(L,Q1);
dLdQ2dot    = diff(L,Q2dot);
dLdQ2       = diff(L,Q2);
% substitute back in so Matlab can do time derivatives
dLdQ1dot    = subs(dLdQ1dot,{Q1,Q2,Q1dot,Q2dot},{q1,q2,diff(q1,t),diff(q2,t)});
dLdQ1       = subs(dLdQ1   ,{Q1,Q2,Q1dot,Q2dot},{q1,q2,diff(q1,t),diff(q2,t)});
dLdQ2dot    = subs(dLdQ2dot,{Q1,Q2,Q1dot,Q2dot},{q1,q2,diff(q1,t),diff(q2,t)});
dLdQ2       = subs(dLdQ2   ,{Q1,Q2,Q1dot,Q2dot},{q1,q2,diff(q1,t),diff(q2,t)});

% now compute the Euler-Lagrange equations
eqn1 = diff(dLdQ1dot,t) - dLdQ1 - T1;
eqn2 = diff(dLdQ2dot,t) - dLdQ2 - T2;
% substitute in position, velocity, acceleration
eqn1 = subs(eqn1,{q1,q2,diff(q1,t),diff(q2,t),diff(q1,t,t),diff(q2,t,t)},{Q1,Q2,Q1dot,Q2dot,Q1ddot,Q2ddot});
eqn2 = subs(eqn2,{q1,q2,diff(q1,t),diff(q2,t),diff(q1,t,t),diff(q2,t,t)},{Q1,Q2,Q1dot,Q2dot,Q1ddot,Q2ddot});

% solve for the accelerations
accel = solve(eqn1,eqn2,Q1ddot,Q2ddot);
accel = [accel.Q1ddot accel.Q2ddot];

% make these symbolic expressions into Matlab functions
matlabFunction(Usave,'file','compute_potential')
matlabFunction(Tsave,'file','compute_kinetic')
matlabFunction(accel,'file','compute_accel')