% Erivelton Gualter dos Santos

clc, clear all, close all

load('pendulumdata.mat')

g = 9.81;

q_min = min(q);
q_max = max(q);
qrange = q_max - q_min; 

qDot_min = min(qDot);
qDot_max = max(qDot);
qDotrange = qDot_max - qDot_min; 

qDDOT_min = min(qDDOT);
qDDOT_max = max(qDDOT);
dDDotrange = qDDOT_max - qDDOT_min;

display(['    q should be between ',num2str(q_min), ' and ',num2str(q_max)]);
display([' qDot should be between ',num2str(qDot_min), ' and ',num2str(qDot_max)]);
display(['qDDot should be between ',num2str(qDDOT_min), ' and ',num2str(qDDOT_max)]);

ranges = [qrange; qDotrange; dDDotrange]
sort_ranges = sort(ranges,'descend')

M = diag(sort_ranges)./sort_ranges(3)

x = [q qDot qDDOT]; 

phi = [g*cos(q)/2 qDDOT/3];
W = pinv(phi)*T;

y = W'*phi;


function out = get_euclid_dist ( x, q)

    out = sqrt((x-q)'*(x-q))
end

function out = get_gaussian_kernal ( d, h)

    out = exp(-d^2/h)
end