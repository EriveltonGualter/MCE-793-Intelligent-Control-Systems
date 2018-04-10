
clear all
clc

load('pendulumdata.mat')


maxq=max(q);
minq=min(q);
maxqDot=max(qDot);
minqDot=min(qDot);
maxqDDot=max(qDDOT);
minqDDot=min(qDDOT);

m1=1/(maxq-minq);
m2=1/(maxqDot-minqDot);
m3=1/(maxqDDot-minqDDot);

M=diag([1 m1 m2 m3]);