close all
clear all
clc

tau=1;
alp_z=4;
bet_z=1;
alp_x=1;
g=1;
Y_0=-5;
Y_0dot=0;
sigma=0.1;
dt=0.01;
t=0:dt:10;

y=Y_0;
z=tau*Y_0dot;
x=1;
c=linspace(0,1,11)';
w=rand(11,1);
for i=1:length(t)
    
f(i) = dmpbasis(x(i),c,sigma,w,g,Y_0);
[ydot,zdot] = dmpderivative(y(i),z(i),g,f(i),alp_z,bet_z,tau);
xdot=-1/tau*alp_x*x(i);
y(i+1)=y(i)+ydot*dt;
z(i+1)=z(i)+zdot*dt;
x(i+1)=x(i)+xdot*dt;
if mod(i,100)==0
dmpPhasePortrait([tau alp_z bet_z g],[-7 7],[-5 5],f(i),1)


pause
end
end
%%
g+f(end)/alp_z/bet_z
figure
plot(t,y(1:end-1))
ylabel('y')
figure
plot(t,z(1:end-1))
ylabel('z')
figure
plot(t,x(1:end-1))
ylabel('x')
figure
plot(t,f)
ylabel('f')