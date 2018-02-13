close all
clear all
clc
load('desiredtrajecotry.mat')
load('pendulumdata.mat')

%%Bayesian Regression
m_pri=[1 1];
s_pri=[10 0;0 10];
B=100;

for k=1:6000
    Phi=[9.81*cos(q(k))/2 qDDOT(k)/3];
    s_post(:,:,k)=inv(inv(s_pri)+B*Phi'*Phi);
    m_post(k,:)=s_post(:,:,k)*((pinv(s_pri)*m_pri')+B*Phi'*T(k)');
    s_pri=s_post(:,:,k);
    m_pri=m_post(k,:); 
    s_va(k,1)=s_post(1,1,k)^0.5;
    s_va(k,2)=s_post(2,2,k)^0.5;
end
plot(m_post)
title('Part (a1)')

ylabel('')
xlabel('Time')
legend('W1','W2')
figure
plot(s_va)
title('Part (a2)')
ylabel('')
xlabel('Time')
legend('Sigma1','Sigma2')
Phi=[9.81*cos(q)/2 qDDOT/3];
W=pinv(Phi)*T;
l=W(2)/W(1);
m=W(1)/l;
display(['Part (a3): m = ',num2str(m)])
display(['Part (a3): l = ',num2str(l)])
%%
%%Computed Torque Control
ts=0.001;
%%0
qN(1)=qDes(1)*0.95;
qDotN(1)=qDotDes(1)*0.95;
for k=1:length(time)
TDes(k)=(m*l^2*qDDotDes(k)/3)+(m*9.81*l*cos(qDes(k))/2);
qDDotN(k)=(3*TDes(k)/(m*l^2))-(3*9.81*cos(qN(k))/(2*l));
k=k+1;
qN(k)=ts*qDotN(k-1);
qN(k)=qN(k)+qN(k-1);
qDotN(k)=ts*qDDotN(k-1);
qDotN(k)=qDotN(k)+qDotN(k-1);
k=k-1;
end

%%1
figure
plot(time,qN(1:6001),'b*',time,qDes,'r*')
title('Part (b1) Angles')
ylabel('Rad')
xlabel('Time')
legend('q','qDes')
%%

%%2
error_total=0;
for k=1:length(time)
    error_temp=(qN(k)-qDes(k))^2;
    error_total=error_total+error_temp;
end

rmserror=((error_total)/length(time))^0.5;
display(['Part (b2): rms error = ',num2str(rmserror)])

%%3
figure
for k=1:10:length(time)
   Px=l*cos(qN(k));
   Py=l*sin(qN(k));
   hold off;
   plot(l*cos(qDes(k)),l*sin(qDes(k)),'o','MarkerSize',10,'MarkerFaceColor','r');
   title('Part (b3) Animation')
   hold on;
   line([0 Px],[0 Py],'LineWidth',4);
   axis([-(1.2*l) 1.2*l -(1.2*l) (1.2*l)]);  
   drawnow;
end
%%
%%Computed Torque Control with Feedback
Kp=60;Kd=15;

%%0
qNN(1)=qDes(1)*0.95;
qDotNN(1)=qDotDes(1)*0.95;
for k=1:length(time)
TDes(k)=(m*l^2*qDDotDes(k)/3)+(m*9.81*l*cos(qDes(k))/2)+Kp*(qDes(k)-qNN(k))+Kd*(qDotDes(k)-qDotNN(k));
qDDotNN(k)=(3*TDes(k)/(m*l^2))-(3*9.81*cos(qNN(k))/(2*l));
k=k+1;
qNN(k)=ts*qDotNN(k-1);
qNN(k)=qNN(k)+qNN(k-1);
qDotNN(k)=ts*qDDotNN(k-1);
qDotNN(k)=qDotNN(k)+qDotNN(k-1);
k=k-1;
end

%%1
figure
plot(time,qNN(1:6001),'b*',time,qDes,'r*')
title('Part (c1) Angles')
ylabel('Rad')
xlabel('Time')
legend('q','qDes')

%%2
error_total=0;
for k=1:length(time)
   
    error_temp=(qNN(k)-qDes(k))^2;
    error_total=error_total+error_temp;
    
end

rmserror=((error_total)/length(time))^0.5;
display(['Part (b2): rms error = ',num2str(rmserror)])

%%3
figure
for k=1:10:length(time)
   Px=l*cos(qNN(k));
   Py=l*sin(qNN(k));
   hold off;
   plot(l*cos(qDes(k)),l*sin(qDes(k)),'o','MarkerSize',10,'MarkerFaceColor','r');
   title('Part (c3) Animation')
   hold on;
   line([0 Px],[0 Py],'LineWidth',4);
   axis([-(1.2*l) 1.2*l -(1.2*l) (1.2*l)]);  
   drawnow;
end