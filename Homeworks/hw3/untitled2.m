% close all
clear all
clc

load('desiredtrajecotry.mat')
load('doublependulumdata.mat')

m=zeros(6,length(q));

p1=20;p2=50;
for i=1:length(q)
            
    x1=[q(i,1);q(i,2);qDot(i,1);qDot(i,2);qDDot(i,1);qDDot(i,2)];
   
    for j=1:length(q)
            
    x2=[q(j,1);q(j,2);qDot(j,1);qDot(j,2);qDDot(j,1);qDDot(j,2)];

    K(i,j)=p1*exp(-(x1-x2)'*(x1-x2)/(2*p2^2));
        
    end
end
for i=1:length(q)
    x1=[q(i,1);q(i,2);qDot(i,1);qDot(i,2);qDDot(i,1);qDDot(i,2)];
            q10=q(:,1);
            q10(i)=[];
            qDot10=qDot(:,1);
            qDot10(i)=[];
            qDDot10=qDDot(:,1);
            qDDot10(i)=[];
            q20=q(:,2);
            q20(i)=[];
            qDot20=qDot(:,2);
            qDot20(i)=[];
            qDDot20=qDDot(:,2);
            qDDot20(i)=[];
            K0=K;
            K0(:,i)=[];
            K0(i,:)=[];
            T0=tau;
            T0(i,:)=[];
   
    for j=1:length(q)-1
 
            x2=[q10(j);q20(j);qDot10(j);qDot20(j);qDDot10(j);qDDot20(j)];
            KS(1,j)=p1*exp(-(x1-x2)'*(x1-x2)/(2*p2^2));
            
    end
             Y_Hat(i,1)=KS*inv(K0(1:500,1:500)+1*eye(size(K0)))*T0(:,1);
             Y_Hat(i,2)=KS*inv(K0(1:500,1:500)+1*eye(size(K0)))*T0(:,2);
end
error_p50=sqrt((sum(tau-Y_Hat).^2)/length(tau));
p1=20;p2=100;
for i=1:length(q)
            
    x1=[q(i,1);q(i,2);qDot(i,1);qDot(i,2);qDDot(i,1);qDDot(i,2)];
   
    for j=1:length(q)
            
    x2=[q(j,1);q(j,2);qDot(j,1);qDot(j,2);qDDot(j,1);qDDot(j,2)];

    K(i,j)=p1*exp(-(x1-x2)'*(x1-x2)/(2*p2^2));
        
    end
end
for i=1:length(q)
    x1=[q(i,1);q(i,2);qDot(i,1);qDot(i,2);qDDot(i,1);qDDot(i,2)];
            q10=q(:,1);
            q10(i)=[];
            qDot10=qDot(:,1);
            qDot10(i)=[];
            qDDot10=qDDot(:,1);
            qDDot10(i)=[];
            q20=q(:,2);
            q20(i)=[];
            qDot20=qDot(:,2);
            qDot20(i)=[];
            qDDot20=qDDot(:,2);
            qDDot20(i)=[];
            K0=K;
            K0(:,i)=[];
            K0(i,:)=[];
            T0=tau;
            T0(i,:)=[];
   
    for j=1:length(q)-1
 
            x2=[q10(j);q20(j);qDot10(j);qDot20(j);qDDot10(j);qDDot20(j)];
            KS(1,j)=p1*exp(-(x1-x2)'*(x1-x2)/(2*p2^2));
            
    end
             Y_Hat(i,1)=KS*inv(K0(1:500,1:500)+1*eye(size(K0)))*T0(:,1);
             Y_Hat(i,2)=KS*inv(K0(1:500,1:500)+1*eye(size(K0)))*T0(:,2);
end
error_p100=sqrt((sum(tau-Y_Hat).^2)/length(tau));
p1=20;p2=30;
for i=1:length(q)
            
    x1=[q(i,1);q(i,2);qDot(i,1);qDot(i,2);qDDot(i,1);qDDot(i,2)];
   
    for j=1:length(q)
            
    x2=[q(j,1);q(j,2);qDot(j,1);qDot(j,2);qDDot(j,1);qDDot(j,2)];

    K(i,j)=p1*exp(-(x1-x2)'*(x1-x2)/(2*p2^2));
        
    end
end
for i=1:length(q)
    x1=[q(i,1);q(i,2);qDot(i,1);qDot(i,2);qDDot(i,1);qDDot(i,2)];
            q10=q(:,1);
            q10(i)=[];
            qDot10=qDot(:,1);
            qDot10(i)=[];
            qDDot10=qDDot(:,1);
            qDDot10(i)=[];
            q20=q(:,2);
            q20(i)=[];
            qDot20=qDot(:,2);
            qDot20(i)=[];
            qDDot20=qDDot(:,2);
            qDDot20(i)=[];
            K0=K;
            K0(:,i)=[];
            K0(i,:)=[];
            T0=tau;
            T0(i,:)=[];
   
    for j=1:length(q)-1
 
            x2=[q10(j);q20(j);qDot10(j);qDot20(j);qDDot10(j);qDDot20(j)];
            KS(1,j)=p1*exp(-(x1-x2)'*(x1-x2)/(2*p2^2));
            
    end
             Y_Hat(i,1)=KS*inv(K0(1:500,1:500)+1*eye(size(K0)))*T0(:,1);
             Y_Hat(i,2)=KS*inv(K0(1:500,1:500)+1*eye(size(K0)))*T0(:,2);
end
error_p30=sqrt((sum(tau-Y_Hat).^2)/length(tau));
LOOCV=(error_p30(1,1)+error_p30(1,2))/2;
display(['LOOCV RMS error for p2=30 is ',num2str(LOOCV)])
LOOCV=(error_p50(1,1)+error_p50(1,2))/2;
display(['LOOCV RMS error for p2=50 is ',num2str(LOOCV)])
LOOCV=(error_p100(1,1)+error_p100(1,2))/2;
display(['LOOCV RMS error for p2=100 is ',num2str(LOOCV)])

%a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a22a2a2a
%Plot Actual Torque and Predictios using p2=30

figure
plot(Y_Hat(:,1),tau(:,1),'ro')
hold on
plot(Y_Hat(:,2),tau(:,2),'b*')
title('Part a2. Actual Torque VS Predicted Torque')
ylabel('Actual Torques')
xlabel('Predicted Torques')
hold off

%bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
%%Computed Torque Control

ts=0.001;
l1=0.8;l2=0.6;
%Initial Conditions
qN(1,1)=q1Des(1);
qDotN(1,1)=q1DotDes(1);
qDDotN(1,1)=q2DDotDes(1);
qN(1,2)=q2Des(1);
qDotN(1,2)=q2DotDes(1);
qDDotN(1,2)=q2DDotDes(1);

for i=1:length(q)
            
    x1=[q(i,1);q(i,2);qDot(i,1);qDot(i,2);qDDot(i,1);qDDot(i,2)];
   
    for j=1:length(q)
            
    x2=[q(j,1);q(j,2);qDot(j,1);qDot(j,2);qDDot(j,1);qDDot(j,2)];

    K(i,j)=p1*exp(-(x1-x2)'*(x1-x2)/(2*p2^2));
        
    end
end
for i=1:length(q1Des)                                                   
x1=[q1Des(i);q2Des(i);q1DotDes(i);q2DotDes(i);q1DDotDes(i);q2DDotDes(i)];

    for j=1:length(q)                                                  
    
    x2=[q(j,1);q(j,2);qDot(j,1);qDot(j,2);qDDot(j,1);qDDot(j,2)];
    KS(i,j)=p1*exp(-(x1-x2)'*(x1-x2)/(2*p2^2));

    end
Y_Hat(i,1)=KS(i,:)*inv(K(:,:)+1*eye(size(K)))*tau(:,1);
Y_Hat(i,2)=KS(i,:)*inv(K(:,:)+1*eye(size(K)))*tau(:,2); 
end
for i=1:length(q1Des)
            
    x1=[q1Des(i);q2Des(i);q1DotDes(i);q2DotDes(i);q1DDotDes(i);q2DDotDes(i)];
   
    for j=1:length(q1Des)
            
    x2=[q1Des(j);q2Des(j);q1DotDes(j);q2DotDes(j);q1DDotDes(j);q2DDotDes(j)];
    KSS(i,j)=p1*exp(-(x1-x2)'*(x1-x2)/(2*p2^2));
        
    end
end
for i=1:length(q)                                                   
x1=[q(i,1);q(i,2);qDot(i,1);qDot(i,2);qDDot(i,1);qDDot(i,2)];

    for j=1:length(q1Des)                                                  

    x2=[q1Des(j);q2Des(j);q1DotDes(j);q2DotDes(j);q1DDotDes(j);q2DDotDes(j)];
    KS2(i,j)=p1*exp(-(x1-x2)'*(x1-x2)/(2*p2^2));
    
    end
end

standart_Y_Hat1=std(Y_Hat(:,1));
standart_Y_Hat2=std(Y_Hat(:,2));
for i=1:length(Y_Hat)
up_95_1(i)=Y_Hat(i,1)+1.66*standart_Y_Hat1;
down_95_1(i)=Y_Hat(i,1)-1.66*standart_Y_Hat1;
up_95_2(i)=Y_Hat(i,2)+1.66*standart_Y_Hat2;
down_95_2(i)=Y_Hat(i,2)-1.66*standart_Y_Hat2;
end

%%Covariance of Y Predicted
Cov_Y_Hat=KSS-KS*inv(K+1*eye(size(K)))*KS2;
Cov_Y_Hat=Cov_Y_Hat(:,1).^0.5;

for i=1:(length(q1Des)-1)             
[qN(i+1,:),qDotN(i+1,:),qDDotN(i+1,:)] = doublependulumstep(qN(i,:),qDotN(i,:),Y_Hat(i,:),ts);
end

%b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1
%%B1.Plot of Predicted Torques

figure
plot(time,Y_Hat(:,1),'b','LineWidth',5)
hold on
plot(time,Y_Hat(:,2),'y','LineWidth',5)
plot(time,up_95_1,'r','LineWidth',1)
plot(time,down_95_1,'g','LineWidth',1)
plot(time,up_95_2,'r','LineWidth',1)
plot(time,down_95_2,'g','LineWidth',1)
title('Part B1. Plot of Predicted Torques and 95% Confidence')
ylabel('N.m')
xlabel('Time')
legend('Torque-Predicted q1','Torque-Predicted q2','95% Limit-Up','95% Limit-Down')
hold off

%b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2
%%B2. Plot of Actual and Desired  Joint Angles
figure
plot(time,qN(:,1),'b')
hold on
plot(time,qN(:,2),'g')
plot(time,q1Des,'y')
plot(time,q2Des,'r')
title('Part B2. Plot of Actual and Desired  Joint Angles')
ylabel('Rad')
xlabel('Time')
legend('Actual q1','Actual q2','Desired q1','Desired q2')
hold off

%b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3
%%B3. Animation
ang=0:0.01:2*pi;
rc=0.25;
xc=0;
yc=1;
xp=rc*cos(ang);
yp=rc*sin(ang);
figure
for i=2:10:length(time)
   Px1=l1*cos(qN(i,1));
   Py1=l1*sin(qN(i,1));
   Px2=l1*cos(qN(i,1))+l2*cos((qN(i,1))+(qN(i,2)));
   Py2=l1*sin(qN(i,1))+l2*sin((qN(i,1))+(qN(i,2)));
   hold off;
   plot(l1*cos(q1Des(i))+l2*cos(q1Des(i)+q2Des(i)),l1*sin(q1Des(i))+l2*sin(q1Des(i)+q2Des(i)),'o','MarkerSize',10,'MarkerFaceColor','r');
   title('Part B.3 Animation of the Solution')
   hold on;   
   plot(xc+xp,yc+yp,'r');
   plot(l1*cos(q1Des(i))+l2*cos(q1Des(i)+q2Des(i)),l1*sin(q1Des(i))+l2*sin(q1Des(i)+q2Des(i)),'o','MarkerSize',10,'MarkerFaceColor','r');
   line([0 Px1],[0 Py1],'LineWidth',4);line([Px1 Px2],[Py1 Py2],'LineWidth',4);
   axis([-1.5*(l1+l2) 1.5*(l1+l2) -1.5*(l1+l2) 1.5*(l1+l2)]);  
   drawnow;
end

%cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
%%Computed Torque Control with Feedback

Kp=500;
Kd=500;
for i=1:(length(q1Des)-1)
Y_HatN(i,1)=Y_Hat(i,1)+Kp*(q1Des(i)-qN(i,1))+Kd*(q1DotDes(i)-qDotN(i,1));
Y_HatN(i,2)=Y_Hat(i,2)+Kp*(q2Des(i)-qN(i,2))+Kd*(q2DotDes(i)-qDotN(i,2));
[qN(i+1,:),qDotN(i+1,:),qDDotN(i+1,:)] = doublependulumstep(qN(i,:),qDotN(i,:),Y_HatN(i,:),ts);

end

%c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1
%%C1.Plot of Torque Predictions                             %%%%%% Two SD95
Y_HatN(length(q1Des),:)=Y_Hat(length(q1Des),:);
figure
plot(time,Y_HatN(:,1),'b')
hold on
plot(time,Y_HatN(:,2),'g')
title('Part C1. Plot of Original and Predicted Torques')
ylabel('N.m')
xlabel('Time')
legend('Torque-Predicted q1','Torque-Predicted q2')
hold off

%c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2
%%C2. Plot of Actual and Desired  Joint Angles
figure
plot(time,qN(:,1),'b')
hold on
plot(time,qN(:,2),'g')
plot(time,q1Des,'y')
plot(time,q2Des,'r')
title('Part C2. Plot of Actual and Desired  Joint Angles')
ylabel('Rad')
xlabel('Time')
legend('Actual q1','Actual q2','Desired q1','Desired q2')
hold off

%c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3
%%C3. Animation
ang=0:0.01:2*pi;
rc=0.25;
xc=0;
yc=1;
xp=rc*cos(ang);
yp=rc*sin(ang);
figure
for i=2:10:length(time)
   Px1=l1*cos(qN(i,1));
   Py1=l1*sin(qN(i,1));
   Px2=l1*cos(qN(i,1))+l2*cos((qN(i,1))+(qN(i,2)));
   Py2=l1*sin(qN(i,1))+l2*sin((qN(i,1))+(qN(i,2)));
   hold off;
   plot(l1*cos(q1Des(i))+l2*cos(q1Des(i)+q2Des(i)),l1*sin(q1Des(i))+l2*sin(q1Des(i)+q2Des(i)),'o','MarkerSize',10,'MarkerFaceColor','r');
   title('Part C.3 Animation of the Solution')
   hold on;   
   plot(xc+xp,yc+yp,'r');
   plot(l1*cos(q1Des(i))+l2*cos(q1Des(i)+q2Des(i)),l1*sin(q1Des(i))+l2*sin(q1Des(i)+q2Des(i)),'o','MarkerSize',10,'MarkerFaceColor','r');
   line([0 Px1],[0 Py1],'LineWidth',4);line([Px1 Px2],[Py1 Py2],'LineWidth',4);
   axis([-1.5*(l1+l2) 1.5*(l1+l2) -1.5*(l1+l2) 1.5*(l1+l2)]);  
   drawnow;
end
