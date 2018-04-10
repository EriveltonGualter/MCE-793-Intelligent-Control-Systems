%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GUALTERhw3.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: Homework 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Erivelton Gualter
% Date created: 3/9/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all

% a) Locally Weighted Regression %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('doublependulumdata.mat'); % Load  q, qDDOT, qDot, tau

x = [q(:,1) q(:,2) qDot(:,1) qDot(:,2) qDDot(:,1) qDDot(:,2)];
y1 = tau(:,1);
y2 = tau(:,2);

id = 1:10:length(q);
x = x(id,:);
y1 = y1(id);
y2 = y2(id);

% Given Matrix 6x6
M = diag([200^0.5 200^0.5 10^0.5 10^0.5 1 1]);

h=100;

for i=1:length(x)
    qq = x(i,:);
    
    xx = x;
    xx(i,:) = [];
    
    for j=1:length(x)-1

        d(i,j) = ((xx(j,:)-qq)*M'*M*(xx(j,:)-qq)')^0.5;
        s = (exp(-((d(i,j)^2)/h)))^0.5;
        
        Z(j,:) = (s*xx(j,:))';
        V1(j,:) = (s*y1(j))';
        V2(j,:) = (s*y2(j))';
    end

    W_hat1 = pinv(((Z')*Z))*((Z')*V1);
    W_hat2 = pinv(((Z')*Z))*((Z')*V2);

    Y_hat(i,1) = W_hat1'*qq';
    Y_hat(i,2) = W_hat2'*qq';
    
end
error= sqrt(median(median((tau(id,:)-Y_hat).^2)))

LOOCV=error;
display(['LOOCV RMS error for h=100 is ',num2str(LOOCV)])

%a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a22a2a2a
%Plot Actual Torque and Predictios using h=100

figure
plot(Y_hat(:,1),tau10(:,1),'ro')
hold on
plot(Y_hat(:,2),tau10(:,2),'b*')
%legend('Actual Torque q1','Predicted Torque q1','Actual Torque q2','Predicted Torque q2')
title('Part a2. Actual Torque VS Predicted Torque')
ylabel('Actual Torques')
xlabel('Predicted Torques')
hold off
% 
% %bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
%Computed Torque Control

% load('desiredtrajecotry.mat')
% 
% ts=0.001;
% l1=0.8;l2=0.6;
% %Initial Conditions
% qN(1,1)=q1Des(1);
% qDotN(1,1)=q1DotDes(1);
% qDDotN(1,1)=q2DDotDes(1);
% qN(1,2)=q2Des(1);
% qDotN(1,2)=q2DotDes(1);
% qDDotN(1,2)=q2DDotDes(1);
% 
% for i=1:length(q1Des)                                                   
% q_q=[q1Des(i);q2Des(i);q1DotDes(i);q2DotDes(i);q1DDotDes(i);q2DDotDes(i)];
% h=100;
% 
%     for j=1:length(q1Des)                                                  
%         x=[q1Des(j);q2Des(j);q1DotDes(j);q2DotDes(j);q1DDotDes(j);q2DDotDes(j)];
%         d=((x-q_q)'*M'*M*(x-q_q))^0.5;
%         k=exp(-((d^2)/h));
%         s(j)=k^0.5;
%         Z(j,:)=(s(j)*[q1Des(j);q2Des(j);q1DotDes(j);q2DotDes(j);q1DDotDes(j);q2DDotDes(j)])';
%         V(j,:)=(s(j)*tau(j,:))';                           
%     end
%     W_hat=pinv(((Z')*Z))*((Z')*V);
%     Y_hat(i,:)=W_hat'*q_q;
% 
%     %%%
%     S_sum=0;
%     for k=1:length(q1Des)
%        r(k,:)=Z(k,:)*W_hat-V(k,:);
%        S_temp=s(k)*s(k);
%        p_lwr(k)=s(k)*s(k)*Z(k,:)*pinv(Z'*Z)*Z(k,:)';
%        S_sum=S_sum+S_temp;
%     end
%     sigma(i,:)=sqrt(sum(r.*r))/sqrt(abs(S_sum-sum(p_lwr)));
%     %%%
% end
% 
% %bb=((Z'*Z)^-1)*(Z'*S)'
% %var_Y_hat=sigma.^2+sigma.^2*bb'*bb;
% mean_Y_hat=mean(Y_hat);
% standart_Y_hat=std(Y_hat);
% up_95=mean_Y_hat+1.66*standart_Y_hat;
% down_95=mean_Y_hat-1.66*standart_Y_hat;
% 
% for i=1:(length(q1Des)-1)           
%     
% [qN(i+1,:),qDotN(i+1,:),qDDotN(i+1,:)] = doublependulumstep(qN(i,:),qDotN(i,:),Y_hat(i,:),ts);
% 
% end
% 
% %b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1b1
% %%B1.Plot of Predicted Torques
% 
% figure
% plot(time,Y_hat(:,1),'b')
% hold on
% plot(time,Y_hat(:,2),'g')
% plot(time,up_95(:,1),'.','MarkerSize',10,'MarkerFaceColor','r')
% plot(time,down_95(:,1),'.','MarkerSize',10,'MarkerFaceColor','r')
% title('Part B1. Plot of Predicted Torques')
% ylabel('N.m')
% xlabel('Time')
% legend('Torque-Predicted q1','Torque-Predicted q2')
% hold off

% %b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2
% %%B2. Plot of Actual and Desired  Joint Angles
% figure
% plot(time,qN(:,1),'b')
% hold on
% plot(time,qN(:,2),'g')
% plot(time,q1Des,'y')
% plot(time,q2Des,'r')
% title('Part B2. Plot of Actual and Desired  Joint Angles')
% ylabel('Rad')
% xlabel('Time')
% legend('Actual q1','Actual q2','Desired q1','Desired q2')
% hold off
% 
% %b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3b3
% %%B3. Animation
% ang=0:0.01:2*pi;
% rc=0.25;
% xc=0;
% yc=1;
% xp=rc*cos(ang);
% yp=rc*sin(ang);
% figure
% for i=2:10:length(time)
%    Px1=l1*cos(qN(i,1));
%    Py1=l1*sin(qN(i,1));
%    Px2=l1*cos(qN(i,1))+l2*cos((qN(i,1))+(qN(i,2)));
%    Py2=l1*sin(qN(i,1))+l2*sin((qN(i,1))+(qN(i,2)));
%    hold off;
%    plot(l1*cos(q1Des(i))+l2*cos(q1Des(i)+q2Des(i)),l1*sin(q1Des(i))+l2*sin(q1Des(i)+q2Des(i)),'o','MarkerSize',10,'MarkerFaceColor','r');
%    title('Part B.3 Animation of the Solution')
%    hold on;   
%    plot(xc+xp,yc+yp,'r');
%    plot(l1*cos(q1Des(i))+l2*cos(q1Des(i)+q2Des(i)),l1*sin(q1Des(i))+l2*sin(q1Des(i)+q2Des(i)),'o','MarkerSize',10,'MarkerFaceColor','r');
%    line([0 Px1],[0 Py1],'LineWidth',4);line([Px1 Px2],[Py1 Py2],'LineWidth',4);
%    axis([-1.5*(l1+l2) 1.5*(l1+l2) -1.5*(l1+l2) 1.5*(l1+l2)]);  
%    drawnow;
% end
% 
% %cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
% %%Computed Torque Control with Feedback
% 
% for i=1:length(q1Des)                                                   
% q_q=[q1Des(i);q2Des(i);q1DotDes(i);q2DotDes(i);q1DDotDes(i);q2DDotDes(i)];
% h=100;
% 
% for j=1:length(q1Des)                                                  
% x=[q1Des(j);q2Des(j);q1DotDes(j);q2DotDes(j);q1DDotDes(j);q2DDotDes(j)];
% d=((x-q_q)'*M'*M*(x-q_q))^0.5;
% k=exp(-((d^2)/h));
% s=k^0.5;
% Z(j,:)=(s*[q1Des(j);q2Des(j);q1DotDes(j);q2DotDes(j);q1DDotDes(j);q2DDotDes(j)])';
% V(j,:)=(s*tau(j,:))';                             %%   CHANGING VALUES
% end
% W_hat=pinv(((Z')*Z))*((Z')*V);
% Y_hat(i,:)=W_hat'*q_q;
% end
% Kp=500;
% Kd=500;
% for i=1:(length(q1Des)-1)                                                    %%%%%%%%%ERRASE
% Y_hat(i,1)=Y_hat(i,1)+Kp*(q1Des(i)-qN(i,1))+Kd*(q1DotDes(i)-qDotN(i,1));
% Y_hat(i,2)=Y_hat(i,2)+Kp*(q2Des(i)-qN(i,2))+Kd*(q2DotDes(i)-qDotN(i,2));
% [qN(i+1,:),qDotN(i+1,:),qDDotN(i+1,:)] = doublependulumstep(qN(i,:),qDotN(i,:),Y_hat(i,:),ts);
% 
% end
% 
% %c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1c1
% %%C1.Plot of Torque Predictions                             %%%%%% Two SD95
% 
% figure
% plot(time,Y_hat(:,1),'b')
% hold on
% plot(time,Y_hat(:,2),'g')
% title('Part C1. Plot of Original and Predicted Torques')
% ylabel('N.m')
% xlabel('Time')
% legend('Torque-Predicted q1','Torque-Predicted q2')
% hold off
% 
% %c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2c2
% %%C2. Plot of Actual and Desired  Joint Angles
% figure
% plot(time,qN(:,1),'b')
% hold on
% plot(time,qN(:,2),'g')
% plot(time,q1Des,'y')
% plot(time,q2Des,'r')
% title('Part C2. Plot of Actual and Desired  Joint Angles')
% ylabel('Rad')
% xlabel('Time')
% legend('Actual q1','Actual q2','Desired q1','Desired q2')
% hold off
% 
% %c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3c3
% %%C3. Animation
% ang=0:0.01:2*pi;
% rc=0.25;
% xc=0;
% yc=1;
% xp=rc*cos(ang);
% yp=rc*sin(ang);
% figure
% for i=2:10:length(time)
%    Px1=l1*cos(qN(i,1));
%    Py1=l1*sin(qN(i,1));
%    Px2=l1*cos(qN(i,1))+l2*cos((qN(i,1))+(qN(i,2)));
%    Py2=l1*sin(qN(i,1))+l2*sin((qN(i,1))+(qN(i,2)));
%    hold off;
%    plot(l1*cos(q1Des(i))+l2*cos(q1Des(i)+q2Des(i)),l1*sin(q1Des(i))+l2*sin(q1Des(i)+q2Des(i)),'o','MarkerSize',10,'MarkerFaceColor','r');
%    title('Part C.3 Animation of the Solution')
%    hold on;   
%    plot(xc+xp,yc+yp,'r');
%    plot(l1*cos(q1Des(i))+l2*cos(q1Des(i)+q2Des(i)),l1*sin(q1Des(i))+l2*sin(q1Des(i)+q2Des(i)),'o','MarkerSize',10,'MarkerFaceColor','r');
%    line([0 Px1],[0 Py1],'LineWidth',4);line([Px1 Px2],[Py1 Py2],'LineWidth',4);
%    axis([-1.5*(l1+l2) 1.5*(l1+l2) -1.5*(l1+l2) 1.5*(l1+l2)]);  
%    drawnow;
% end