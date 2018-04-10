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

% Problem 1
% a) Locally Weighted Regression %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('doublependulumdata.mat'); % Load  q, qDDOT, qDot, tau

% id = 1:100:length(q);
% x = [q(id') qDot(id') qDDot(id')];
% y = tau(id');

x = [q(:,1) q(:,2) qDot(:,1) qDot(:,2) qDDot(:,1) qDDot(:,2)];
y1 = tau(:,1);
y2 = tau(:,2);

id = 1:100:length(q);
x = x(id,:);
y1 = y1(id);
y2 = y2(id);

% Given Matrix 6x6
M = diag([sqrt(200), sqrt(200), sqrt(10), sqrt(10), 1, 1]);

% Parameter to tune: (Bandwidth hyperparameter) 
harray = [100 1000 10000];

for hi=1: length(harray)
    h = harray(hi);
    for i=1:length(x)
        qq = x(i,:);
        xq = x;
        xq(i,:) = [];

        for j=1:length(x)-1        

            d = ((xq(j,:)-qq)*M'*M*(xq(j,:)-qq)')^0.5;
            k = exp(-((d^2)/h));
            s = k^0.5;
            Z(j,:)  = (s*x(j,:))';
            V1(j,:) = (s*y1(j))';
            V2(j,:) = (s*y2(j))';
        end

        what1 = pinv(Z'*Z)*Z'*V1;   % what1 = ((Z')*Z)\(Z'*V1);
        what2 = pinv(Z'*Z)*Z'*V2;   % what2 = ((Z')*Z)\(Z'*V2);

        yhat1(i,hi) = what1'*qq';
        yhat2(i,hi) = what2'*qq';
    end
    LOOCV(hi,:) = rms([yhat1(:,hi); yhat2(:,hi)] - [y1; y2]);
%     LOOCV(hi) = sqrt(mean(([yhat1(:,hi); yhat2(:,hi)] - [y1; y2]).^2));
end

%%
display(['LOOCV RMS error for h = 100 is ',num2str(LOOCV(1))]);
display(['LOOCV RMS error for h = 1000 is ',num2str(LOOCV(2))]);
display(['LOOCV RMS error for h = 10000 is ',num2str(LOOCV(3))]);
display(['Processing simulation ...']);

[ml,il] = min(LOOCV);

figure; hold on;
plot(yhat1(:,hi), y1, 'o'); plot(yhat2(:,hi), y2, '*');
title('Part a2','Fontsize',14,'interpreter','latex');
ylabel('Actual Torques','Fontsize',14,'interpreter','latex');
xlabel('Predicted Torques','Fontsize',14,'interpreter','latex');
axis equal

%% b) Computed Torque Control %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('desiredtrajecotry.mat');
load('doublependulumdata.mat'); % Load  q, qDDOT, qDot, tau

% Resampling
% id = 1:10:length(q1Des);
% 
% q1DDotDes = q1DDotDes(id);
% q1Des = q1Des(id);
% q1DotDes = q1DotDes(id);
% q2DDotDes = q2DDotDes(id);
% q2Des = q2Des(id);
% q2DotDes = q2DotDes(id);
% time = time(id);
% xDes = xDes(id);
% yDes = yDes(id);

% Given Matrix 6x6
M = diag([sqrt(200), sqrt(200), sqrt(10), sqrt(10), 1, 1]);

dt = 0.001;

%Initial Conditions
qN(1,1) = q1Des(1);
qDotN(1,1) = q1DotDes(1);
qDDotN(1,1) = q2DDotDes(1);
qN(1,2) = q2Des(1);
qDotN(1,2) = q2DotDes(1);
qDDotN(1,2) = q2DDotDes(1);

x = [q1Des q2Des q1DotDes q2DotDes q1DDotDes q2DDotDes];

h = 100;

for i=1:length(x)                                                   
    qq = x(i,:)';

    for j=1:length(x)                                                  
        xx = x(j,:)';
        
        d = ((xx-qq)'*M'*M*(xx-qq))^0.5;
        k = exp(-((d^2)/h));
        s(j) = k^0.5;
        Z(j,:) = (s(j)*xx)';
        V(j,:) = (s(j)*tau(j,:))';                           
    end
    
    what = pinv(((Z')*Z))*((Z')*V);
    yhat(i,:) = what'*qq;

end

yhat_mean   = mean(yhat);
yhat_std    = std(yhat);
top         = yhat_mean + yhat_std;
botton      = yhat_mean - yhat_std;

for i=1 : length(x)-1
    [qN(i+1,:),qDotN(i+1,:),qDDotN(i+1,:)] = doublependulumstep(qN(i,:),qDotN(i,:),yhat(i,:),dt);
end


%% b) Computed Torque Control %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; hold on
plot(time,yhat(:,1));
plot(time,yhat(:,2));
plot(time,top(:,1),'b.');
plot(time,botton(:,1),'b.');
title('Part ($B_1$)','Fontsize',14,'interpreter','latex');
ylabel('N.m','Fontsize',14,'interpreter','latex');
xlabel('Time','Fontsize',14,'interpreter','latex');
legend('q1','q2');
xlim([time(1) time(end)]);

figure; hold on
plot(time,qN(:,1),'LineWidth',2);
plot(time,qN(:,2),'LineWidth',2);
plot(time,q1Des,'LineWidth',2);
plot(time,q2Des,'LineWidth',2);
xlim([time(1) time(end)]);
title('Part $B_2$','Fontsize',14,'interpreter','latex');
ylabel('Rad','Fontsize',14,'interpreter','latex');
xlabel('Time','Fontsize',14,'interpreter','latex');
legend('actual-q1','actual-q2','desired-q1','desired-q2');

%% Animation - Computed Torque Control
clear global
figure;
title('Computed Torque','Fontsize',14,'interpreter','latex');


l1 = 0.8;
l2 = 0.6;

q1 = qN(:,1);
q2 = qN(:,2);

x1 = l1*cos(q1);
y1 = l1*sin(q1);
x2 = l1*cos(q1) + l2*cos(q1+q2);
y2 = l1*sin(q1) + l2*sin(q1+q2);
    
pos = [x1'; y1'; x2'; y2'];

for i=1:length(time)
    px1 = pos(1,i);
    py1 = pos(2,i);
    px2 = pos(3,i);
    py2 = pos(4,i);
    
    draw(px1, py1, px2, py2);
end 

%% c) Computed Torque Control with Feedback %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kp = 600;
Kd = 50;

for i=1:(length(q1Des)-1)
    yhat(i,1) = yhat(i,1) + Kp*(q1Des(i)-qN(i,1)) + Kd*(q1DotDes(i)-qDotN(i,1));
    yhat(i,2) = yhat(i,2) + Kp*(q2Des(i)-qN(i,2)) + Kd*(q2DotDes(i)-qDotN(i,2));
    [qN(i+1,:),qDotN(i+1,:),qDDotN(i+1,:)] = doublependulumstep(qN(i,:),qDotN(i,:),yhat(i,:),dt);
end

%% Animation - Computed Torque Control with Feedback
clear global
figure;
title('Computed Torque Control with Feedback','Fontsize',14,'interpreter','latex');

l1 = 0.8;
l2 = 0.6;

q1 = q1Des;
q2 = q2Des;

x1 = l1*cos(q1);
y1 = l1*sin(q1);
x2 = l1*cos(q1) + l2*cos(q1+q2);
y2 = l1*sin(q1) + l2*sin(q1+q2);
    
pos = [x1'; y1'; x2'; y2'];

for i=1:length(time)
    px1 = pos(1,i);
    py1 = pos(2,i);
    px2 = pos(3,i);
    py2 = pos(4,i);
    
    draw(px1, py1, px2, py2);

end

%% c1
figure; hold on;
plot(time, yhat(:,1),'LineWidth',2);
plot(time, yhat(:,2),'LineWidth',2);
title('Part $C_2$','Fontsize',14,'interpreter','latex');
ylabel('Rad','Fontsize',14,'interpreter','latex');
xlabel('Time','Fontsize',14,'interpreter','latex');
legend('actual-q1','actual-q2','desired-q1','desired-q2');

%% c3
figure; hold on;
plot(time, qN(:,1),'LineWidth',2);
plot(time, qN(:,2),'LineWidth',2);
plot(time, q1Des,'LineWidth',2);
plot(time, q2Des,'LineWidth',2);
title('Part $C_3$','Fontsize',14,'interpreter','latex');
ylabel('Rad','Fontsize',14,'interpreter','latex');
xlabel('Time','Fontsize',14,'interpreter','latex');
legend('actual-q1','actual-q2','desired-q1','desired-q2');

display(['Done Part 1']);

% PROBLEM 2 --------------------------------------------------------------
%% a) Gaussian process regression %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

load('desiredtrajecotry.mat')
load('doublependulumdata.mat')

id = 1:1:501;
q       = q(id,:);
qDDot   = qDDot(id,:);
qDot    = qDot(id,:);
tau     = tau(id,:);

p1 = 20;
p2a = [50 100 30];

for k=1:length(p2a)
    p2 = p2a(k);
    
    for i=1:length(q)
        x1 = [q(i,1); q(i,2); qDot(i,1); qDot(i,2); qDDot(i,1); qDDot(i,2)];
        for j=1:length(q)
            x2 = [q(j,1); q(j,2); qDot(j,1); qDot(j,2); qDDot(j,1); qDDot(j,2)];
            K(i,j)=p1*exp(-(x1-x2)'*(x1-x2)/(2*p2^2));
        end
    end
    for i=1:length(q)
        x1 = [q(i,1); q(i,2); qDot(i,1); qDot(i,2); qDDot(i,1); qDDot(i,2)];
        q10     = q(:,1);           q10(i) = [];
        qDot10  = qDot(:,1);     qDot10(i) = [];
        qDDot10 = qDDot(:,1);   qDDot10(i) = [];
        q20     = q(:,2);           q20(i) = [];
        qDot20  = qDot(:,2);     qDot20(i) = [];
        qDDot20 = qDDot(:,2);   qDDot20(i) = [];

        K0 = K;
        K0(:,i) = [];
        K0(i,:) = [];
        T0 = tau;
        T0(i,:) = [];
        
        for j=1:length(q)-1
            x2 = [q10(j); q20(j); qDot10(j); qDot20(j); qDDot10(j); qDDot20(j)];
            KS(1,j)=p1*exp(-(x1-x2)'*(x1-x2)/(2*p2^2));
        end
        yhat(i,1) = KS*inv(K0+1*eye(size(K0)))*T0(:,1);
        yhat(i,2) = KS*inv(K0+1*eye(size(K0)))*T0(:,2);
    end

    switch p2
        case 30
            error_p30 = sqrt((sum(tau-yhat).^2)/length(tau));
        case 50
            error_p50 = sqrt((sum(tau-yhat).^2)/length(tau));
        case 100
            error_p100 = sqrt((sum(tau-yhat).^2)/length(tau));
    end
end
LOOCV = (error_p30(1,1)+error_p30(1,2))/2;
display(['LOOCV RMS error for p2 = 30 is ',num2str(LOOCV)])
LOOCV = (error_p50(1,1)+error_p50(1,2))/2;
display(['LOOCV RMS error for p2 = 50 is ',num2str(LOOCV)])
LOOCV = (error_p100(1,1)+error_p100(1,2))/2;
display(['LOOCV RMS error for p2 = 100 is ',num2str(LOOCV)])

display(['Processing simulation ...']);
% Plot Actual Torque and Predictios for p2 = 30

figure; hold on
plot(yhat(:,1),tau(:,1),'ro');
plot(yhat(:,2),tau(:,2),'b*');
title('Part a2. Actual Torque VS Predicted Torque')
ylabel('Actual Torques','Fontsize',14,'interpreter','latex');
xlabel('Predicted Torques','Fontsize',14,'interpreter','latex');

%%Computed Torque Control

ts = 0.001;
l1 = 0.8;
l2 = 0.6;

%Initial Conditions
qN(1,1)     = q1Des(1);
qDotN(1,1)  = q1DotDes(1);
qDDotN(1,1) = q2DDotDes(1);
qN(1,2)     = q2Des(1);
qDotN(1,2)  = q2DotDes(1);
qDDotN(1,2) = q2DDotDes(1);

for i=1:length(q)
    x1 = [q(i,1);q(i,2);qDot(i,1);qDot(i,2);qDDot(i,1);qDDot(i,2)];
    
    for j=1:length(q)
        x2 = [q(j,1);q(j,2);qDot(j,1);qDot(j,2);qDDot(j,1);qDDot(j,2)];
        K(i,j) = p1*exp(-(x1-x2)'*(x1-x2)/(2*p2^2));
    end
end

for i=1:length(q1Des)                                                   
    x1 = [q1Des(i);q2Des(i);q1DotDes(i);q2DotDes(i);q1DDotDes(i);q2DDotDes(i)];

    for j=1:length(q)                                                  
        x2 = [q(j,1);q(j,2);qDot(j,1);qDot(j,2);qDDot(j,1);qDDot(j,2)];
        KS(i,j) = p1*exp(-(x1-x2)'*(x1-x2)/(2*p2^2));
    end
    yhat(i,1) = KS(i,:)*inv(K(:,:)+1*eye(size(K)))*tau(:,1);
    yhat(i,2) = KS(i,:)*inv(K(:,:)+1*eye(size(K)))*tau(:,2); 
end

for i=1:length(q1Des)
            
    x1 = [q1Des(i);q2Des(i);q1DotDes(i);q2DotDes(i);q1DDotDes(i);q2DDotDes(i)];

    for j=1:length(q1Des)
        x2 = [q1Des(j);q2Des(j);q1DotDes(j);q2DotDes(j);q1DDotDes(j);q2DDotDes(j)];
        KSS(i,j) = p1*exp(-(x1-x2)'*(x1-x2)/(2*p2^2));
    end
end

for i=1:length(q)                                                   
    x1 = [q(i,1);q(i,2);qDot(i,1);qDot(i,2);qDDot(i,1);qDDot(i,2)];

    for j=1:length(q1Des)                                                  
        x2 = [q1Des(j);q2Des(j);q1DotDes(j);q2DotDes(j);q1DDotDes(j);q2DDotDes(j)];
        KS2(i,j) = p1*exp(-(x1-x2)'*(x1-x2)/(2*p2^2));
    end
end

standart_yhat1 = std(yhat(:,1));
standart_yhat2 = std(yhat(:,2));
for i=1:length(yhat)
    up_95_1(i)  =yhat(i,1)+1.66*standart_yhat1;
    down_95_1(i)=yhat(i,1)-1.66*standart_yhat1;
    up_95_2(i)  =yhat(i,2)+1.66*standart_yhat2;
    down_95_2(i)=yhat(i,2)-1.66*standart_yhat2;
end

%%Covariance of Y Predicted
covYhat = KSS-KS*inv(K+eye(size(K)))*KS2;
covYhat = covYhat(:,1).^0.5;

for i=1:(length(q1Des)-1)             
    [qN(i+1,:),qDotN(i+1,:),qDDotN(i+1,:)] = doublependulumstep(qN(i,:),qDotN(i,:),yhat(i,:),ts);
end

%% b1 Plot of Predicted Torques

figure; hold on
plot(time,yhat(:,1),'b','LineWidth',2)
plot(time,yhat(:,2),'y','LineWidth',2)
plot(time,up_95_1,'r','LineWidth',1)
plot(time,down_95_1,'g','LineWidth',1)
plot(time,up_95_2,'r','LineWidth',1)
plot(time,down_95_2,'g','LineWidth',1)
title('Part B1. Plot of Predicted Torques and 95% Confidence')
ylabel('N.m','Fontsize',14,'interpreter','latex');
xlabel('Time','Fontsize',14,'interpreter','latex');
legend('Torque-Predicted q1','Torque-Predicted q2','95% Limit-Up','95% Limit-Down')

%% b2 Plot of Actual and Desired  Joint Angles
figure; hold on
plot(time,qN(:,1),'b','LineWidth',2)
plot(time,qN(:,2),'g','LineWidth',2)
plot(time,q1Des,'y','LineWidth',2)
plot(time,q2Des,'r','LineWidth',2)
title('Part B2. Plot of Actual and Desired  Joint Angles')
ylabel('Rad','Fontsize',14,'interpreter','latex');
xlabel('Time','Fontsize',14,'interpreter','latex');
legend('Actual q1','Actual q2','Desired q1','Desired q2')
hold off

%% Animation - Computed Torque Control
clear global
figure;
title('Computed Torque','Fontsize',14,'interpreter','latex');

l1 = 0.8;
l2 = 0.6;

q1 = qN(:,1);
q2 = qN(:,2);

x1 = l1*cos(q1);
y1 = l1*sin(q1);
x2 = l1*cos(q1) + l2*cos(q1+q2);
y2 = l1*sin(q1) + l2*sin(q1+q2);
    
pos = [x1'; y1'; x2'; y2'];

for i=1:length(time)
    px1 = pos(1,i);
    py1 = pos(2,i);
    px2 = pos(3,i);
    py2 = pos(4,i);
    
    draw(px1, py1, px2, py2);
end

%% Computed Torque Control with Feedback

Kp = 300;
Kd = 350;
for i=1:(length(q1Des)-1)
    yhatN(i,1) = yhat(i,1)+Kp*(q1Des(i)-qN(i,1))+Kd*(q1DotDes(i)-qDotN(i,1));
    yhatN(i,2) = yhat(i,2)+Kp*(q2Des(i)-qN(i,2))+Kd*(q2DotDes(i)-qDotN(i,2));
    [qN(i+1,:),qDotN(i+1,:),qDDotN(i+1,:)] = doublependulumstep(qN(i,:),qDotN(i,:),yhatN(i,:),ts);
end

%% c1 Plot of Torque Predictions                             
yhatN(length(q1Des),:)=yhat(length(q1Des),:);
figure; hold on;
plot(time,yhatN(:,1),'b','LineWidth',2)
plot(time,yhatN(:,2),'g','LineWidth',2)
title('Part C1. Plot of Original and Predicted Torques')
ylabel('N.m','Fontsize',14,'interpreter','latex');
xlabel('Time','Fontsize',14,'interpreter','latex');
legend('Torque-Predicted q1','Torque-Predicted q2')

%% c2 Plot of Actual and Desired  Joint Angles
figure; hold on;
plot(time,qN(:,1),'b','LineWidth',2)
plot(time,qN(:,2),'g','LineWidth',2)
plot(time,q1Des,'y','LineWidth',2)
plot(time,q2Des,'r','LineWidth',2)
title('Part C2. Plot of Actual and Desired  Joint Angles')
ylabel('Rad','Fontsize',14,'interpreter','latex');
xlabel('Time','Fontsize',14,'interpreter','latex');
legend('Actual q1','Actual q2','Desired q1','Desired q2')
hold off

%% c) Computed Torque Control with Feedback %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear global
figure;
title('Computed Torque','Fontsize',14,'interpreter','latex');

l1 = 0.8;
l2 = 0.6;

q1 = qN(:,1);
q2 = qN(:,2);

x1 = l1*cos(q1);
y1 = l1*sin(q1);
x2 = l1*cos(q1) + l2*cos(q1+q2);
y2 = l1*sin(q1) + l2*sin(q1+q2);
    
pos = [x1'; y1'; x2'; y2'];

for i=1:length(time)
    px1 = pos(1,i);
    py1 = pos(2,i);
    px2 = pos(3,i);
    py2 = pos(4,i);
    
    draw(px1, py1, px2, py2);
end 

display(['Done Part 2']);

% Functions --------------------------------------------------------------
function draw(px1, py1, px2, py2)

    extents = [-2 2 -2 2];

    global link1Handle link2Handle joint1Handle joint2Handle circleHandle

    cx = 0;
    cy = 1;
    r = 0.25;
    
    if isempty(circleHandle)
        circleHandle = rectangle(...
            'Position',[cx-r, cy-r, 2*r, 2*r],...
            'Curvature',[1,1],...   % <-- Draws a circle...
            'LineWidth',2,...
            'FaceColor',[0.3, 0.2, 0.7],...
            'EdgeColor',0.5*[0.3, 0.2, 0.7]);  %Darker version of color
    else
        set(circleHandle,...
            'Position', [cx-r, cy-r, 2*r, 2*r]);
    end
    
        
    if isempty(link1Handle)
        link1Handle = line([0 px1], [0 py1], ...
            'LineWidth',4, ...
            'Color',[0.2, 0.7, 0.2]);
    else
        set(link1Handle, ...
            'xData',[0 px1], ...
            'yData',[0 py1]);
    end
        
    if isempty(link2Handle)
        link2Handle = line([px1 px2], [py1 py2], ...
            'LineWidth',4, ...
            'Color',[0.2, 0.7, 0.2]);
    else
        set(link2Handle, ...
            'xData',[px1 px2], ...
            'yData',[py1,py2]);
    end
    
    r = 0.05;
    if isempty(joint1Handle)
        joint1Handle = rectangle(...
            'Position',[px1-r, py1-r, 2*r, 2*r],...
            'Curvature',[1,1],...   % <-- Draws a circle...
            'LineWidth',2,...
            'FaceColor',[0.3, 0.2, 0.7],...
            'EdgeColor',0.5*[0.3, 0.2, 0.7]);  %Darker version of color
    else
        set(joint1Handle,...
            'Position',[px1-r, py1-r, 2*r, 2*r]);
    end
    
    if isempty(joint2Handle)
        joint2Handle = rectangle(...
            'Position',[px2-r, py2-r, 2*r, 2*r],...
            'Curvature',[1,1],...   % <-- Draws a circle...
            'LineWidth',2,...
            'FaceColor',[0.3, 0.2, 0.7],...
            'EdgeColor',0.5*[0.3, 0.2, 0.7]);  %Darker version of color
    else
        set(joint2Handle,...
            'Position',[px2-r, py2-r, 2*r, 2*r]);
    end

    axis equal; axis(extents); axis on;      %  <-- Order is important here
    drawnow;
end

%% Covariance function
function k = cov_func(x, xp, p1, p2)
    k = p1*exp( -(x-xp)'*(x-xp) / (2*p2^2) );
end

% Kernel function
function k = kernal(x,q,M,h)
    d = ((x-q)*M'*M*(x-q)')^0.5;
    k = exp(-((d^2)/h));
end

