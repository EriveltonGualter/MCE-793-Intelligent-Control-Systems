clear
clc

%% b) Computed Torque Control %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('desiredtrajecotry.mat');
load('doublependulumdata.mat'); % Load  q, qDDOT, qDot, tau

% Resampling

id = 1:10:length(q1Des);

q1DDotDes = q1DDotDes(id);
q1Des = q1Des(id);
q1DotDes = q1DotDes(id);
q2DDotDes = q2DDotDes(id);
q2Des = q2Des(id);
q2DotDes = q2DotDes(id);
time = time(id);
xDes = xDes(id);
yDes = yDes(id);

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

%% B1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
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
title('Part $B_2$','Fontsize',14,'interpreter','latex');
ylabel('Rad','Fontsize',14,'interpreter','latex');
xlabel('Time','Fontsize',14,'interpreter','latex');
legend('actual-q1','actual-q2','desired-q1','desired-q2');
xlim([time(1) time(end)]);

%% B3. Animation
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

%% c. Animation
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

%%
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