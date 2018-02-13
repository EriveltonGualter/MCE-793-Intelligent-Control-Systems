%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GUALTERhw1.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: Homework 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Erivelton Gualter
% Date created: 1/24/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

% Parameters
m1 = 1;     % kg
m2 = 1;     % kg
l1 = 1;     % m
l2 = 1;     % m
r1 = 0.5;   % m
r2 = 0.5;   % m
I1 = m1*l1^2/12; % m^3
I2 = m2*l2^2/12; % m^3
g = 9.81;

% Simulation parameters
tf = 5;     % Final time [s]
dt = 1e-3;  % Integration step size [s]
t = 0:dt:tf; % array time

% Initial Conditional
X10 = pi/4;
X20 = 0;
X30 = 0;
X40 = 0;

X = [X10; X20; X30; X40];

for i=1 : length(t)-1
    
    x1 = X(1,i);
    x2 = X(2,i);
    x3 = X(3,i);
    x4 = X(4,i);
    
    Q1=x1; Q2=x3; Q1dot=x2; Q2dot = x4;
    accel = compute_accel(I1,I2,x1,x3,x2,x4,0,0,g,l1,m1,m2,r1,r2);
    
    xd1 = x2;
    xd2 = accel(1);  
    xd3 = x4;
    xd4 = accel(2);
    
    XDOT = [xd1; xd2; xd3; xd4];
    
    X(:,i+1) = X(:,i) + XDOT* dt;
    
    T(i) = compute_kinetic(I1,I2,Q2,Q1dot,Q2dot,l1,m1,m2,r1,r2);
    U(i) = compute_potential(Q1,Q2,g,l1,m1,m2,r1,r2);
end

%%
close all
q1 = X(1,:);
q2 = X(3,:);

x1 = l1*cos(q1);
y1 = l1*sin(q1);
x2 = l1*cos(q1) + l2*cos(q1+q2);
y2 = l1*sin(q1) + l2*sin(q1+q2);
    
pos = [x1; y1; x2; y2];

i = 1;
time = 0;
tic;
while time < t(end)

    % Compute the position of the system at the current real world time
    posDraw = interp1(t',pos',time')';
    px1 = posDraw(1);
    py1 = posDraw(2);
    px2 = posDraw(3);
    py2 = posDraw(4);    
    
    if not(isnan(posDraw))
       draw(px1, py1, px2, py2);
    end
    
    % Update current time
    time = toc;
end 

%%
close all
figure; plot(t,X(1,:), t,X(3,:), 'LineWidth', 2)
        title('Joint Positions', 'Interpreter','Latex', 'FontSize',14);
        xlabel('Time [s]', 'Interpreter','Latex', 'FontSize',14);
        ylabel('[rad]', 'Interpreter','Latex', 'FontSize',14);
        legend('Q1','Q2');

figure; plot(t,X(2,:), t,X(4,:), 'LineWidth', 2)
        title('Joint Velocities', 'Interpreter','Latex', 'FontSize',14);
        xlabel('Time [s]', 'Interpreter','Latex', 'FontSize',14);
        ylabel('[rad/s]', 'Interpreter','Latex', 'FontSize',14);
        legend('dQ1','dQ2');

figure; plot(t(2:end), T, t(2:end), U, t(2:end), T+U, 'LineWidth', 2)
        title('Energy', 'Interpreter','Latex', 'FontSize',14);
        xlabel('Time [s]', 'Interpreter','Latex', 'FontSize',14);
        ylabel('Energy [J]', 'Interpreter','Latex', 'FontSize',14);
        legend('Kinetic Energy','Potential Energy','Total Energy');

        %legend('','','');
        %legend('','','');

%%
function draw(px1, py1, px2, py2)

    extents = [-2 2 -2 2];

    global link1Handle link2Handle joint1Handle joint2Handle

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