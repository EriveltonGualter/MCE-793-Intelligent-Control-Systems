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

% Bayesian Regression %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('pendulumdata.mat'); % Load  g, q, qDDOT, qDot, T

g= 9.81;   % gravitational acceleration

wo = [1; 1];        % Prior Estimate (prior mean)
S0 = 10*eye(2,2);   % Prior Covariance

binv = inv(0.01);   % Sensor Variance

t = [T];

SN = zeros(2,2,length(T));
MN = zeros(length(T),2);

for k=1:size(T)
    
    if k>1
        phi = vertcat(phi, [g*cos(q(k))/2 qDDOT(k)/3]); 
        Tk = vertcat(Tk, t(k));
    else
        phi = [g*cos(q(k))/2 qDDOT(k)/3];
        Tk = t(k);
    end
    
    SN(:,:,k) = inv( inv(S0) + binv*phi'*phi );
    MN(k,:) = SN(:,:,k) * (inv(S0)*wo + binv*phi'*Tk);  
    
    var(k,1) = SN(1,1,k).^0.5;
    var(k,2) = SN(2,2,k).^0.5;
       
    % inv(SN): Posterior Confidence
    % MN : Posterior Mean% Posterior Mean
end

figure; plot(MN, 'LineWidth', 2); 
title('Part (a1)', 'interpreter', 'latex', 'FontSize', 14);
xlim([0 length(T)]);

figure; hold on; plot(var(:,1), 'LineWidth', 4); plot(var(:,2), 'LineWidth', 2); 
title('Part (a2)', 'interpreter', 'latex', 'FontSize', 14);
legend('var1', 'var2'); 
xlim([0 length(T)]);

%%
Phi = [g*cos(q)/2 qDDOT/3];
Wml = pinv(Phi)*T;

l = Wml(2)/Wml(1);
m = Wml(1)/l;

display(['Part (a3): m = ',num2str(m)])
display(['Part (a3): l = ',num2str(l)])

%% Computed Torque Control %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('desiredtrajecotry.mat');

m = 0.6; % kg
l = 4.0; % m

g = 9.81;   % gravitational acceleration

tstep = 1e-3; % Time Step [s] 

%plot(time, q(1:6001), time, qDes, 'LineWidth', 2); 
   
X10 = qDes(1);
X20 = qDotDes(1);

rmserror = 0;
X = [X10; X20];
for k=1 : length(time)-1
        
    rmserror = rmserror + (X(1,k)-qDes(k))^2;
    
    x1 = X(1,k);
    x2 = X(2,k);
    
    TDes(k) = m*l^2*qDDotDes(k)/3 + m*g*l*cos(x1)/2;
        
    xd1 = x2;
    xd2 = (TDes(k) -  m*g*l*cos(x1)/2) * 3/(m*l^2);
    
    XDOT = [xd1; xd2];
    
    X(:,k+1) = X(:,k) + XDOT*tstep;
end

figure; hold on;
plot(time, X(1,:)', 'LineWidth', 2);
plot(time, qDes, '--', 'LineWidth', 2);

xlabel('Time [s]', 'interpreter', 'latex', 'FontSize', 14);
ylabel('Angle', 'interpreter', 'latex', 'FontSize', 14);
title('Part (b1)', 'interpreter', 'latex', 'FontSize', 14);
legend('Computed', 'Desired');

rmserror = sqrt(rmserror);
display(['Part (b2): rms error = ',num2str(rmserror)])

%% Simulation 1
figure;

tsim = time ; 

xact = l*cos(X(1,:)');
yact = l*sin(X(1,:)');
xdes = l*cos(qDes);
ydes = l*sin(qDes);
    
pos = [xact'; yact'; xdes'; ydes'];

timesim = 0;
tic;
while timesim < tsim(end)

    % Compute the position of the system at the current real world time
    posDraw = interp1(tsim',pos',timesim')';
    pxdes = posDraw(1);
    pydes = posDraw(2);
    pxact = posDraw(3);
    pyact = posDraw(4);
    
    if not(isnan(posDraw))
       draw(pxact, pyact, pxdes, pydes, 1);
    end
    
    % Update current time
    timesim = toc;
end

%% Computed Torque Control with Feedback %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X10 = qDes(1);
X20 = qDotDes(1);

Kp = 120; Kd = 15;

X = [X10; X20];
rmserror = 0;
for k=1 : length(time)-1
    
    rmserror = rmserror + (X(1,k)-qDes(k))^2;
    
    x1 = X(1,k);
    x2 = X(2,k);
    
    Tcon(k) = TDes(k) + Kp*(qDes(k)-x1) + Kd*(qDotDes(k)-x2);
    
    xd1 = x2;
    xd2 = (Tcon(k) -  m*g*l*cos(x1)/2) * 3/(m*l^2);
    
    XDOT = [xd1; xd2];
    
    X(:,k+1) = X(:,k) + XDOT*tstep;
end

figure; hold on;
plot(time, X(1,:)', 'LineWidth', 2);
plot(time, qDes, '--', 'LineWidth', 2);

xlabel('Time [s]', 'interpreter', 'latex', 'FontSize', 14);
ylabel('Angle', 'interpreter', 'latex', 'FontSize', 14);
title('Part (c1)', 'interpreter', 'latex', 'FontSize', 14);
legend('Computed', 'Desired');

rmserror = sqrt(rmserror);
display(['Part (c2): rms error = ',num2str(rmserror)])

%% Simulation 2
clear global
figure

tsim = time ; 

xact = l*cos(X(1,:)');
yact = l*sin(X(1,:)');
xdes = l*cos(qDes);
ydes = l*sin(qDes);
    
pos = [xact'; yact'; xdes'; ydes'];

timesim = 0;
tic;
while timesim < tsim(end)

    % Compute the position of the system at the current real world time
    posDraw = interp1(tsim',pos',timesim')';
    pxdes = posDraw(1);
    pydes = posDraw(2);
    pxact = posDraw(3);
    pyact = posDraw(4);
    
    if not(isnan(posDraw))
       draw(pxact, pyact, pxdes, pydes, 2);
    end
    
    % Update current time
    timesim = toc;
end 

%% Draw function
function draw(pxact, pyact, pxdes, pydes, flag)

    extents = [-5 5 -5 5];
	
    global link1Handle link2Handle jointHandle desHandle

    if isempty(link1Handle)
        link1Handle = line([0 pxact], [0 pyact], ...
            'LineWidth',4, ...
            'Color',[0.3, 0.2, 0.7]);
    else
        set(link1Handle, ...
            'xData',[0 pxact], ...
            'yData',[0 pyact]);
    end
           
    r = 0.25;
    if isempty(jointHandle)
        jointHandle = rectangle(...
            'Position',[pxact-r, pyact-r, 2*r, 2*r],...
            'Curvature',[1,1],...   % <-- Draws a circle...
            'LineWidth',2,...
            'FaceColor',[0.3, 0.2, 0.7],...
            'EdgeColor',0.5*[0.3, 0.2, 0.7]);  %Darker version of color
    else
        set(jointHandle,...
            'Position',[pxact-r, pyact-r, 2*r, 2*r]);
    end
    
    if isempty(link2Handle)
        link2Handle = line([0 pxdes], [0 pydes], ...
            'LineStyle','--', ...
            'LineWidth',3, ...
            'Color',[1, 0, 0]);
    else
        set(link2Handle, ...
            'xData',[0 pxdes], ...
            'yData',[0 pydes]);
    end
    
    r = 0.25;
    if isempty(desHandle)
        desHandle = rectangle(...
            'Position',[pxdes-r, pydes-r, 2*r, 2*r],...
            'Curvature',[1,1],...   % <-- Draws a circle...
            'LineWidth',2,...
            'FaceColor',[1, 0, 0],...
            'EdgeColor',0.5*[1, 0, 0]);  %Darker version of color
    else
        set(desHandle,...
            'Position',[pxdes-r, pydes-r, 2*r, 2*r]);
    end
    
    legend('Computed','Desired');
    
    if flag == 1
        title('Computed Torque Control', 'interpreter', 'latex', 'FontSize', 14);
    else
        title('Computed Torque Control with Feedback', 'interpreter', 'latex', 'FontSize', 14);
    end
    
    axis equal; axis(extents); axis on;      %  <-- Order is important here
    drawnow;
end


