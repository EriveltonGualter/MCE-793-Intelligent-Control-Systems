% Erivelton Gualter dos Santos

clc, clear all, close all

% Bayesian Regression %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('pendulumdata.mat'); % Load  g, q, qDDOT, qDot, T

g= 9.81;   % gravitational acceleration

wo = [1; 1];        % Prior Estimate (prior mean)
S0 = 10*eye(2,2);   % Prior Covariance

binv = inv(0.01);   % Sensor Variance

% phi = g*cos(q)/2 qDDOT/3

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

figure; plot(var, 'LineWidth', 2); 
title('Part (a1)', 'interpreter', 'latex', 'FontSize', 14);
xlim([0 length(T)]);

%%
Phi = [g*cos(q)/2 qDDOT/3];
Wml = pinv(Phi)*T;

l = Wml(2)/Wml(1);
m = Wml(1)/l;

display(['Part (a3): m = ',num2str(m)])
display(['Part (a3): l = ',num2str(l)])

%% Computed Torque Control %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
load('desiredtrajecotry.mat');

m = 0.6; % kg
l = 4.0; % m

g= 9.81;   % gravitational acceleration

tstep = 1e-3; % Time Step [s] 

%plot(time, q(1:6001), time, qDes, 'LineWidth', 2); 

% qDes qDotDes qDDotDes
for k=1 : length(time)
    
    TDes(k) = m*l^2*qDDotDes(k)/3 + m*g*l*cos(qDes(k))/2;
    
end

plot(time, qDes, '--', 'LineWidth', 2);
xlabel('Time [s]', 'interpreter', 'latex', 'FontSize', 14);
ylabel('Angle', 'interpreter', 'latex', 'FontSize', 14);
title('Part (b1)', 'interpreter', 'latex', 'FontSize', 14);

legend('Desired')

% display(['Part (b2): rms error = ',num2str(rmserror)])



%% Simulation
clear all
close all

load('pendulumdata.mat'); % Load  g, q, qDDOT, qDot, T
load('desiredtrajecotry.mat');

g= 9.81;   % gravitational acceleration

m = 0.6; % kg
l = 4.0; % m

display(['Part (a3): m = ',num2str(m)])
display(['Part (a3): l = ',num2str(l)])

q1 = qDes;
t = time ; 

x1 = l*cos(q1);
y1 = l*sin(q1);
    
pos = [x1'; y1'];

time = 0;
tic;
while time < t(end)

    % Compute the position of the system at the current real world time
    posDraw = interp1(t',pos',time')';
    px = posDraw(1);
    py = posDraw(2);
    
    if not(isnan(posDraw))
       draw(px, py);
    end
    
    % Update current time
    time = toc;
end 

function draw(px, py)

    extents = [-5 5 -5 5];

    global linkHandle jointHandle 

    if isempty(linkHandle)
        linkHandle = line([0 px], [0 py], ...
            'LineWidth',4, ...
            'Color',[0.2, 0.7, 0.2]);
    else
        set(linkHandle, ...
            'xData',[0 px], ...
            'yData',[0 py]);
    end
           
    r = 0.25;
    if isempty(jointHandle)
        jointHandle = rectangle(...
            'Position',[px-r, py-r, 2*r, 2*r],...
            'Curvature',[1,1],...   % <-- Draws a circle...
            'LineWidth',2,...
            'FaceColor',[0.3, 0.2, 0.7],...
            'EdgeColor',0.5*[0.3, 0.2, 0.7]);  %Darker version of color
    else
        set(jointHandle,...
            'Position',[px-r, py-r, 2*r, 2*r]);
    end
    
    axis equal; axis(extents); axis on;      %  <-- Order is important here
    drawnow;
end


