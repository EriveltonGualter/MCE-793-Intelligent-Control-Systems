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

