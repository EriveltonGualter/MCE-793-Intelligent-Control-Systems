function mySinglePendulumAnimation( t, q, l )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%q = X(1,2:end);
x1 = 0; x2 = l*sin(q);
y1 = 0; y2 = l*cos(q);
    
pos = [x2; y2];

extents = [-1 1 -1 1];

i = 1;
time = 0;
tic;
while time < t(end)

    % Compute the position of the system at the current real world time
    posDraw = interp1(t',pos',time')';
    
    if not(isnan(posDraw))
        pole = plot([x1,posDraw(1,1)], [y1,posDraw(2,1)], 'LineWidth',4, 'Color',[0.2, 0.7, 0.2]);
        axis equal; axis(extents); axis on;      %  <-- Order is important here
        drawnow;
        pause(1/40);
    end
    
    % Update current time
    time = toc;
end 
