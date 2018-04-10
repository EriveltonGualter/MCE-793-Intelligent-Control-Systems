%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dmpPhasePortrait.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: This function plots a phase portrait for transformation
% system of a DMP.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%   parameters  = [tau alphaZ betaZ g]
%   yrange      = [yMin yMax]
%   zrange      = [zMin zMax]
%   f           = value of nonlinear function
%   figNum      = number of figure in which to plot the phase portrait
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Eric Schearer
% Date created: 4/18/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dmpPhasePortrait(yprev,ydotprev,parameters,yrange,zrange,f,figNum)
close all
% parameters
tau     = parameters(1);
alphaZ  = parameters(2);
betaZ   = parameters(3);
g       = parameters(4);

% set up vectors of the position and velocity vectors
y = linspace(yrange(1),yrange(2),21);
z = linspace(zrange(1),zrange(2),21);

% this makes a grid of all the combos of positions and velocities
[y1,y2] = meshgrid(y,z);

% initialize the vectors for the phase plot
u = zeros(size(y1));
v = zeros(size(y1));

% compute the derivatives
for i = 1:numel(y1)
    [u(i),v(i)] = dmpderivative(y1(i),y2(i),g,f,alphaZ,betaZ,tau);
end

% plot the phase portrait
hold off
figure(figNum)
hold on
quiver(y1,y2,u,v,'r');
plot(yprev,ydotprev)
plot(yprev(1),ydotprev(1),'go','MarkerFaceColor','g')
plot(yprev(end),ydotprev(end),'ko','MarkerFaceColor','k')
xlabel('y')
ylabel('z')
axis tight equal;

