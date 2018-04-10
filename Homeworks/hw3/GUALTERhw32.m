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

% id = 1:100:length(q);
% x = [q(id') qDot(id') qDDot(id')];
% y = tau(id');

x = [q(:,1) q(:,2) qDot(:,1) qDot(:,2) qDDot(:,1) qDDot(:,2)];
y1 = tau(:,1);
y2 = tau(:,2);

id = 1:10:length(q);
x = x(id,:);
y1 = y1(id);
y2 = y2(id);

% Given Matrix 6x6
M = diag([sqrt(200), sqrt(200), sqrt(10), sqrt(10), 1, 1]);

% Parameter to tune: (Bandwidth hyperparameter) 
%harray = [100 1000 10000];
h = 100;

for i=1:length(q)
    qq = x(i,:);

    for j=1:length(q)-1
        xq = x;
        xq(i,:) = [];
        
        d(i,j)=((xq-qq)'*M'*M*(xq-qq))^0.5;
        
        s=(exp(-((d(i,j)^2)/h)))^0.5;
        Z(j,:)=(s*xq(j,:))';
        V1(j,:)=(s*tau10(j,1))';
        V2(j,:)=(s*tau10(j,2))';
    end

    W_hat1 = pinv(((Z')*Z))*((Z')*V1);
    W_hat2 = pinv(((Z')*Z))*((Z')*V2);

    Y_hat(i,1) = W_hat1'*qq;
    Y_hat(i,2) = W_hat2'*qq;

end
error= sqrt(median(median((tau10-Y_hat).^2)))

LOOCV=error;
display(['LOOCV RMS error for h=100 is ',num2str(LOOCV)])

plot(Y_hat(:,1),tau10(:,1),'ro')

%%
% display(['LOOCV RMS error for h = 100 is ',num2str(LOOCV(1))]);
% display(['LOOCV RMS error for h = 1000 is ',num2str(LOOCV(2))]);
% display(['LOOCV RMS error for h = 10000 is ',num2str(LOOCV(3))]);
% LOOCV
% error_1000=sqrt(median(median(([yhat1(:,hi) yhat2(:,hi)] - [y1 y2]).^2)));
% 
% tau10 = tau(id,:);
% error= sqrt(median(median((tau10-Y_hat).^2)))


[M,I] = min(LOOCV);

figure; hold on;
plot(yhat1(:,hi), y1, 'o'); plot(yhat2(:,hi), y2, '*');
title('Part a2','Fontsize',14,'interpreter','latex');
ylabel('Actual Torques','Fontsize',14,'interpreter','latex');
xlabel('Predicted Torques','Fontsize',14,'interpreter','latex');

% b) Computed Torque Control %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% c) Computed Torque Control with Feedback %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% a) Gaussian process regression %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% b) Computed Torque Control %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% c) Computed Torque Control with Feedback %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Covariance function
function k = cov_func(x, xp, p1, p2)
    k = p1*exp( -(x-xp)'*(x-xp) / (2*p2^2) );
end

% Kernel function
function k = kernal(x,q,M,h)
    d = ((x-q)*M'*M*(x-q)')^0.5;
    k = exp(-((d^2)/h));
end

