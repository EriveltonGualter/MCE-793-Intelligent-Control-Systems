% Erivelton Gualter dos Santos

clc, clear all, close all

load('pendulumdata.mat')

g = 9.81;

% Decreasing the data with a 100 sample step
id = 1:500:length(q);
q       = q(id');
qDot    = qDot(id');
qDDOT   = qDDOT(id');
T       = T(id');

q_min = min(q); q_max = max(q);
qrange = q_max - q_min; 

qDot_min = min(qDot); qDot_max = max(qDot);
qDotrange = qDot_max - qDot_min; 

qDDOT_min = min(qDDOT); qDDOT_max = max(qDDOT);
dDDotrange = qDDOT_max - qDDOT_min;

% Display the range just to have an ideia what's going on :D
display(['    q should be between ',num2str(q_min), ' and ',num2str(q_max)]);
display([' qDot should be between ',num2str(qDot_min), ' and ',num2str(qDot_max)]);
display(['qDDot should be between ',num2str(qDDOT_min), ' and ',num2str(qDDOT_max)]);

ranges = [qrange; qDotrange; dDDotrange];
ranges = sort(ranges)

% Probably it should be 4x4 bc we have 1 in the first pos, but it is 3x3
M = diag([1; 1./ranges])

% Parameter to tune: (Bandwidth hyperparameter) 
h = 0.002;

X = [1; q; qDot; qDDOT];

error = 0;
for ii=1:length(q)
	qq=[1;q(ii);qDot(ii);qDDOT(ii)];

	for i=1:length(q)
		s = (kernal([1;q(i);qDot(i);qDDOT(i)],qq,M,h))^0.5;
		Z(i,:) = (s*[1;q(i);qDot(i);qDDOT(i)])';
		V(i,:) = (s*T(i))';
	end

	what=((Z')*Z)\(Z'*V);
	yhat(ii)=what'*qq;
    
    error = error + (T(ii)-yhat(ii))^2;
end
error


figure; plot(id,T,id,yhat,'+'); legend('T','yhat');
figure; plot(T,yhat, '+')


function k = kernal(x, q, M, h)
    d = ((x-q)'*M'*M*(x-q))^0.5;
    k = exp(-((d^2)/h));
end

% function out = get_gaussian_kernal( x, M, q, h)
%     d = sqrt((x-q)*M*M*(x-q)');
%     out = exp(-d^2/h);
% end

% function out = get_euclid_dist ( x, q)
%     out = sqrt((x-q)'*(x-q))
% end
