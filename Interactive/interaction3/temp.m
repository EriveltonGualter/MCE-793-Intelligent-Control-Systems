clear all
clc
close all

load('pendulumdata.mat')

zz=1:100:length(q);
qDDot=qDDOT(zz');
q=q(zz');
qDot=qDot(zz');
T=T(zz');

maxq=max(q);
minq=min(q);
maxqDot=max(qDot);
minqDot=min(qDot);
maxqDDot=max(qDDot);
minqDDot=min(qDDot);

m1=1/(maxq-minq);
m2=1/(maxqDot-minqDot);
m3=1/(maxqDDot-minqDDot);

M=diag([1 m1 m2 m3]);

for ii=1:length(q)
	qq=[1;q(ii);qDot(ii);qDDot(ii)];
	h=0.1;

	for i=1:length(q)
		s=(kernal([1;q(i);qDot(i);qDDot(i)],qq,M,h))^0.5;
		Z(i,:)=(s*[1;q(i);qDot(i);qDDot(i)])';
		V(i,:)=(s*T(i))';
	end

	what=((Z')*Z)\(Z'*V);
	yhat(ii)=what'*qq;
end

figure
plot(zz,T,zz,yhat,'*')
legend('T','yhat')

figure,scatter(T,yhat)
