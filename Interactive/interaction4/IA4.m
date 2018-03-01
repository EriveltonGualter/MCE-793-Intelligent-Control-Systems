clear all
clc
close all

%(a)
load pendulumdata

%(b)
sam=200;
t=1:sam:length(qDDOT);
qDDot=qDDOT(1:sam:length(qDDOT));
q=q(1:sam:length(q));
qDot=qDot(1:sam:length(qDot));
T=T(1:sam:length(T));

p1=1;p2=30;
m=zeros(3,length(q));
for i=1:length(q)
    x1=[q(i);qDot(i);qDDot(i)];
   
    for j=1:length(q)
            x2=[q(j);qDot(j);qDDot(j)];
        K(i,j)=p1*exp(-(x1-x2)'*(x1-x2)/(2*p2^2));
        
    end
end

f=mvnrnd(m,K);
scatter(q,f(1,:))
hold on

%(c)
p1=1;p2=1;
m=zeros(3,length(q));
for i=1:length(q)
    x1=[q(i);qDot(i);qDDot(i)];
    for j=1:length(q)
            x2=[q(j);qDot(j);qDDot(j)];
        Ks(1,j)=p1*exp(-(x1-x2)'*(x1-x2)/(2*p2^2));
        
    end
    yhat(i)=Ks*inv(K+1*eye(size(K)))*T;
end
scatter(T,yhat)
plot(yhat,'r')
hold on
plot(T)

% (d)
%p1=1;p2=0.001;
m=zeros(3,length(q));
for i=1:length(q)
    x1=[q(i);qDot(i);qDDot(i)];
    for j=1:length(q)-2
        qt=q;
        qt(i)=[];
                qtDot=qDot;
        qtDot(i)=[];
                qtDDot=qDDot;
        qtDDot(i)=[];
            x2=[qt(j);qtDot(j);qtDDot(j)];
            
            Kt=K;
            Kt(:,i)=[];
            Kt(i,:)=[];
            Tt=T;
            Tt(i)=[];
        Ks(1,j)=p1*exp(-(x1-x2)'*(x1-x2)/(2*p2^2));
        
    end
    yhat(i)=Ks(1:250)*inv(Kt+0*eye(size(Kt)))*Tt;
    i
end
scatter(T,yhat)

Err=sqrt((sum(T'-yhat).^2)/length(T))
Err2=sqrt(median(median((T'-yhat).^2)))
