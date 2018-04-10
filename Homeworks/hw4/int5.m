close all
clear
clc

In = 2;
Hid= 100;
Out= 1;

x = 0:pi/20:4*pi;
x = x';
y = sin(x)+0.05*randn(size(x));
plot(x,y)

x =[x,ones(length(x))];

%%
W1 = rand(Hid,In);
W2 = rand(Hid,Out);
grad = zeros(2*Hid,In)+2;
yhat = zeros(size(x,1),1);
err = Inf;

while err>0.1
    
    for i=1:Hid
        a(:,i) = x(:,1)*W1(i,1) + x(:,2)*W1(i,2);
        z(:,i) = 1./(1+exp(-a(:,i)));
        yhat(:,1) = yhat(:,1)+z(:,i)*W2(i,1);
    end
    
    Delta = yhat-y;
    
    for i=1:Hid
        Delta1(:,i) = exp(-a(:,i))./((1+exp(-a(:,i))).^2).*Delta(:,1).*W2(i,1);
    end
    
    for i=1:Hid
        dEndw1(:,i) = Delta(:,1).*z(:,i);
        dEndw11(:,i)= Delta1(:,i).*x(:,1);
        dEndw12(:,i)= Delta1(:,i).*x(:,2);
    end
        
    grad = [sum(dEndw11) ; sum(dEndw12) ; sum(dEndw1)]';
        
    eta =0.001;
    W1 = W1-eta*grad(1:Hid,:);
    W2 = W2-eta*grad(Hid+1:end,:);
    err = (sum(Delta.^2)/size(x,1)).^0.5
end

hold on

plot(x(:,1),yhat)