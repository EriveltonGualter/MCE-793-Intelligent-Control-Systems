
clc
clear

x = 0:pi/20:4*pi;
x = x';
y = sin(x);%+0.05*randn(size(x));
plot(x,y)

x=[x,ones(length(x),1)];

nn = 10;
W1 = rand(Hid,In);
W2 = rand(Hid,Out);
grad = 1;
yhat = zeros(size(x,1),1);
err = Inf;

while err>0.1

    for i=1:nn
        a(:,i)= x(:,1)*W1(i,1) + x(:,2)*W1(i,2) ;
        z(:,i)=1./(1+exp(-a(:,i)));
        yhat=yhat+z(:,i)*W2(i);
    end

    delta21 = yhat-y;
    
    for i=1:nn
        delta1(:,i) = exp(-a(:,i))./((1+exp(-a(:,i))).^2).*W2(i).*delta21;
        dEndw2(:,i) = delta21.*z(:,i);
        dEndw11(:,i)= delta1(:,i).*x(:,1); 
        dEndw12(:,i)= delta1(:,i).*x(:,2); 
    end
    
    grad1 = [sum(dEndw11);sum(dEndw12)]';
    grad2 = sum(dEndw2)';

    eta=0.001;
    W1 = W1-eta*grad1;
    W2 = W2-eta*grad2;
    err= (sum(delta21.^2)/size(x,1)).^0.5

end

hold on; plot(x(:,1),yhat)