
clc
clear all
close all

x = linspace(0,4*pi,100);
x = x';
y = sin(x)+0.05*randn(size(x));

x = [x, ones(length(x),1)];

m = 5;

[l,b] = size(x);
[n,a] = size(y);

W1 = rand(m,b); 
W2 = rand(m,a);

New_Error  =  Inf;
buffer = 0;
count = 1;

while New_Error>0.1
    
    a = x*W1';
    z = 1./(1+exp(-a));
    yhat = z*W2;

    delta21 = yhat-y;
    for i=1:m
        delta1(:,i)=exp(-a(:,i))./((1+exp(-a(:,i))).^2).*delta21(:,1).*W2(i,1);
    end
%     delta1 = exp(-a)./((1+exp(-a)).^2).*delta21.*W2;
    
    dEndw2 = delta21.*z;
    dEndw11 = delta1.*x(:,1);
    dEndw12 = delta1.*x(:,2);

    grad1 = [sum(dEndw11); sum(dEndw12)]';
    grad2 = sum(dEndw2)';

    eta = 0.005;
    W1 = W1-eta*grad1; W1a(count,:,:) = W1;
    W2 = W2-eta*grad2; W2a(count,:,:) = W2;

    count = count + 1;
%     New_Error = sum(sqrt((sum(y-yhat).^2)/length(y)))/2
    New_Error = rms(yhat-y)
    id = mod(count, 20);
    if id == 0
        if std(buffer) < 0.00001
            break
        end
    else
        buffer(id) = New_Error;
    end
end
%%
close all
% W1plot(:,:) = W1a(:,:,1);
% subplot(211); plot(W1plot,'LineWidth',2); xlim([0 count]);
% subplot(212); plot(W2a,'LineWidth',2);    xlim([0 count]);

figure
plot(x(:,1),y,x(:,1),yhat,'LineWidth',2); xlim([x(1,1) x(end,1)]);
