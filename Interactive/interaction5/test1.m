
close all
clear all
clc

x = 0:pi*0.001:4*pi;
x = x';
y = sin(x)+0.05*randn(size(x)); % trainning data set
plot(x,y)

W1=[1;1;1];
W2=[1;1;1];
grad=1;
count = 1;
while norm(grad)>0.001

    a1=x*W1(1);
    a2=x*W1(2);
    a3=x*W1(3);

    z1=1./(1+exp(-a1));
    z2=1./(1+exp(-a2));
    z3=1./(1+exp(-a3));

    yhat=z1*W2(1)+z2*W2(2)+z3*W2(3);

    delta21=yhat-y;
    delta11=exp(-a1)./((1+exp(-a1)).^2).*W2(1).*delta21;
    delta12=exp(-a2)./((1+exp(-a2)).^2).*W2(2).*delta21;
    delta13=exp(-a3)./((1+exp(-a3)).^2).*W2(3).*delta21;

    dEndw21=delta21.*z1;
    dEndw22=delta21.*z2;
    dEndw23=delta21.*z3;
    
    dEndw11=delta11.*x;
    dEndw12=delta12.*x;
    dEndw13=delta13.*x;

    grad=[sum(dEndw11);sum(dEndw12);sum(dEndw13);sum(dEndw21);sum(dEndw22);sum(dEndw23)];

    eta=0.001;
    W1=W1-eta*grad(1:3);
    W2=W2-eta*grad(4:6);

    norm(grad);
    count= count+1;
end

hold on; plot(x,yhat); ylim([-pi/2 pi/2]);