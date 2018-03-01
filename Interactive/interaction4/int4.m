% Erivelton Gualter dos Santos

clear all, close all

load('pendulumdata.mat')

p1 = 1;
p2 = 20.65;

id = 1:100:length(q);
x = [q(id') qDot(id') qDDOT(id')];
y = T(id');

K = zeros(length(x),length(x));
for i=1:length(x)
    for j=1:length(x)
        K(i,j) = cov_func(x(i,:)', x(j,:)', p1, p2);
    end
end

%%
sigma = var(x');
hold on
Ks1 = zeros(1,length(x));
Ks2 = zeros(length(x),1);
for i=1:length(x)
    for ii=1:length(x)
        Ks1(1,ii) = cov_func(x(ii,:)', x(i,:)', p1, p2);
        Ks2(ii,1) = cov_func(x(i,:)', x(ii,:)', p1, p2);
    end
    fhat(1,i) = Ks1*inv(K+sigma*eye(size(K)))*y;
    
    Kss = cov_func(x(i,:)', x(i,:)', p1, p2);
    
    cov(i) = Kss - Ks1*inv(K+sigma*eye(size(K))) * Ks2; 
end
%%
error = LOOCV(fhat', y)

%%
close all
plot(fhat); hold on; plot(y)

figure; plot(cov)