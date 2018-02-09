% Erivelton Gualter dos Santos

%clc, clear all, close all

load('pendulumdata.mat'); % Load data 

g= 9.81;   % gravitational acceleration

t = [qDDOT];
phi = [3*T -3*g*cos(q)/2];

%Moore–Penrose inverse
PHI = inv(phi'*phi)*phi'; % or PHI =  pinv(phi); % Pseudoinverse function

WML = PHI*t;
%%

B = 0.01^-1;
M0 = [1 ; 1];

S0 = 10*eye(2,2);

for k=1:size(T)
    
    if k>1
        phi = vertcat(phi, [3*T(k) -3*g*cos(q(k))/2]); 
        Tk = vertcat(Tk, t(k));
    else
        phi = [3*t(k) -3*g*cos(q(k))/2];
        Tk = t(k);
    end
    
    SN(:,:,k) = inv( inv(S0) + B*phi'*phi );
    MN(k,:) = SN(:,:,k) * (inv(S0)*M0 + B*phi'*Tk);
end

plot(MN)
