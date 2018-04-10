close all
clear
clc


In = 6;
Hid= 10;
Out= 2;

W1=rand(Hid,In);
W2=rand(Hid,Out);
grad=zeros(2*Hid,In)+2;
New_Error=20;            

x = 0:pi/2:4*pi;
tau = sin(x);
%  x=[q(:,1) q(:,2) qDot(:,1) qDot(:,2) qDDot(:,1) qDDot(:,2)];

while New_Error>0.5
Y_Hat = zeros(size(tau));
clear a z

    for i=1:Hid
    a(:,i)=x(:,1)*W1(i,1)+x(:,2)*W1(i,2)+x(:,3)*W1(i,3)+x(:,4)*W1(i,4)+x(:,5)*W1(i,5)+x(:,6)*W1(i,6);
    z(:,i)=1./(1+exp(-a(:,i)));
    Y_Hat(:,1)=z(:,i)*W2(i,1);
    end

    Delta(:,1)=Y_Hat(:,1)-tau(:,1);

    for i=1:Hid
    Delta1(:,i)=exp(-a(:,i))./((1+exp(-a(:,i))).^2).*Delta(:,1).*W2(i,1);
    end

    for i=1:Hid
    dEndw2y1(:,i)=Delta(:,1).*z(:,i);
    dEndw1x1(:,i)=Delta1(:,i).*x(:,i);
    dEndw1x3(:,i)=Delta1(:,i).*x(:,1);
    dEndw1x5(:,i)=Delta1(:,i).*x(:,2);
    end

    grad1=[sum(dEndw1x3);sum(dEndw1x5)]';
    grad2=sum(dEndw1x1)';
%     grad=[sum(dEndw1x1(:,1)) sum(dEndw1x2(:,1)) sum(dEndw1x3(:,1)) sum(dEndw1x4(:,1)) sum(dEndw1x5(:,1)) sum(dEndw1x6(:,1))];
% 
%     for i=2:Hid
%     grad=vertcat(grad,[sum(dEndw1x1(:,i)) sum(dEndw1x2(:,i)) sum(dEndw1x3(:,i)) sum(dEndw1x4(:,i)) sum(dEndw1x5(:,i)) sum(dEndw1x6(:,i))]);
%     end
% 
%     for i=1:Hid
%     grad=vertcat(grad,[sum(dEndw2y1(:,i)) sum(dEndw2y2(:,i)) 0 0 0 0]);
%     end

eta=0.0001;
% W1=W1-eta*grad(1:Hid,:);
% W2=W2-eta*grad(Hid+1:2*Hid,1:2);
eta=0.001;
W1=W1-eta*grad1;
W2=W2-eta*grad2;
New_Error=sum(sqrt((sum(tau-Y_Hat).^2)/length(tau)))/2

end

% function hw5(Input, Output)
%     %Assigning the number of hidden neurons in hidden layer
%     m = 5;
% 
%     %Find the size of Input and Output Vectors
%     [l,b] = size(Input);
%     [n,a] = size(Output);
%     
%     %Initialize the weight matrices with random weights
%     W1 = rand(l,m); % Weight matrix from Input to Hidden
%     W2 = rand(m,n); % Weight matrix from Hidden to Output
%     
%     while errorValue > 1
%         for i=1:m
%             Input_of_HiddenLayer = W1' * Input;
%             Input_of_OutputLayer = W2'*1./(1+exp(-Input_of_HiddenLayer));
%         end
%     end
%     
%     Delta = Input_of_OutputLayer-Input;
%     
%     for i=1:m
%         Delta1 = exp(-Input_of_HiddenLayer)./((1+exp(-Input_of_HiddenLayer)).^2).*Delta(:,1).*W2(i,1);
%     end
%     
%     errorValue = sqrt(sum(sqrt(Output - Output_of_OutputLayer)));
% end