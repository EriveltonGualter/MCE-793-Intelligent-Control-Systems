% Function to train the network
%Created By : Anoop.V.S & Lekshmi B G
%Created On : 18-09-2013
%Description : Function to train the network
 
function [errorValue delta_W1 delta_W2] = trainNeuralNet(Input, Output, W1, W2, delta_W1, delta_W2)

    Output_of_InputLayer = Input;
    Input_of_HiddenLayer = W1' * Output_of_InputLayer;
       
    Output_of_HiddenLayer = 1./(1+exp(-Input_of_HiddenLayer));
    Input_of_OutputLayer = W2'*Output_of_HiddenLayer;

    clear m n;

    [m n] = size(Input_of_OutputLayer);

    Output_of_OutputLayer = 1./(1+exp(-Input_of_OutputLayer));

    difference = Output - Output_of_OutputLayer;
    square = difference.*difference;
    errorValue = sqrt(sum(square(:)));

    clear m n

    [n a] = size(Output);
 
    for i = 1 : n
        for j = 1 : a
            d(i,j) =(Output(i,j)-Output_of_OutputLayer(i,j))*Output_of_OutputLayer(i,j)*(1-Output_of_OutputLayer(i,j));
        end
    end

    %Now, calculate the Y matrix

    Y = Output_of_HiddenLayer * d'; %STEP 11

    if nargin == 4
        delta_W2=zeros(size(W2));
        delta_W1=zeros(size(W1));
    end

    %Initializing eetta with 0.6 and alpha with 1
    etta=100;alpha=1;
    
    %Calculating delta W
    delta_W2= alpha.*delta_W2 + etta.*Y;%STEP 12

    %STEP 13
    %Calculating error matrix
    error = W2*d;

    %Calculating d*
    clear m n
    [m n] = size(error);
    
    for i = 1 : m
        for j = 1 :n
            d_star(i,j)= error(i,j)*Output_of_HiddenLayer(i,j)*(1-Output_of_HiddenLayer(i,j));
        end
    end

    %Now find matrix, X (Input * transpose of d_star)

    X = Input * d_star';

    %STEP 14

    %Calculating delta V

    delta_W1=alpha*delta_W1+etta*X;
 
end