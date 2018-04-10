function k=kernal(x,q,M,h)
d=((x-q)'*M'*M*(x-q))^0.5;
k=exp(-((d^2)/h));
end