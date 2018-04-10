
clear al
close all
clc

x = 0:pi/20:4*pi;
in = x;sin(x);
ou = x;sin(x);

in = in;
ou = ou;

BackPropAlgo(in,ou)