function [ result ] = cost_func( I,I0,dx,dy)
%COST_FUNC Summary of this function goes here
%   

I=I/(sum(sum(I)));
[a b]=size(I);
meanI=sum(sum(I))/(a*b);
meanI0=sum(sum(I0))/(a*b);
deltaI=I-meanI;
deltaI0=I0-meanI0;
result=abs(sum(sum(deltaI.*deltaI0))/(sqrt(sum(sum(deltaI.^2))*sum(sum(deltaI0.^2)))));
end

