function y=stochastic(m,n)
error(nargchk(1,2,nargin));
if nargin==1
    n=m;
end
y=rand(m,n);
[row col]=find(y>0.5);
len=length(row);
for i=1:len
    y(row(i),col(i))=1;
end
[row1 col1]=find(y<=0.5);
len1=length(row1);
for i=1:len1
    y(row1(i),col1(i))=-1;
end 
end