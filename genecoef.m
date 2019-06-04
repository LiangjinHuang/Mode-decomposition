
function [r1 r2 r3 r4 r5 r6 p2 p3 p4 p5 p6 ]=geneCoef
r1=rand(1);
r2=rand(1);
r3=rand(1);
r4=rand(1);
r5=rand(1);
r6=rand(1);
rt=r1+r2+r3+r4+r5+r6;
r1=r1/rt;
r2=r2/rt;
r3=r3/rt;
r4=r4/rt;
r5=r5/rt;
r6=r6/rt;
p2=2*pi*(rand(1)-0.5);
p3=2*pi*(rand(1)-0.5);
p4=2*pi*(rand(1)-0.5);
p5=2*pi*(rand(1)-0.5);
p6=2*pi*(rand(1)-0.5);
end

