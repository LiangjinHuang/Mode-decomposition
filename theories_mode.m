
clear;
clc;
clear;
clc;
close all;
Z0=377.3233;
pix=4.4e-6;   % the size of pixel
dx=pix; 
dy=dx;
N=220;
M=50; % magnification factor
lambda=1.073e-6;
k0=2*pi/lambda; 
a=11.8e-6;  % core radius of the investigated fiber
n1=1.458;   % core refractive index
NA=0.064;
n2=sqrt(n1^2-NA^2);   % cladding refractive index

n_01=1.457727432734001;
n_11=1.457323122892004;
n_21=1.456825922970006;
n_02=1.456711429140005; % refractive index of eigenmodes

n=[n_01 n_11 n_21 n_02];
U1=a*k0*sqrt(n1^2-n.^2);
W1=a*k0*sqrt(n.^2-n2^2);
fai=0*pi;
%% Acquiring the eigenmodes 
dx_0=pix/M;
x=-N*dx_0:dx_0:N*dx_0;
y=x;
leng=length(x);                  % the amount of points in x direction
[X,Y]=meshgrid(x,y);             
[THETA,RHO] = cart2pol(X,Y);     
mode1=zeros(size(RHO));               %Initializing the electric field
mode2=mode1;  % LP11e
mode3=mode1;  % LP11o
mode4=mode1;  % LP21e
mode5=mode1;  % LP21o
mode6=mode1;  % LP02

orientation='even';              % even -> cos£¬odd -> sin
for ii=1:6
    if ii==1
        U=U1(ii);
        W=W1(ii);
        m=0;
    elseif ii==2
        U=U1(ii);
        W=W1(ii);
        m=1;
    elseif ii==3
        U=U1(ii-1);
        W=W1(ii-1);
        m=1;
        orientation='odd'; 
     elseif ii==4
        U=U1(ii-1);
        W=W1(ii-1);
        orientation='even'; 
        m=2;
    elseif ii==5
        U=U1(ii-2);
        W=W1(ii-2);
        orientation='odd';
        m=2;
    elseif ii==6
        U=U1(ii-2);
        W=W1(ii-2);
        orientation='even';
        m=0;
    elseif ii==7
        U=U1(ii-2);
        W=W1(ii-2);
        orientation='even';
        m=3;
    elseif ii==8
        U=U1(ii-3);
        W=W1(ii-3);
        orientation='odd';
        m=3;
    end
%% 3.1 R<a
[row col] = find(RHO <= a);  % finding the coordinates satisfying R<a, then assigning new value to them.
        for i=1:size(row)
            if strcmp(orientation,'even') 
                E(row(i), col(i)) =besselj(m,U/a*RHO(row(i), col(i))).*cos(m*THETA(row(i), col(i))+fai); 
            else  
                E(row(i), col(i)) =besselj(m,U/a*RHO(row(i), col(i))).*sin(m*THETA(row(i), col(i))+fai); 
            end
        end
%% 3.2 a<R<b
[row col] = find(RHO>a);    % finding the coordinates satisfying a<R<b, then assigning new value to them.
        for i=1:size(row)
            if strcmp(orientation,'even') 
                E(row(i), col(i)) =besselj(m,U)/besselk(m,W)*besselk(m,W/a*RHO(row(i), col(i))).*cos(m*THETA(row(i), col(i))+fai); 
            else  
                E(row(i), col(i)) =besselj(m,U)/besselk(m,W)*besselk(m,W/a*RHO(row(i), col(i))).*sin(m*THETA(row(i), col(i))+fai); 
            end
        end
        eval(strcat('mode',num2str(ii),'=E;'));
         
end

mode1=mode1./sqrt((sum(sum(abs(mode1).^2))));   % normalization
mode2=mode2./sqrt((sum(sum(abs(mode2).^2))));
mode3=mode3./sqrt((sum(sum(abs(mode3).^2))));
mode4=mode4./sqrt((sum(sum(abs(mode4).^2))));  
mode5=mode5./sqrt((sum(sum(abs(mode5).^2))));
mode6=mode6./sqrt((sum(sum(abs(mode6).^2))));

clearvars variables -except mode1 mode2 mode3 mode4 mode5 mode6  N dx dy dx_0