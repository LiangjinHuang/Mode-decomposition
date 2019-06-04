
theories_mode    %calculate the eigenmodes
x=-N*dx_0:dx_0:N*dx_0; 
y=x;
%% ------------------------Preparing measured beam intensity------------------------%

[r01,r11e,r11o,r21e,r21o, r02, dphi11e,dphi11o,dphi21e,dphi21o,dphi02]=genecoef; %randomly generating the mode coefficients
Er=sqrt(r01)*mode1 + sqrt(r11e)*exp(1j*dphi11e)*mode2 + sqrt(r11o)*exp(1j*dphi11o)*mode3+...
sqrt(r21e)*exp(1j*dphi21e)*mode4 + sqrt(r21o)*exp(1j*dphi21o)*mode5+sqrt(r02)*exp(1j*dphi02)*mode6; %simulating the measured near-field. 
I0=Er.*conj(Er);  %beam intensity
P=dx*dy*sum(sum(I0));
I0=I0/P;   % normalizating the measured beam intensity
%---------------------------------------------------------------------------%
NN=1;
loopnumber=1000;  % loopnumber. increaing this value can increase the accuracy but will cost more time.
Jnn=zeros(NN,loopnumber);
amp_coe=zeros(NN,6);
ttt=zeros(1,NN);
jjj=0;
while jjj<NN
 %%  generating the initial values 
 jjj=jjj+1 
amp1=rand;   
amp2=rand;
amp3=rand;
amp4=rand;    
amp5=rand;
amp6=rand;
deta_p12=rand;
deta_p13=rand;
deta_p14=rand;
deta_p15=rand;
deta_p16=rand;
amp0=[amp1,amp2,amp3,amp4,amp5,amp6];  
deta_p0=[deta_p12 deta_p13  deta_p14  deta_p15  deta_p16];       
%% Alternating correction
sign0=0;         
sign00=0;  
tic
while sign00<1   
    sign00=sign00+1
    
    amp=amp0;
    deta_p=deta_p0;    %Initialization
    sign0=0;

    amp_w=1000;  
    deta_p_w=1200;
    temp=0.042;   
    temp2=0.018;  % These four values can be adjusted. 
    
    iteration1=loopnumber;

    while sign0<1   
        sign0=sign0+1;
           sign=0;
            while sign<iteration1   
                sign=sign+1;
                pert=stochastic(1,11);     % perturbation
                d_amp=temp*pert(1:6);     
                d_phase=temp2*pert(7:11);    
                amp_f=amp+d_amp;  
                
                deta_p_f=deta_p+d_phase;
                coe_f=coe_trans_for6mode(amp_f,deta_p_f);
                I1=coe_f(1)*mode1+coe_f(3)*mode2+coe_f(5)*mode3+...
                    coe_f(7)*mode4+coe_f(9)*mode5+coe_f(11)*mode6;
                I1=n1/(2*Z0)*abs(I1).^2;

                amp_b=amp-d_amp;  
                deta_p_b=deta_p-d_phase;
                coe_b=coe_trans_for6mode(amp_b,deta_p_b);
                I2=coe_b(1)*mode1+coe_b(3)*mode2+coe_b(5)*mode3+...
                    coe_b(7)*mode4+coe_b(9)*mode5+coe_b(11)*mode6;
                I2=n1/(2*Z0)*abs(I2).^2;
               dJ(sign)=(cost_func(I1,I0,dx,dy)-cost_func(I2,I0,dx,dy)); 
                amp=amp+amp_w*dJ(sign).*d_amp;   
                deta_p=deta_p+deta_p_w*dJ(sign).*d_phase;
                coe=coe_trans_for6mode(amp,deta_p);
                I=coe(1)*mode1+coe(3)*mode2+coe(5)*mode3+...
                    coe(7)*mode4+coe(9)*mode5+coe(11)*mode6;
                I=n1/(2*Z0)*abs(I).^2;
                J(sign)=cost_func(I,I0,dx,dy);    
            end
toc

    end
    
J(end);
end

[amp_cons deta_p_cons]=inver_trans_for6mode(coe);
amp_cons =amp_cons/sqrt(sum(abs(amp_cons).^2 ));

amp_coe(jjj,:)=amp_cons;
phase_coe(jjj,:)=deta_p_cons;

end
amp_coe.^2
%% residual of measured and reconstructed beam intensity
I=I/(sum(sum(I))); %reconstructed beam intensity
I0=I0/(sum(sum(I0)));%measured beam intensity
residual_I=abs(I0-I);
P_residual=dx*dy*sum(sum(residual_I));
residual_error=P_residual;
%% plotting the measured, reconstructed and residual beam intensity
x=1e6*x;
y=1e6*y;
figure;
colormap('jet')
I0=I0/(max(max(I0)));
axes1 = axes('YDir','reverse','Layer','top','FontSize',14);
imagesc(x(1,:),y(1,:),I0);
title('Measured')
axis off

figure;
colormap('jet')
I=I/(max(max(I)));
axes1 = axes('YDir','reverse','Layer','top','FontSize',14);
imagesc(x(1,:),y(1,:),I);
title('Reconstructed')
axis off
colorbar;
residual_I=abs(I0-I);

figure;
colormap('jet')
axes1 = axes('YDir','reverse','Layer','top','FontSize',14);
imagesc(x(1,:),y(1,:),residual_I);
title('Residual')
axis off
colorbar;
caxis([0 1]);
%%
[amp_cons deta_p_cons]=inver_trans_for6mode(coe);
amp_cons =amp_cons/sqrt(sum(abs(amp_cons).^2 ));
J(end) %correlation between the measured and reconstructed beam intensity.
J_avr=sum(Jnn)/NN;