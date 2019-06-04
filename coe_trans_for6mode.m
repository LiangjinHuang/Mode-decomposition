function result=coe_trans_for6mode(amp,deta_p)

result=zeros(1,11);
amp=amp/sqrt(sum(abs(amp).^2 ));
result(1)=amp(1);
result(3)=amp(2)*exp(1i*deta_p(1));
result(5)=amp(3)*exp(1i*deta_p(2));
result(7)=amp(4)*exp(1i*deta_p(3));
result(9)=amp(5)*exp(1i*deta_p(4));
result(11)=amp(6)*exp(1i*deta_p(5));
end
