function [amp deta_p]=inver_trans_for6mode(coe)
amp(1,1)=abs(coe(1));
amp(1,2)=abs(coe(3));
amp(1,3)=abs(coe(5));
amp(1,4)=abs(coe(7));
amp(1,5)=abs(coe(9));
amp(1,6)=abs(coe(11));


deta_p(1)=angle(coe(3))-angle(coe(1));
deta_p(2)=angle(coe(5))-angle(coe(1));
deta_p(3)=angle(coe(7))-angle(coe(1));
deta_p(4)=angle(coe(9))-angle(coe(1));
deta_p(5)=angle(coe(11))-angle(coe(1));




end