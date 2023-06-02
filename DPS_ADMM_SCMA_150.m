function [SCMA_MU_CW_det]=DPS_ADMM_SCMA_150(H,y,J,rho,al,T,d_v,al_R,be_I)
I=eye(J*d_v)/(H'*H*J+rho*eye(J*d_v));
Ma=H'*y;
%Initializing
x_11=zeros(J*d_v,1)+zeros(J*d_v,1)*1i;
x_12=zeros(J*d_v,1)+zeros(J*d_v,1)*1i;
x_13=zeros(J*d_v,1)+zeros(J*d_v,1)*1i;
x_14=zeros(J*d_v,1)+zeros(J*d_v,1)*1i;
x_15=zeros(J*d_v,1)+zeros(J*d_v,1)*1i;
x_16=zeros(J*d_v,1)+zeros(J*d_v,1)*1i;
x_1_bar=zeros(J*d_v,1)+zeros(J*d_v,1)*1i;
x_0_bar=zeros(J*d_v,1)+zeros(J*d_v,1)*1i;
u=zeros(J*d_v,1)+zeros(J*d_v,1)*1i;
%SCMA_MU_CW_det=zeros(J*d_v,1)+zeros(J*d_v,1)*1i;

for t=1:T
   x_11=(rho/(rho+al))*(x_11+x_0_bar-u-x_1_bar);
   x_12=(rho/(rho+al))*(x_12+x_0_bar-u-x_1_bar);
   x_13=(rho/(rho+al))*(x_13+x_0_bar-u-x_1_bar);
   x_14=(rho/(rho+al))*(x_14+x_0_bar-u-x_1_bar);
   x_15=(rho/(rho+al))*(x_15+x_0_bar-u-x_1_bar);
   x_16=(rho/(rho+al))*(x_16+x_0_bar-u-x_1_bar);
   %Uesr-1
   x_11R=real(x_11);x_11I=imag(x_11);
   x_11R(x_11R<-al_R(1))=-al_R(1); x_11R(x_11R>al_R(1))=al_R(1);
   x_11I(x_11I<-be_I(1))=-be_I(1); x_11I(x_11I>be_I(1))=be_I(1);
   %User-2
   x_12R=real(x_12);x_12I=imag(x_12);
   x_12R(x_12R<-al_R(2))=-al_R(2); x_12R(x_12R>al_R(2))=al_R(2);
   x_12I(x_12I<-be_I(2))=-be_I(2);x_12I(x_12I>be_I(2))=be_I(2);
   %User-3
   x_13R=real(x_13);x_13I=imag(x_13);
   x_13R(x_13R<-al_R(3))=-al_R(3); x_13R(x_13R>al_R(3))=al_R(3);
   x_13I(x_13I<-be_I(3))=-be_I(3);x_13I(x_13I>be_I(3))=be_I(3);
   %user-4
   x_14R=real(x_14);x_14I=imag(x_14);
   x_14R(x_14R<-al_R(4))=-al_R(4); x_14R(x_14R>al_R(4))=al_R(4);
   x_14I(x_14I<-be_I(4))=-be_I(4);x_14I(x_14I>be_I(4))=be_I(4);
   %User-5
   x_15R=real(x_15);x_15I=imag(x_15);
   x_15R(x_15R<-al_R(5))=-al_R(5); x_15R(x_15R>al_R(5))=al_R(5);
   x_15I(x_15I<-be_I(5))=-be_I(5);x_15I(x_15I>be_I(5))=be_I(5);
   %User-6
   x_16R=real(x_16);x_16I=imag(x_16);
   x_16R(x_16R<-al_R(6))=-al_R(6); x_16R(x_16R>al_R(6))=al_R(6);
   x_16I(x_16I<-be_I(6))=-be_I(6);x_16I(x_16I>be_I(6))=be_I(6);
   x_11=x_11R+x_11I*1j;x_12=x_12R+x_12I*1j;x_13=x_13R+x_13I*1j;
   x_14=x_14R+x_14I*1j;x_15=x_15R+x_15I*1j;x_16=x_16R+x_16I*1j;  
    x_11(3:end)=0;x_12(1:2)=0;x_12(5:end)=0;x_13(1:4)=0;x_13(7:end)=0;
    x_14(1:6)=0;x_14(9:end)=0;x_15(1:8)=0;x_15(11:end)=0;x_16(1:10)=0;
   x_1_bar=sum(x_11+x_12+x_13+x_14+x_15+x_16)/J;  
   x_0_bar=I*(Ma+rho*(x_1_bar+u));
   u=u+(x_1_bar-x_0_bar);      
end
SCMA_MU_CW_det=J*x_0_bar;



