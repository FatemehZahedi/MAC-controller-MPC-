clear
close all
clc
% q=100; V=100; Cas=.0882; dH=2e5; ro=1e3; Cp=1; roc=1e3; Cpc=1; qc=100; ha=7e5; Ts=441.2; K0=7.2e10; J=1e4; Ks=K0*exp(-J/Ts); Ca0=1; T0=350; Tc0=350; Ks_=K0*(exp(-J/Ts))*(J/(Ts^2));
% a11=-q/V-Ks;
% a12=-Cas*Ks_;
% a21=-(dH/(ro*Cp))*Ks;
% a22=-q/V+(-dH*Cas/(ro*Cp))*Ks_+(-roc*Cpc/(ro*Cp*V))*qc+(roc*Cpc/(ro*Cp*V))*qc*exp(-ha/(qc*ro*Cp));
% b11=(Ca0-Cas)/V;
% b12=0;
% b21=(T0-Ts)/V;
% b22=((roc*Cpc)/(ro*Cp*V))*(Tc0-Ts)*(qc*(-exp(-ha/(qc*roc*Cpc))*(ha/((qc^2)*roc*Cpc)))+(1-exp(-ha/(qc*roc*Cpc))));
% A=[a11 a12; a21 a22];
% B=[b11 b12; b21 b22];
% C=[0 1];
% D=[0 0];
% [a1,b1]=ss2tf(A,B,C,D,1);
% G1=tf(a1,b1);
% [a2,b2]=ss2tf(A,B,C,D,2);
[a1,b1,a2,b2]=Inputsys(1);
G1=tf(a1,b1);
G2=tf(a2,b2);
figure(1);
step(G1);
hold on
[y0,t]=step(G1);
ts_matlab1=t(20,1)-t(19,1);  %Ts=0.0305  Ineligible
step(G2,'r');
[y00,t]=step(G2);
ts_matlab2=t(20,1)-t(19,1);   %Ts=0.0305  Ineligible
sys_info = stepinfo(G1);
tsc1 = sys_info.SettlingTime;  %ts1= 2.6054
sys_info = stepinfo(G2);
tsc2 = sys_info.SettlingTime; %ts2=2.5427
%  ==>  Ts=0.25;
Ts=0.25;
Gd1=c2d(G1,0.25,'impulse');
figure(2);
step(Gd1);
hold on
sys_info1 = stepinfo(Gd1);
ts1 = sys_info1.SettlingTime;
tr1=sys_info1.RiseTime;
N1=ts1/Ts;
P1=tr1/Ts;
Gd2=c2d(G2,0.25,'impulse');
step(Gd2,'r');
sys_info2 = stepinfo(Gd2);
ts2 = sys_info2.SettlingTime;
tr2=sys_info2.RiseTime;
N2=ts2/Ts;
P2=tr2/Ts;
if N2~=N1
    N2=N1;
end
if P2~=P1
    P2=P1;
end
% In this method M=P
M1=P1; M2=P2;
 [num1,den1]=tfdata(Gd1,'v');
 [num2,den2]=tfdata(Gd2,'v');
%% Design of MAC.............................................................................................................................................................
t = 0:0.25:30;
[h1,t1] =lsim(Gd1,[1 zeros(1,120)],t);
figure(3);
impulse(Gd1,t);
[h2,t2] =lsim(Gd2,[1 zeros(1,120)],t);
figure(4);
impulse(Gd2,t);
%...............Toeplitz Matrix............................................
b01=zeros(1,P1);
b01(1,1)=h1(2);
a01=h1(2:P1+1);
H1=toeplitz(a01,b01);
b02=zeros(1,P2);
b02(1,1)=h2(2);
a02=h2(2:P2+1);
H2=toeplitz(a02,b02);
H=[H1 H2];
%.....................Hankek Matrix........................................
c1=h1(3:P1+2);
r1=[h1(P1+2:N1+1)' zeros(1,P1-1)];
H1_=hankel(c1,r1);
c2=h2(3:P2+2);
r2=[h2(P2+2:N2+1)' zeros(1,P2-1)];
H2_=hankel(c2,r2);
H_=[H1_ H2_];
%........................................................................................
gamma=1/2;
% gain_DC1=(num1(1)+num1(2)+num1(3))/(den1(1)+den1(2)+den1(3));
% gain_DC2=(num2(1)+num2(2)+num2(3))/(den2(1)+den2(2)+den2(3));
gain_DC1=(num1(1)+num1(2))/(den1(1)+den1(2));
gain_DC2=(num2(1)+num2(2))/(den2(1)+den2(2));
% gain_DC1=(-0.0912 - 0.02606)/(1-0.9801+0.1035);
% gain_DC2=(-0.09053 + 0.04442 )/(1-0.9801+0.1035);
Q=eye(P1);
R1=gamma*(gain_DC1)^2*eye(M1);
R2=gamma*(gain_DC2)^2*eye(M2);
R=[R1 zeros(M1); zeros(M2) R2];
Kmac=((H'*Q*H)+R)\(H'*Q);

alpha=0.5;
y(1)=0;
r=ones(150+P1+1,1);
U=zeros(M1+M2,1);
U_=zeros(N1+N2-2,1);
ypast=zeros(150+P1,1);
ym=zeros(P1,1);
y_1=0; y_2=0;
d(1)=0;

for t=1:150
%     y(t)=-0.0912*U(1)-0.09053*U(M1+1)-0.02606*U_(1)+0.04442*U_(N1)+0.9801*y_1-0.1035*y_2;
% y(t)=num1(1)*U(1)+num2(1)*U(M1+1)+num1(2)*U_(1)+num2(2)*U_(N1)+num1(3)*U_(2)+num2(3)*U(N1+1)-den1(2)*y_1-den1(3)*y_2;
%  y(t)=num1(1)*U(1)+num2(1)*U(M1+1)+num1(2)*U_(1)+num2(2)*U_(N1)-den1(2)*y_1;
 y(t+1)=num1(1)*U(1)+num2(1)*U(M1+1)+num1(2)*U_(2)+num2(2)*U_(N1+1)-den1(2)*y_1;
    %y_2=y_1;
    y_1=y(t+1);
    yd(t+1)=y(t+1);
    for i=1:P1
        yd(t+i+1)=alpha*yd(t+i)+(1-alpha)*r(t+i+1);  % Programmed
        %yd(t+i)=alpha*yd(t+i-1)+(1-alpha)*r(t);   % Unprogrammed
    end
    
    d(t+1)=y(t+1)-ym(1);
     E=(yd(t+2:t+P1+1)'-ypast(t+1:t+P1)-d(t+1)*ones(P1,1));
     U=Kmac*E;
     ym=H*U+H_*U_;
     U_(2:N1-1)=U_(1:N1-2);
     U_(N1+1:N1+N2-2)=U_(N1:N1+N2-3);
     U_(1)=U(1); U_(N1)=U(M1+1);
     u1(t)=U(1);
     u2(t)=U(M1+1);
     ypast(t+1:t+P1)=H_*U_;
end 
 figure(5);
 plot(r,'b');
 hold on
 plot(y,'r');
 %axis([0 150 0 3]);
 figure(6);
 plot(u1);
 hold on
 plot(u2,'r');
%  %.................................................................................................................
%  %...................square wave..........................................
% %  [r,t] = gensig('square',8,636,0.1);
% %  r=0.1*r;
% %  for l = 1:length(r)
% % if(r(l)==0)
% % r(l) = -0.1;
% % end
% % end
% %r1=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -0.01 -0.01 -0.01 -0.01 -0.01 -0.01 -0.01 -0.01 -0.01 -0.01 -0.01 -0.01 -0.01 -0.01 -0.01 -0.01 -0.01 -0.01 -0.01 -0.01 -0.01 -0.01 -0.01 -0.01]';
% %r1=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -0.02 -0.02 -0.02 -0.02 -0.02 -0.02 -0.02 -0.02 -0.02 -0.02 -0.02 -0.02 -0.02 -0.02 -0.02 -0.02 -0.02 -0.02 -0.02 -0.02 -0.02 -0.02 -0.02 -0.02]';
% r1=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3]';
% r=[r1; r1; r1; r1; r1; r1; r1; r1; r1; r1; r1; r1; r1; r1; r1];
%  U=zeros(M1+M2,1);
% U_=zeros(N1+N2-2,1);
%  ypast=zeros(150+P1,1);
% ym=zeros(P1,1);
% %d(1)=0;
% y_1=0; y_2=0;
% 
% for t=1:150
%     y(t)=-0.0912*U(1)-0.09053*U(M1+1)-0.02606*U_(1)+0.04442*U_(N1)+0.9801*y_1-0.1035*y_2;
%     y_2=y_1;
%     y_1=y(t);
%     yd(t)=y(t);
%     for i=1:P1
%         yd(t+i)=alpha*yd(t+i-1)+(1-alpha)*r(t+i);  % Programmed
%         %yd(t+i)=alpha*yd(t+i-1)+(1-alpha)*r(t);   % Unprogrammed
%     end
%     
%     d(t)=y(t)-ym(1);
%      E=(yd(t+1:t+P1)'-ypast(t+1:t+P1)-d(t)*ones(P1,1));
%      U=Kmac*E;
%      ym=H*U+H_*U_;
%      U_(2:N1-1)=U_(1:N1-2);
%      U_(N1+1:N1+N2-2)=U_(N1:N1+N2-3);
%      U_(1)=U(1); U_(N1)=U(M1+1);
%      u(t)=U(1);
%      ypast(t+1:t+P1)=H_*U_;
% end 
%  figure(9);
%  plot(r,'b');
%  hold on
%  plot(y,'m');
%  axis([0 150 -2 6]);
%  figure(10);
%  plot(u);
%   %............................sine wave...........................................................   
%    [r,t] = gensig('sine',8,636,0.1);
%  U=zeros(M1+M2,1);
% U_=zeros(N1+N2-2,1);
%  ypast=zeros(150+P1,1);
% ym=zeros(P1,1);
% %d(1)=0;
% y_1=0; y_2=0;
% 
% for t=1:150
%     y(t)=-0.0912*U(1)-0.09053*U(M1+1)-0.02606*U_(1)+0.04442*U_(N1)+0.9801*y_1-0.1035*y_2;
%     y_2=y_1;
%     y_1=y(t);
%     yd(t)=y(t);
%     for i=1:P1
%         yd(t+i)=alpha*yd(t+i-1)+(1-alpha)*r(t+i);  % Programmed
%         %yd(t+i)=alpha*yd(t+i-1)+(1-alpha)*r(t);   % Unprogrammed
%     end
%     
%     d(t)=y(t)-ym(1);
%      E=(yd(t+1:t+P1)'-ypast(t+1:t+P1)-d(t)*ones(P1,1));
%      U=Kmac*E;
%      ym=H*U+H_*U_;
%      U_(2:N1-1)=U_(1:N1-2);
%      U_(N1+1:N1+N2-2)=U_(N1:N1+N2-3);
%      U_(1)=U(1); U_(N1)=U(M1+1);
%      u(t)=U(1);
%      ypast(t+1:t+P1)=H_*U_;
% end 
%  figure(11);
%  plot(r,'b');
%  hold on
%  plot(y,'m');
%  axis([0 150 -1 1]);
%  figure(12);
%  plot(u);
%  %.............................................. Disturbance............................................................
%   r=ones(180+P1+1,1);
% U=zeros(M1+M2,1);
% U_=zeros(N1+N2-2,1);
% ypast=zeros(180+P1,1);
% ym=zeros(P1,1);
% %d(1)=0;
% y_1=0; y_2=0;
% dist=zeros(180+P1,1);
% dist(16:20)=ones(5,1); dist(100:104)=ones(5,1); dist(65:69)=ones(5,1); 
% 
% for t=1:180
%     y(t)=-0.0912*U(1)-0.09053*U(M1+1)-0.02606*U_(1)+0.04442*U_(N1)+0.9801*y_1-0.1035*y_2+dist(t);
%     y_2=y_1;
%     y_1=y(t);
%     yd(t)=y(t);
%     for i=1:P1
%         yd(t+i)=alpha*yd(t+i-1)+(1-alpha)*r(t+i);  % Programmed
%         %yd(t+i)=alpha*yd(t+i-1)+(1-alpha)*r(t);   % Unprogrammed
%     end
%     
%     d(t)=y(t)-ym(1);
%      E=(yd(t+1:t+P1)'-ypast(t+1:t+P1)-d(t)*ones(P1,1));
%      U=Kmac*E;
%      ym=H*U+H_*U_;
%      U_(2:N1-1)=U_(1:N1-2);
%      U_(N1+1:N1+N2-2)=U_(N1:N1+N2-3);
%      U_(1)=U(1); U_(N1)=U(M1+1);
%      u(t)=U(1);
%      ypast(t+1:t+P1)=H_*U_;
% end 
%  figure(13);
%  plot(r,'b');
%  hold on
%  plot(y,'m');
%  axis([0 180 0 4]);
%  figure(14);
%  plot(u);
% %................................noise...........................................................................
%  r=ones(180+P1+1,1);
% U=zeros(M1+M2,1);
% U_=zeros(N1+N2-2,1);
% ypast=zeros(180+P1,1);
% ym=zeros(P1,1);
% %d(1)=0;
% y_1=0; y_2=0;
% noise=0.1*rand(180+P1,1); 
% 
% for t=1:180
%     y(t)=-0.0912*U(1)-0.09053*U(M1+1)-0.02606*U_(1)+0.04442*U_(N1)+0.9801*y_1-0.1035*y_2+noise(t);
%     y_2=y_1;
%     y_1=y(t);
%     yd(t)=y(t);
%     for i=1:P1
%         yd(t+i)=alpha*yd(t+i-1)+(1-alpha)*r(t+i);  % Programmed
%         %yd(t+i)=alpha*yd(t+i-1)+(1-alpha)*r(t);   % Unprogrammed
%     end
%     
%     d(t)=y(t)-ym(1);
%      E=(yd(t+1:t+P1)'-ypast(t+1:t+P1)-d(t)*ones(P1,1));
%      U=Kmac*E;
%      ym=H*U+H_*U_;
%      U_(2:N1-1)=U_(1:N1-2);
%      U_(N1+1:N1+N2-2)=U_(N1:N1+N2-3);
%      U_(1)=U(1); U_(N1)=U(M1+1);
%      u(t)=U(1);
%      ypast(t+1:t+P1)=H_*U_;
% end 
%  figure(15);
%  plot(r,'b');
%  hold on
%  plot(y,'m');
%  %axis([0 180 0 4]);
%  figure(16);
%  plot(u);
%  %................................uncertainty.........................................................................
%  % 5%
%  [a1,b1]=ss2tf(A,B,C,D,1);
% G1=tf(a1,b1+[0 0.05*b1(2) 0]);
% [a2,b2]=ss2tf(A,B,C,D,2);
% G2=tf(a2,b2+[0 0.05*b2(2) 0]);
% Gd1=c2d(G1,0.1,'impulse');
% % figure(3);
% % step(Gd1);
% sys_info1 = stepinfo(Gd1);
% ts1 = sys_info1.SettlingTime;
% tr1=sys_info1.RiseTime;
% N1=floor(ts1/Ts);
% P1=floor(tr1/Ts);
% Gd2=c2d(G2,0.1,'impulse');
% % figure(4);
% % step(Gd2);
% sys_info2 = stepinfo(Gd2);
% ts2 = sys_info2.SettlingTime;
% tr2=sys_info2.RiseTime;
% N2=floor(ts2/Ts);
% P2=P1;
% %P2=floor(tr2/Ts);
% % In this method M=P
% M1=P1; M2=P2;
% 
%  t = 0:0.1:20;
% [h1,t1] =lsim(Gd1,[1 zeros(1,200)],t);
% % figure(5);
% % impulse(Gd1,t);
% [h2,t2] =lsim(Gd2,[1 zeros(1,200)],t);
% % figure(6);
% % impulse(Gd2,t);
% %...............Toeplitz Matrix............................................
% b1=zeros(1,P1);
% b1(1,1)=h1(2);
% a1=h1(2:P1+1);
% H1=toeplitz(a1,b1);
% b2=zeros(1,P2);
% b2(1,1)=h2(2);
% a2=h2(2:P2+1);
% H2=toeplitz(a2,b2);
% H=[H1 H2];
% %.....................Hankek Matrix........................................
% c1=h1(3:P1+2);
% r1=[h1(P1+2:N1+1)' zeros(1,P1-1)];
% H1_=hankel(c1,r1);
% c2=h2(3:P2+2);
% r2=[h2(P2+2:N2+1)' zeros(1,P2-1)];
% H2_=hankel(c2,r2);
% H_=[H1_ H2_];
% %........................................................................................
% gamma=1;
% % gain_DC1=(-0.0912 - 0.02606)/(1-0.9801+0.1035);
% % gain_DC2=(-0.09053 + 0.04442 )/(1-0.9801+0.1035);
% gain_DC1=(-0.0912 - 0.02198 - 1.831e-18)/( 1 - 0.9734  + 0.09243);
% gain_DC2=( -0.09053 + 0.04598)/( 1 - 0.9734  + 0.09243);
% Q=7*eye(P1);
% R1=gamma*(gain_DC1)^2*eye(M1);
% R2=gamma*(gain_DC2)^2*eye(M2);
% R=[R1 zeros(M1); zeros(M2) R2];
% Kmac=((H'*Q*H)+R)\(H'*Q);
% 
% alpha=0.5;
% r=ones(150+P1+1,1);
% U=zeros(M1+M2,1);
% U_=zeros(N1+N2-2,1);
% ypast=zeros(150+P1,1);
% ym=zeros(P1,1);
% y_1=0; y_2=0;
% 
% for t=1:150
%     y(t)=-0.0912*U(1)-0.09053*U(M1+1)-0.02198*U_(1)-(1.831e-18)*U_(2)+0.04442*U_(N1)+0.9734*y_1-0.09243*y_2;
%     y_2=y_1;
%     y_1=y(t);
%     yd(t)=y(t);
%     for i=1:P1
%         yd(t+i)=alpha*yd(t+i-1)+(1-alpha)*r(t+i);  % Programmed
%         %yd(t+i)=alpha*yd(t+i-1)+(1-alpha)*r(t);   % Unprogrammed
%     end
%     
%     d(t)=y(t)-ym(1);
%      E=(yd(t+1:t+P1)'-ypast(t+1:t+P1)-d(t)*ones(P1,1));
%      U=Kmac*E;
%      ym=H*U+H_*U_;
%      U_(2:N1-1)=U_(1:N1-2);
%      U_(N1+1:N1+N2-2)=U_(N1:N1+N2-3);
%      U_(1)=U(1); U_(N1)=U(M1+1);
%      u(t)=U(1);
%      ypast(t+1:t+P1)=H_*U_;
% end 
%  figure(17);
%  plot(r,'b');
%  hold on
%  plot(y,'m');
%  axis([0 150 0 1.5]);
%  figure(18);
%  plot(u);
% %.................................................................................................
%  % 10%
%  [a1,b1]=ss2tf(A,B,C,D,1);
% G1=tf(a1,b1+[0 0.1*b1(2) 0]);
% [a2,b2]=ss2tf(A,B,C,D,2);
% G2=tf(a2,b2+[0 0.1*b2(2) 0]);
% Gd1=c2d(G1,0.1,'impulse');
% % figure(3);
% % step(Gd1);
% sys_info1 = stepinfo(Gd1);
% ts1 = sys_info1.SettlingTime;
% tr1=sys_info1.RiseTime;
% N1=floor(ts1/Ts);
% P1=floor(tr1/Ts);
% Gd2=c2d(G2,0.1,'impulse');
% % figure(4);
% % step(Gd2);
% sys_info2 = stepinfo(Gd2);
% ts2 = sys_info2.SettlingTime;
% tr2=sys_info2.RiseTime;
% N2=floor(ts2/Ts);
% P2=P1;
% %P2=floor(tr2/Ts);
% % In this method M=P
% M1=P1; M2=P2;
% 
%  t = 0:0.1:20;
% [h1,t1] =lsim(Gd1,[1 zeros(1,200)],t);
% % figure(5);
% % impulse(Gd1,t);
% [h2,t2] =lsim(Gd2,[1 zeros(1,200)],t);
% % figure(6);
% % impulse(Gd2,t);
% %...............Toeplitz Matrix............................................
% b1=zeros(1,P1);
% b1(1,1)=h1(2);
% a1=h1(2:P1+1);
% H1=toeplitz(a1,b1);
% b2=zeros(1,P2);
% b2(1,1)=h2(2);
% a2=h2(2:P2+1);
% H2=toeplitz(a2,b2);
% H=[H1 H2];
% %.....................Hankek Matrix........................................
% c1=h1(3:P1+2);
% r1=[h1(P1+2:N1+1)' zeros(1,P1-1)];
% H1_=hankel(c1,r1);
% c2=h2(3:P2+2);
% r2=[h2(P2+2:N2+1)' zeros(1,P2-1)];
% H2_=hankel(c2,r2);
% H_=[H1_ H2_];
% %........................................................................................
% gamma=1;
% % gain_DC1=(-0.0912 - 0.02606)/(1-0.9801+0.1035);
% % gain_DC2=(-0.09053 + 0.04442 )/(1-0.9801+0.1035);
% gain_DC1=(-0.0912- 0.01814)/( 1- 0.9677 + 0.08252);
% gain_DC2=(-0.09053 + 0.04745)/( 1 - 0.9677 + 0.08252);
% Q=7*eye(P1);
% R1=gamma*(gain_DC1)^2*eye(M1);
% R2=gamma*(gain_DC2)^2*eye(M2);
% R=[R1 zeros(M1); zeros(M2) R2];
% Kmac=((H'*Q*H)+R)\(H'*Q);
% 
% alpha=0.5;
% r=ones(150+P1+1,1);
% U=zeros(M1+M2,1);
% U_=zeros(N1+N2-2,1);
% ypast=zeros(150+P1,1);
% ym=zeros(P1,1);
% y_1=0; y_2=0;
% 
% for t=1:150
%     y(t)=-0.0912*U(1)-0.09053*U(M1+1)- 0.01814*U_(1)+ 0.04745*U_(N1)+ 0.9677*y_1-0.08252*y_2;
%     y_2=y_1;
%     y_1=y(t);
%     yd(t)=y(t);
%     for i=1:P1
%         yd(t+i)=alpha*yd(t+i-1)+(1-alpha)*r(t+i);  % Programmed
%         %yd(t+i)=alpha*yd(t+i-1)+(1-alpha)*r(t);   % Unprogrammed
%     end
%     
%     d(t)=y(t)-ym(1);
%      E=(yd(t+1:t+P1)'-ypast(t+1:t+P1)-d(t)*ones(P1,1));
%      U=Kmac*E;
%      ym=H*U+H_*U_;
%      U_(2:N1-1)=U_(1:N1-2);
%      U_(N1+1:N1+N2-2)=U_(N1:N1+N2-3);
%      U_(1)=U(1); U_(N1)=U(M1+1);
%      u(t)=U(1);
%      ypast(t+1:t+P1)=H_*U_;
% end 
%  figure(19);
%  plot(r,'b');
%  hold on
%  plot(y,'m');
%  axis([0 150 0 1.5]);
%  figure(20);
%  plot(u);
%  %.............................................................................................................
%   % 30%
%  [a1,b1]=ss2tf(A,B,C,D,1);
% G1=tf(a1,b1+[0 0.3*b1(2) 0]);
% [a2,b2]=ss2tf(A,B,C,D,2);
% G2=tf(a2,b2+[0 0.3*b2(2) 0]);
% Gd1=c2d(G1,0.1,'impulse');
% % figure(3);
% % step(Gd1);
% sys_info1 = stepinfo(Gd1);
% ts1 = sys_info1.SettlingTime;
% tr1=sys_info1.RiseTime;
% N1=floor(ts1/Ts);
% P1=floor(tr1/Ts);
% Gd2=c2d(G2,0.1,'impulse');
% % figure(4);
% % step(Gd2);
% sys_info2 = stepinfo(Gd2);
% ts2 = sys_info2.SettlingTime;
% tr2=sys_info2.RiseTime;
% N2=floor(ts2/Ts);
% P2=P1;
% %P2=floor(tr2/Ts);
% % In this method M=P
% M1=P1; M2=P2;
% 
%  t = 0:0.1:20;
% [h1,t1] =lsim(Gd1,[1 zeros(1,200)],t);
% % figure(5);
% % impulse(Gd1,t);
% [h2,t2] =lsim(Gd2,[1 zeros(1,200)],t);
% % figure(6);
% % impulse(Gd2,t);
% %...............Toeplitz Matrix............................................
% b1=zeros(1,P1);
% b1(1,1)=h1(2);
% a1=h1(2:P1+1);
% H1=toeplitz(a1,b1);
% b2=zeros(1,P2);
% b2(1,1)=h2(2);
% a2=h2(2:P2+1);
% H2=toeplitz(a2,b2);
% H=[H1 H2];
% %.....................Hankek Matrix........................................
% c1=h1(3:P1+2);
% r1=[h1(P1+2:N1+1)' zeros(1,P1-1)];
% H1_=hankel(c1,r1);
% c2=h2(3:P2+2);
% r2=[h2(P2+2:N2+1)' zeros(1,P2-1)];
% H2_=hankel(c2,r2);
% H_=[H1_ H2_];
% %........................................................................................
% gamma=1;
% % gain_DC1=(-0.0912 - 0.02606)/(1-0.9801+0.1035);
% % gain_DC2=(-0.09053 + 0.04442 )/(1-0.9801+0.1035);
% gain_DC1=(-0.0912  - 0.004749  - 1.648e-20)/(   1- 0.952  + 0.05243);
% gain_DC2=( -0.09053  + 0.05258 )/( 1 - 0.952  + 0.05243);
% Q=7*eye(P1);
% R1=gamma*(gain_DC1)^2*eye(M1);
% R2=gamma*(gain_DC2)^2*eye(M2);
% R=[R1 zeros(M1); zeros(M2) R2];
% Kmac=((H'*Q*H)+R)\(H'*Q);
% 
% alpha=0.5;
% r=ones(150+P1+1,1);
% U=zeros(M1+M2,1);
% U_=zeros(N1+N2-2,1);
% ypast=zeros(150+P1,1);
% ym=zeros(P1,1);
% y_1=0; y_2=0;
% 
% for t=1:150
%     y(t)=-0.0912*U(1)-0.09053*U(M1+1)-  0.004749*U_(1)-(1.648e-20)*U_(2)+ 0.05258*U_(N1)+ 0.952*y_1- 0.05243*y_2;
%     y_2=y_1;
%     y_1=y(t);
%     yd(t)=y(t);
%     for i=1:P1
%         yd(t+i)=alpha*yd(t+i-1)+(1-alpha)*r(t+i);  % Programmed
%         %yd(t+i)=alpha*yd(t+i-1)+(1-alpha)*r(t);   % Unprogrammed
%     end
%     
%     d(t)=y(t)-ym(1);
%      E=(yd(t+1:t+P1)'-ypast(t+1:t+P1)-d(t)*ones(P1,1));
%      U=Kmac*E;
%      ym=H*U+H_*U_;
%      U_(2:N1-1)=U_(1:N1-2);
%      U_(N1+1:N1+N2-2)=U_(N1:N1+N2-3);
%      U_(1)=U(1); U_(N1)=U(M1+1);
%      u(t)=U(1);
%      ypast(t+1:t+P1)=H_*U_;
% end 
%  figure(21);
%  plot(r,'b');
%  hold on
%  plot(y,'m');
%  axis([0 150 0 1.5]);
%  figure(22);
%  plot(u);
% 
