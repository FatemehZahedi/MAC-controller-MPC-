% [a1,b1,a2,b2]=Inputsys(3);
% G1=tf(a1,b1);
% Gd1=c2d(G1,0.25,'impulse');
% figure(3);
% step(Gd1);
% sys_info1 = stepinfo(Gd1);
% ts1 = sys_info1.SettlingTime;
% tr1=sys_info1.RiseTime;
% Ts=0.25;
% N1=ts1/Ts;
% P1=tr1/Ts;
% M1=P1;
% [num1,den1]=tfdata(Gd1,'v');
% %% Design of MAC.............................................................................................................................................................
% t = 0:0.25:30;
% [h1,t1] =lsim(Gd1,[1 zeros(1,120)],t);
% figure(5);
% impulse(Gd1,t);
% %...............Toeplitz Matrix............................................
% b1=zeros(1,P1);
% b1(1,1)=h1(2);
% a1=h1(2:P1+1);
% H1=toeplitz(a1,b1);
% %.....................Hankek Matrix........................................
% c1=h1(3:P1+2);
% r1=[h1(P1+2:N1+1)' zeros(1,P1-1)];
% H1_=hankel(c1,r1);
% %........................................................................................
% H=H1;
% H_=H1_;
% 
% gamma=1/10;
% gain_DC1=(num1(1)+num1(2))/(den1(1)+den1(2));
% Q=eye(P1);
% R1=gamma*((gain_DC1)^2)*eye(M1);
% R=R1;
% Kmac=((H'*Q*H)+R)\(H'*Q);
% y(1) = 0;
% alpha=0.5;
% r=ones(150+P1+1,1);
% U=zeros(M1,1);
% U_=zeros(N1-1,1);
% ypast=zeros(150+P1,1);
% d(1)=0;
% ym=zeros(P1,1);
% y_1=0; 
% for t=1:150
%      y(t+1)=num1(1)*U(1)+num1(2)*U_(2)-den1(2)*y_1;
%      y_1=y(t+1);
%     yd(t+1)=y(t+1);
%     for i=1:P1
%         yd(t+i+1)=alpha*yd(t+i)+(1-alpha)*r(t+i+1);
%     end
%     d(t+1)=y(t+1)-ym(1);
%      E=(yd(t+2:t+P1+1)'-ypast(t+1:t+P1)-d(t+1)*ones(P1,1));
%      U=Kmac*E;
%      ym=H*U+H_*U_;
%      U_(2:N1-1)=U_(1:N1-2);
%      U_(1)=U(1);
%      u(t)=U(1);
%      ypast(t+2:t+P1+1)=H_*U_;
% end 
% figure(7);
%  plot(r,'b');
%  hold on
%  plot(y,'m');
%  %axis([0 150 0 3]);
%  figure(8);
%  plot(u);
% % %..................................................................
% % close all
% % clear all
% % clc
% % G_s = tf(1,[1 1]);
% % G_d = c2d(G_s,0.25,'impulse');
% % [num1,den1]=tfdata(G_d,'v');
% % figure
% % step(G_d)
% % sys_info = stepinfo(G_d);
% % ts = sys_info.SettlingTime;
% % tr=sys_info.RiseTime;
% % Ts = 0.25;
% % N = (ts)/Ts;
% % t = 0:0.25:30;
% % [h,t1] =lsim(G_d,[1 zeros(1,120)],t);
% % figure
% % impulse(G_d,t)
% % P = tr/Ts; M = P;
% % % toeplitz matrix
% % b = zeros(1,P); b(1,1)= h(2);
% % a = h(2:P+1);
% % H = toeplitz(a,b);
% % H(:,M)=H(:,M:P)*ones(P-M+1,1);
% % H=H(:,1:M);
% % % hankel matrix
% % c = [h(3:P+2)];
% % r = [(h(P+2:N+1))' zeros(1,P-1)];
% % H_minus = hankel(c,r);
% % %......................................................................
% % gamma = 1/2;
% % gain_DC=(num1(1)+num1(2))/(den1(1)+den1(2));
% % Q = eye(P,P);
% % R = gamma*gain_DC^2*eye(M,M);
% % K_mac = inv(H'*Q*H+R)*H'*Q;
% % y(1) = 0;
% % alpha = 0.5;
% % U_minus = zeros(N-1,1);
% % U=zeros(M,1);
% % yr = ones(150+P+1,1);
% % ypast = zeros(150+P,1);
% % d(1)= 0;
% % y_m=zeros(P,1);
% % y1 = 0;
% % y2 = 0;
% % for t = 1:150
% % y(t+1) =num1(1)*U(1) +num1(2)*U_minus(2) -den1(2)*y1;
% % % y2 = y1;
% % y1 = y(t+1);
% % yd(t+1) =y(t+1);
% % 
% % for i = 1:P
% % yd(t+i+1) = alpha*yd(t+i)+(1-alpha)*yr(t+i+1);% prog
% % %yd(t+i) = alpha*yd(t+i-1)+(1-alpha)*yr(t);% unprog
% % end
% % d(t+1) = y(t+1) - y_m(1);
% % E = (yd(t+2:t+P+1))' - ypast(t+1:t+P) - d(t+1)*ones(P,1);
% % U = K_mac * E;
% % y_m = H*U + H_minus*U_minus;
% % U_minus(2:N-1) = U_minus(1:N-2);
% % U_minus(1) = U(1);
% % u(t) = U(1);
% % ypast(t+2:t+P+1) = H_minus*U_minus;
% % end
% % %
% % figure
% % plot(yr)
% % hold on
% % plot(y,'r');
% 
% close all
% clear all
% clc
% G_s = tf(1,[1 1]);
% G_d = c2d(G_s,0.25,'impulse');
% [num1,den1]=tfdata(G_d,'v');
% G_s2 = tf(1,[1 1]);
% G_d2 = c2d(G_s2,0.25,'impulse');
% [num2,den2]=tfdata(G_d2,'v');
% figure
% step(G_d)
% sys_info = stepinfo(G_d);
% ts = sys_info.SettlingTime;
% tr=sys_info.RiseTime;
% Ts = 0.25;
% N = (ts)/Ts;
% t = 0:0.25:30;
% [h,t1] =lsim(G_d,[1 zeros(1,120)],t);
% figure
% impulse(G_d,t)
% P = tr/Ts; M = P;
% % toeplitz matrix
% b = zeros(1,P); b(1,1)= h(2);
% a = h(2:P+1);
% H = toeplitz(a,b);
% H(:,M)=H(:,M:P)*ones(P-M+1,1);
% H=H(:,1:M);
% H=[H H];
% % hankel matrix
% c = [h(3:P+2)];
% r = [(h(P+2:N+1))' zeros(1,P-1)];
% H_minus = hankel(c,r);
% H_minus=[H_minus H_minus];
% %......................................................................
% gamma = 1/7;
% gain_DC=(num1(1)+num1(2))/(den1(1)+den1(2));
% Q = eye(P,P);
% R = gamma*gain_DC^2*eye(M,M);
% R=[R zeros(M); zeros(M) R];
% K_mac = inv(H'*Q*H+R)*H'*Q;
% y(1) = 0;
% alpha = 0.5;
% U_minus = zeros(2*N-2,1);
% U=zeros(2*M,1);
% yr = ones(150+P+1,1);
% ypast = zeros(150+P,1);
% d(1)= 0;
% y_m=zeros(P,1);
% y1 = 0;
% y2 = 0;
% for t = 1:150
% y(t+1) =num1(1)*U(1) +num1(2)*U_minus(2) -den1(2)*y1+num2(1)*U(M+1)+num2(2)*U_minus(N+1);
% % y2 = y1;
% y1 = y(t+1);
% yd(t+1) =y(t+1);
% 
% for i = 1:P
% yd(t+i+1) = alpha*yd(t+i)+(1-alpha)*yr(t+i+1);% prog
% %yd(t+i) = alpha*yd(t+i-1)+(1-alpha)*yr(t);% unprog
% end
% d(t+1) = y(t+1) - y_m(1);
% E = (yd(t+2:t+P+1))' - ypast(t+1:t+P) - d(t+1)*ones(P,1);
% U = K_mac * E;
% y_m = H*U + H_minus*U_minus;
% U_minus(2:N-1) = U_minus(1:N-2);
% U_minus(1) = U(1);
% U_minus(N+1:2*N-2) = U_minus(N:2*N-3);
% U_minus(N) = U(M+1);
% u(t) = U(1);
% ypast(t+2:t+P+1) = H_minus*U_minus;
% end
% %
% figure
% plot(yr)
% hold on
% plot(y,'r');
% 
clear
clc
V=100;  dH=2e5; ro=1e3; Cp=1; roc=1e3; Cpc=1; ha=7e5;  K0=7.2e10; J=1e4; Ca0=1; T0=350; Tc0=350; 
U1=ones(101,1);
U2=ones(101,1);
x1(1)=0.0882; x2(1)=441.2;
y(1)=441.2;
for t=1:100
 x1(t+1)=x1(t)+0.25*(((U1(t)+100)/V)*(Ca0-x1(t))-K0*x1(t)*exp(-J/x2(t)));
 x2(t+1)=x2(t)+0.25*(((U1(t)+100)/V)*(T0-x2(t))-((-dH/(ro*Cp))*K0*x1(t)*exp(-J/x2(t)))+((roc*Cpc)/(ro*Cp*V))*(U2(t)+100)*(1-exp(-ha/((U2(t)+100)*roc*Cpc)))*(Tc0-x2(t)));
 y(t+1)=x2(t+1);
end
plot(y);
[n1,d1,n2,d2]=Inputsys(1);
Gs1 = tf(n1,d1);
Gd1 = c2d(Gs1,0.25,'impulse');
[num1,den1]=tfdata(Gd1,'v');
Gs2 = tf(n2,d2);
Gd2 = c2d(Gs2,0.25,'impulse');
[num2,den2]=tfdata(Gd2,'v');
U1(1)=0; U2(1)=0;
y_1=0; y_2=0;
for t=1:100
    yl(t+1) =num1(1)*U1(t+1) +num1(2)*U1(t) -den1(2)*y_1+num2(1)*U2(t+1)+num2(2)*U2(t)-den1(3)*y_2;
 y_2 = y_1;
y_1 = yl(t+1);
end
figure(2)
plot(yl)