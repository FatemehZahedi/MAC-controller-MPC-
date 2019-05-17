 close all
 clear
 clc
V=100;  dH=2e5; ro=1e3; Cp=1; roc=1e3; Cpc=1; ha=7e5;  K0=7.2e10; J=1e4; Ca0=1; T0=350; Tc0=350; 
[n1,d1,n2,d2]=Inputsys(1);
Gs1 = tf(n1,d1);
Gd1 = c2d(Gs1,0.25,'impulse');
[num1,den1]=tfdata(Gd1,'v');
Gs2 = tf(n2,d2);
Gd2 = c2d(Gs2,0.25,'impulse');
[num2,den2]=tfdata(Gd2,'v');
figure(1)
step(Gd1)
hold on
step(Gd2,'r');
sys_info = stepinfo(Gd1);
ts = sys_info.SettlingTime;
tr=sys_info.RiseTime; 
Ts = 0.25;
N = ts/Ts;
t = 0:0.25:30;
[h1,t1] =lsim(Gd1,[1 zeros(1,120)],t);
[h2,t2] =lsim(Gd2,[1 zeros(1,120)],t);
figure(2)
impulse(Gd1,t);
hold on
impulse(Gd2,t,'r');
P = tr/Ts; M = P; % In this method Control Horizon is equal to Predictive Horizon
% ......................................................Toeplitz Matrix..........................................................
b = zeros(1,P); 
b(1,1)= h1(2);
a = h1(2:P+1);
H1 = toeplitz(a,b);
H1(:,M)=H1(:,M:P)*ones(P-M+1,1);
H1=H1(:,1:M);
b2 = zeros(1,P); 
b2(1,1)= h2(2);
a2 = h2(2:P+1);
H2 = toeplitz(a2,b2);
H2(:,M)=H2(:,M:P)*ones(P-M+1,1);
H2=H2(:,1:M);
H=[H1 H2];
% ...........................................................Hankel Matrix............................................................
c = h1(3:P+2);
r = [(h1(P+2:N+1))' zeros(1,P-1)];
H_1 = hankel(c,r);
c2 = h2(3:P+2);
r2 = [(h2(P+2:N+1))' zeros(1,P-1)];
H_2 = hankel(c2,r2);
H_=[H_1 H_2];
%...................................................................................................................................................
gamma =1;
gain_DC=(num1(1)+num1(2)+num1(3))/(den1(1)+den1(2)+den1(3));
gain_DC2=(num2(1)+num2(2)+num2(3))/(den2(1)+den2(2)+den2(3));
Q = eye(P);
R1 = gamma*gain_DC^2*eye(M);
R2=((2.5)^2)*gamma*gain_DC2^2*eye(M);
R=[R1 zeros(M); zeros(M) R2];
Kmac = (H'*Q*H+R)\(H'*Q);

y(1) = 441.2;
alpha = 0.5;
U_ = zeros(2*N-2,1);
U=zeros(2*M,1);
r =44.12*ones(150+P+1,1);
ypast = zeros(150+P,1);
ym=zeros(P,1);
y_1 = 0;
y_2 = 0;
x1(1)=0.0882; x2(1)=441.2;

for t = 1:150
% ..................................Nonlinear System.........................................................
for i=1:P
 x1(t+i)=x1(t+i-1)+0.25*(((U(1)+100)/V)*(Ca0-x1(t+i-1))-K0*x1(t+i-1)*exp(-J/x2(t+i-1)));
 x2(t+i)=x2(t+i-1)+0.25*(((U(1)+100)/V)*(T0-x2(t+i-1))-((-dH/(ro*Cp))*K0*x1(t+i-1)*exp(-J/x2(t+i-1)))+((roc*Cpc)/(ro*Cp*V))*(U(M+1)+100)*(1-exp(-ha/((U(M+1)+100)*roc*Cpc)))*(Tc0-x2(t+i-1)));
 y(t+i)=x2(t+i);
end
yd(t+1) =y(end)-441.2;

for i = 1:P
yd(t+i+1) = alpha*yd(t+i)+(1-alpha)*r(t+i+1); % Programmed
%yd(t+i) = alpha*yd(t+i-1)+(1-alpha)*r(t);% Unprogrammed
end

d(t+1) = y(end)-441.2- ym(1);
E = (yd(t+2:t+P+1))' - ypast(t+1:t+P) - d(t+1)*ones(P,1);
U = Kmac * E;
ym = H*U + H_*U_;
%..........................................................................................................................
U_(2:N-1) = U_(1:N-2);
U_(1) = U(1);
U_(N+1:2*N-2) = U_(N:2*N-3);
U_(N) = U(M+1);
u(t) = U(1);
ypast(t+2:t+P+1) = H_*U_;
end
%
figure(3)
plot(r+441.2)
hold on
plot(y,'r');
% U1minus = zeros(N-1,length(t));
% U2minus = zeros(N-1,length(t));
% Uminus=[U1minus ; U2minus];
% d=zeros(1,length(t));
% y1=[];
% y=0;
% Yd=zeros(P,length(t));
% Ypast=zeros(P,length(t));
% Ym=zeros(P,length(t));
% D=zeros(P,length(t));
% Ebar=zeros(P,length(t));
% U1=zeros(M,length(t));
% U2=zeros(M,length(t));
% y_r = ones(length(t),1);
% U=[U1;U2];
% alpha = 0.5;
% 
% for i=1:length(t)-1 
% 
% 
% for j=1:P
%   Yd(j,i)=(alpha^j)*y+(1-(alpha^j))*y_r(i);
% end 
%     
% Ypast(:,i)=H_*Uminus(:,i);
% D(:,i)=d(i)*ones(P,1);
%     
% Ebar(:,i)=Yd(:,i)-Ypast(:,i)-D(:,i);
% U(:,i)=Kmac*Ebar(:,i);
% U1(:,i)=U(1:M,i);
% U2(:,i)=U(M+1:2*M,i);
% 
% U(:,i)=[U1(:,i);U2(:,i)];
%     
% Ym(:,i)=H*U(:,i)+Ypast(:,i);
% 
% 
% U1minus(2:N-1,i+1) = U1minus(1:N-2,i+1);
% U1minus(1,i+1)=U1(1,i);
% U2minus(2:N-1,i+1) = U2minus(1:N-2,i+1);
% U2minus(1,i+1)=U2(1,i);
% 
% Uminus(:,i+1)=[U1minus(:,i+1);U2minus(:,i+1)];
% 
% u1(i)=U1(1,i);
% u2(i)=U2(1,i);
% x1(1)=0.0882;
% x2(1)=441.2;
%  x1(i+1)=x1(i)+0.25*(((u1(i)+100)/V)*(Ca0-x1(i))-K0*x1(i)*exp(-J/x2(i)));
% x2(i+1)=x2(i)+0.25*(((u1(i)+100)/V)*(T0-x2(i))-((-dH/(ro*Cp))*K0*x1(i)*exp(-J/x2(i)))+((roc*Cpc)/(ro*Cp*V))*(U(M+1)+100)*(1-exp(-ha/((U(M+1)+100)*roc*Cpc)))*(Tc0-x2(i)));
% y=x2(i+1)-441.2;
% y1=[y1; x2(i+1)];
% d(i+1)=y-Ym(1,i);
% end
% figure(3)
% plot(y1)
