close all
clear
clc
[n1,d1,n2,d2]=Inputsys(1);
Gs1 = tf(n1,d1);
Gd1 = c2d(Gs1,0.25,'impulse');
[num1,den1]=tfdata(Gd1,'v');
Gs2 = tf(n2,d2);
Gd2 = c2d(Gs2,0.25,'impulse');
[num2,den2]=tfdata(Gd2,'v');
figure(1);
step(Gd1);
hold on
step(Gd2,'r')
grid on
sys_info = stepinfo(Gd1);
ts = sys_info.SettlingTime;
tr=sys_info.RiseTime; 
Ts = 0.25;
N = ts/Ts;
t = 0:0.25:30;
[h1,t1] =lsim(Gd1,[1 zeros(1,120)],t);
[h2,t2] =lsim(Gd2,[1 zeros(1,120)],t);
figure(2)
impulse(Gd1,t)
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
gamma =1/10;
gain_DC=(num1(1)+num1(2)+num1(3))/(den1(1)+den1(2)+den1(3));
gain_DC2=(num2(1)+num2(2)+num2(3))/(den2(1)+den2(2)+den2(3));
Q = eye(P);
R1 = gamma*gain_DC^2*eye(M);
R2=((2.5)^2)*gamma*gain_DC2^2*eye(M);
R=[R1 zeros(M); zeros(M) R2];
Kmac = inv(H'*Q*H+R)*H'*Q;
y(1) = 0;
alpha = 0.5;
U_ = zeros(2*N-2,1);
U=zeros(2*M,1);
%r = ones(150+P+1,1);
 [r,t1]= gensig('sine',20,40,0.25);
ypast = zeros(150+P,1);
d(1)= 0;
ym=zeros(150+P,1);
y_1 = 0;
y_2 = 0;
x1(1)=0; x2(1)=0;

for t = 1:150
y(t+1) =num1(1)*U(1) +num1(2)*U_(2) -den1(2)*y_1+num2(1)*U(M+1)+num2(2)*U_(N+1)-den1(3)*y_2;
 y_2 = y_1;
y_1 = y(t+1);
yd(t+1) =y(t+1);

for i = 1:P
yd(t+i+1) = alpha*yd(t+i)+(1-alpha)*r(t+i+1); % Programmed
%yd(t+i) = alpha*yd(t+i-1)+(1-alpha)*r(t);% Unprogrammed
end
%...............................................................................
d(t+1) = y(t+1) - ym(1);
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
plot(r)
hold on
plot(y,'r');
grid on
