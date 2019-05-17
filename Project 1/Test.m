%close all
 clear
 clc
[n1,d1,n2,d2]=Inputsys(1);
Gs1 = tf(n1,d1);
Gd1 = c2d(Gs1,0.1,'impulse');
[num1,den1]=tfdata(Gd1,'v');
Gs2 = tf(n2,d2);
Gd2 = c2d(Gs2,0.1,'impulse');
[num2,den2]=tfdata(Gd2,'v');
% figure(1)
% step(Gd1)
% hold on
% step(Gd2,'r');
% grid on
% legend('u1','u2');
sys_info = stepinfo(Gd1);
ts = sys_info.SettlingTime;
tr=sys_info.RiseTime; 
sys_info = stepinfo(Gd2);
ts2 = sys_info.SettlingTime;
tr2=sys_info.RiseTime; 
%Ts=0.05;
Ts = 0.1;
%Ts=0.25;
N =floor( ts/Ts);
%N=30;
%N=3
t = 0:0.1:4;
[h1,t1] =lsim(Gd1,[1 zeros(1,40)],t);
[h2,t2] =lsim(Gd2,[1 zeros(1,40)],t);
% figure(2)
% impulse(Gd1,t);
% hold on
% impulse(Gd2,t,'r');
% grid on
% legend('u1','u2');
%P=3;
P =floor(tr/Ts); 
%P=10;
M = P; % In this method Control Horizon is equal to Predictive Horizon
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
gamma =1/2;
gain_DC=(num1(1)+num1(2)+num1(3))/(den1(1)+den1(2)+den1(3));
gain_DC2=(num2(1)+num2(2)+num2(3))/(den2(1)+den2(2)+den2(3));
Q = eye(P);
R1 =((1.4)^2)*gamma*gain_DC^2*eye(M);
R2=gamma*gain_DC2^2*eye(M);
R=[R1 zeros(M); zeros(M) R2];
Kmac = (H'*Q*H+R)\(H'*Q);
%alpha=0.3;
alpha=0.5;
%alpha=0.8;
x01=0.0882;
x02=441.2;

U1_ = zeros(N-1,length(t));
U2_ = zeros(N-1,length(t));
U_=[U1_ ; U2_];
d=zeros(1,length(t));
%y1=0; %linear
y1=441.2;
y=0;
Y_d=zeros(P,length(t));
Y_past=zeros(P,length(t));
Y_m=zeros(P,length(t));
D=zeros(P,length(t));
E_bar=zeros(P,length(t));
U1=zeros(M,length(t));
U2=zeros(M,length(t));
%..................step...........................
r =ones(length(t),1);
%....................pulse............................
%  [r,t1]= gensig('square',8,10,0.25);
% r=r;
%  for l = 1:length(r)
%  if (r(l)==0)
% r(l) = -1;
%  end
% end
%..................sine..................................
%[r,t1]= gensig('sine',8,10,0.25);
%......................step with various jump..............................
%  r1=[0 0 0 0 0 0 1 1 1 1 1 1 2 2 2 2 2 2 5 5 5 5 5 5 -1 -1 -1 -1 -1 -1 3 3 3 3 3 3]';
%  r=[r1; r1];
%......................................Step................................................
% r1=5*4.412*[1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0]';
% r=[r1];
%...........................................................................................................
U=[U1;U2];


for i=1:length(t)-1 
    
for j=1:P
  Y_d(j,i)=(alpha^j)*y+(1-(alpha)^j)*r(i); % Programmed
end 
    
Y_past(:,i)=H_*U_(:,i);
D(:,i)=d(i)*ones(P,1);
    
E_bar(:,i)=Y_d(:,i)-Y_past(:,i)-D(:,i);
U(:,i)=Kmac*E_bar(:,i);
U1(:,i)=U(1:M,i);
U2(:,i)=U(M+1:2*M,i);

U(:,i)=[U1(:,i);U2(:,i)];
    
Y_m(:,i)=H*U(:,i)+Y_past(:,i);


U1_(2:N-1,i+1) = U1_(1:N-2,i+1);
U1_(1,i+1)=U1(1,i);
U2_(2:N-1,i+1) = U2_(1:N-2,i+1);
U2_(1,i+1)=U2(1,i);


U_(:,i+1)=[U1_(:,i+1);U2_(:,i+1)];

u1=U1(1,i);
u2=U2(1,i);
sim('Model')
%d(i+1)=yl(end)-Y_m(1,i); %linear
d(i+1)=y(end)-Y_m(1,i);
%y1=[y1;yl(end)];  % linear
y1=[y1; y(end)+441.2];
x01=x1(end);
x02=x2(end);
%y=yl(end); % linear
y=y(end);    % nonlinear
end
figure(3);
plot(y1,'b');
hold on
plot(r+441.2,'r');
grid on
%axis([0 45 439 447]);
legend('y','r');
title('Response of the nonlinear system');
xlabel('sample');
%.............. Linear..............................
% plot(y1,'b');
% hold on
% % plot(r,'r');
% grid on
% title('Response of the linear system');
% xlabel('sample');
