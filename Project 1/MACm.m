clear all
% close all
clc
format long
s=tf('s');
u1=11.9;
u2=2130;

sim('nonlinear1')
w_r1=w_r(end);
w_g1=w_g(end);
S1=S(end);
Teta1=Teta(end);
T_g1=T_g(end);

clear w_r w_g S Teta T_g
Ds=8e+4;
Jr=90000;
Ng=24.6;
Ks=8e+6;
Jg=10;
Tou_teta=0.15;
Tou_T=0.1;
Rou=1.2;
R=14.5;
v=14;
syms w_r w_g S Teta T_g
Landa=(v/R)/w_r;
Landa_t_inv=(1/(inv(Landa)+0.12*Teta))-(0.035/((1.5*Teta)^3+1));
C_p=0.22*((116*Landa_t_inv)-0.6*Teta-5)*exp(-12.5*Landa_t_inv);
P_r=0.5*Rou*pi*(R^2)*(v^3)*C_p;
f1=(P_r/(Jr*w_r))-((Ds*w_r)/Jr)+((Ds*w_g)/(Jr*Ng))-((Ks*S)/Jr);
f2=((Ds*w_r)/(Jg*Ng))-((Ds*w_g)/(Jg*Ng*Ng))+((Ks*S)/(Jg*Ng))-(T_g/Jg);
f3=w_r-(w_g/Ng);
f4=-Teta/Tou_teta;
f5=-T_g/Tou_T;

f11=diff(f1,w_r);f12=diff(f1,w_g);f13=diff(f1,S);f14=diff(f1,Teta);f15=diff(f1,T_g);

f21=diff(f2,w_r);f22=diff(f2,w_g);f23=diff(f2,S);f24=diff(f2,Teta);f25=diff(f2,T_g);

f31=diff(f3,w_r);f32=diff(f3,w_g);f33=diff(f3,S);f34=diff(f3,Teta);f35=diff(f3,T_g);

f41=diff(f4,w_r);f42=diff(f4,w_g);f43=diff(f4,S);f44=diff(f4,Teta);f45=diff(f4,T_g);

f51=diff(f5,w_r);f52=diff(f5,w_g);f53=diff(f5,S);f54=diff(f5,Teta);f55=diff(f5,T_g);

F=[f11 f12 f13 f14 f15 ; f21 f22 f23 f24 f25 ; f31 f32 f33 f34 f35 ; f41 f42 f43 f44 f45 ; f51 f52 f53 f54 f55];



A=double(subs(F,[w_r w_g S Teta T_g],[w_r1 w_g1 S1 Teta1 T_g1]));
B=[0 0 ; 0 0 ; 0 0 ; 1/Tou_teta 0 ; 0 1/Tou_T];
C=[1 0 0 0 0];

G=C*inv((s*eye(5)-A))*B;

save('G.mat','G');
save('w_r1.mat','w_r1');
save('w_g1.mat','w_g1');
save('S1.mat','S1');
save('Teta1.mat','Teta1');
save('T_g1.mat','T_g1');
save('u1.mat','u1');
save('u2.mat','u2');
clear all
clc
load('G.mat')
load('w_r1.mat');
load('w_g1.mat');
load('S1.mat');
load('Teta1.mat');
load('T_g1.mat');
load('u1.mat');
load('u2.mat');


T=10;
sampletime=0.001;
Ts=0.1;
M1=8;M2=8;P=5;
G1_d=c2d(G(1),Ts,'impulse');
G2_d=c2d(G(2),Ts,'impulse');
G1_info = stepinfo(G1_d);
G2_info = stepinfo(G2_d);
ts1=G1_info.SettlingTime;
ts2=G2_info.SettlingTime;
N1=ts1/Ts;
N2=ts2/Ts;
t1 = 0:Ts:ceil(20*ts1);
t2 = 0:Ts:ceil(20*ts2);
t  = 0:Ts:T;
time=0:sampletime:T;
[h1,t1] =lsim(G1_d,[1 zeros(1,length(t1)-1)],t1);
[h2,t1] =lsim(G2_d,[1 zeros(1,length(t2)-1)],t2);

%..................................................

% Toeplitz Matrix

b1 = zeros(1,P); b1(1,1)= h1(2);
a1 = h1(2:P+1);
H1 = toeplitz(a1,b1);
H1=H1*eye(P,M1);

b2 = zeros(1,P); b2(1,1)= h2(2);
a2 = h2(2:P+1);
H2 = toeplitz(a2,b2);
H2=H2*eye(P,M2);

H=[H1 H2];  

%..................................................
    
% Hankel Matrix

c1 = h1(3:P+2);
r1 = [(h1(P+2:N1+1))' zeros(1,P-1)];
H1_minus = hankel(c1,r1);

c2 = h2(3:P+2);
r2 = [(h2(P+2:N2+1))' zeros(1,P-1)];
H2_minus = hankel(c2,r2);

H_minus=[H1_minus H2_minus];

%..................................................

gamma = 2;
alpha = 0.9;
gain1_DC = dcgain(G1_d);
gain2_DC = dcgain(G2_d);
Q = eye(P,P);
R1 = gamma*gain1_DC^2*eye(M1);
R2 = gamma*gain2_DC^2*eye(M2);
R=[R1 zeros(M2);zeros(M1) R2];

K_mac = inv(H'*Q*H+R)*H'*Q;

%..................................................

U1_minus = zeros(N1-1,length(t));
U2_minus = zeros(N2-1,length(t));
U_minus=[U1_minus ; U2_minus];
d=zeros(1,length(t));
y1=[];
y=0;
Y_d=zeros(P,length(t));
Y_past=zeros(P,length(t));
Y_m=zeros(P,length(t));
D=zeros(P,length(t));
E_bar=zeros(P,length(t));
U1=zeros(M1,length(t));
U2=zeros(M2,length(t));
y_r = ones(length(t),1);
U=[U1;U2];


for i=1:length(t)-1 
    
for j=1:P
  Y_d(j,i)=(alpha^j)*y+(1-(alpha^j))*y_r(i);
end 
    
Y_past(:,i)=H_minus*U_minus(:,i);
D(:,i)=d(i)*ones(P,1);
    
E_bar(:,i)=Y_d(:,i)-Y_past(:,i)-D(:,i);
U(:,i)=K_mac*E_bar(:,i);
U1(:,i)=U(1:M1,i);
U2(:,i)=U(M1+1:M1+M2,i);

U(:,i)=[U1(:,i);U2(:,i)];
    
Y_m(:,i)=H*U(:,i)+Y_past(:,i);


U1_minus(2:N1-1,i+1) = U1_minus(1:N1-2,i+1);
U1_minus(1,i+1)=U1(1,i);
U2_minus(2:N2-1,i+1) = U2_minus(1:N2-2,i+1);
U2_minus(1,i+1)=U2(1,i);

% U1_minus(:,i+1)=[zeros(1,N1-1);[eye(N1-2) zeros(N1-2,1)]]*U1_minus(:,i)+U1(1,i)*[1;zeros(N1-2,1)];
% U2_minus(:,i+1)=[zeros(1,N2-1);[eye(N2-2) zeros(N2-2,1)]]*U2_minus(:,i)+U2(1,i)*[1;zeros(N2-2,1)];
U_minus(:,i+1)=[U1_minus(:,i+1);U2_minus(:,i+1)];

u11=U1(1,i);
u22=U2(1,i);
sim('nonlinear2')
d(i+1)=y(end)-Y_m(1,i);
y1=[y1;y(1:end-1)];
w_r1=w_r(end);
w_g1=w_g(end);
S1=S(end);
Teta1=Teta(end);
T_g1=T_g(end);
y=y(end);
end


    
  
    


