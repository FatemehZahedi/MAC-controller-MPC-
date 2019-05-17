close all
clear all
clc
A=[-1.327*10^-5 1.617;-0.4 -3.434];
B=[-5 0;14.8 0.2];
C=[0 1];
D=[0 0];
[num11,den11]=ss2tf(A,B,C,D,1);
[num12,den12]=ss2tf(A,B,C,D,2);
g11=tf(num11,den11);
g12=tf(num12,den12);
gz11=c2d(g11,0.25,'impulse');
gz12=c2d(g12,0.25,'impulse');
T=0:0.25:25;
[h1,t1]=impulse(gz11,T);
h1=h1';
[h2,t2]=impulse(gz12,T);
h2=h2';
P=9;
M=9;
N=60;
%toeplitz H1
fr=[h1(1) zeros(1,P-1)]; %first row
fc=h1(1:P);              %first colomn
H1=toeplitz(fc,fr);
if M~=P
    if M<P
    H1(:,M)=H1(:,M:P)*ones(P-M+1,1);
    H1=H1(:,1:M);
    else
        M=P;
    end
    
end
%toeplitz H2
fr=[h2(1) zeros(1,P-1)]; %first row
fc=h2(1:P);              %first colomn
H2=toeplitz(fc,fr);
if M~=P
    if M<P
    H2(:,M)=H2(:,M:P)*ones(P-M+1,1);
    H2=H2(:,1:M);
    else
        M=P;
    end
    
end

H=zeros(P,2*M);
H(:,1:M)=H1;
H(:,M+1:2*M)=H2;

%hankel Hminus1
lr=[h1(P+1:N) zeros(1,P-1)]; %last row
fc=h1(2:P+1);                %first colomn
Hminus1=hankel(fc,lr);

lr=[h2(P+1:N) zeros(1,P-1)]; %last row
fc=h2(2:P+1);                %first colomn
Hminus2=hankel(fc,lr);


Hminus=zeros(P,2*(N-1));
Hminus(:,1:N-1)=Hminus1;
Hminus(:,N:2*(N-1))=Hminus2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q=eye(P,P);
gamma1=0.5;gamma2=0.5;
R1=gamma1*(dcgain(gz11))^2*eye(M,M);
R2=gamma2*(dcgain(gz12))^2*eye(M,M);
R=zeros(2*M,2*M);
R(1:M,1:M)=R1;
R(M+1:2*M,M+1:2*M)=R2;
a=H'*Q*H+R;
b=H'*Q;
K_mac = a\b; %%%%%%%constant
t=0:0.01:1;
U1minus = zeros(N-1,length(t));
U2minus = zeros(N-1,length(t));
Uminus=[U1minus ; U2minus];
d=zeros(1,length(t));
y1=[];
y=0;
Yd=zeros(P,length(t));
Ypast=zeros(P,length(t));
Ym=zeros(P,length(t));
D=zeros(P,length(t));
Ebar=zeros(P,length(t));
U1=zeros(M,length(t));
U2=zeros(M,length(t));
y_r = ones(length(t),1);
U=[U1;U2];


 alfa = .5;
%  rs >>> output
%  dil and rsin >>>> input


for i=1:length(t)-1 

rs(1)=0;
for j=1:P
  Yd(j,i)=(alfa^j)*rs(i)+(1-(alfa^j))*y_r(i);
end 
    
Ypast(:,i)=Hminus*Uminus(:,i);
D(:,i)=d(i)*ones(P,1);
    
Ebar(:,i)=Yd(:,i)-Ypast(:,i)-D(:,i);
U(:,i)=K_mac*Ebar(:,i);
U1(:,i)=U(1:M,i);
U2(:,i)=U(M+1:2*M,i);

U(:,i)=[U1(:,i);U2(:,i)];
    
Ym(:,i)=H*U(:,i)+Ypast(:,i);


U1minus(2:N-1,i+1) = U1minus(1:N-2,i+1);
U1minus(1,i+1)=U1(1,i);
U2minus(2:N-1,i+1) = U2minus(1:N-2,i+1);
U2minus(1,i+1)=U2(1,i);

Uminus(:,i+1)=[U1minus(:,i+1);U2minus(:,i+1)];

dil(i)=U1(1,i);
rsin(i)=U2(1,i);
rx(1)=5;
rs(1)=0;
% dil(1)=0.2;
% rsin(1)=0.15;
ki=40;
ks=0.1;
miomax=0.31;
yxs=0.5;
mio(i)=(miomax*rs(i))/(ks+rs(i)+(rs(i)^2)/ki);
rx(i+1)=mio(i)*rx(i)-dil(i)*rx(i);
rs(i+1)=-(mio(i)/yxs)*rx(i)+dil(i)*(rsin(i)-rs(i));
d(i+1)=rs(i)-Ym(1,i);

end


