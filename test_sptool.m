clc;
clear;
%%
b=[0.069735;0.388726;0.360530;0.388726;0.069735];
[sos,g]=tf2sos(b,1);
V1=2^ceil(log2(max(abs(b))));%E=log2(max(abs(b)))
bq=fi(b,1,8);%(-1)^S*M*2^(E),log2(0.5<M<=1)+E,-1<m<=0
dBgain=6;
N=256;
Q=2;
Fs=8000;
f0=1024;
A=10^(dBgain/40);
w0=2*pi*697/Fs;
alpha=sin(w0)/(2*Q);
d=[1+alpha*A;-2*cos(w0);1-alpha*A];
c=[1+alpha/A;-2*cos(w0);1-alpha/A];
K=1+alpha/A;
b=d./K;
a=c./K;

[f,w]=freqz(b,a,N);
h=remez(N,w/pi,abs(f));
