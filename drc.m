%%
clc;
close all;
clear all;
%%
addpath('./CHINESE(MANDARIN)/');
% [nearEnds,fs]=wavread('Ch_f5');
% nearEnds=resample(nearEnds,1,2);
t=[0:1/8000:5]';
nearEnds=chirp(t,0,5,4000);
[h0,h1,g0,g1]=firpr2chfb(127,.45);
% fid=fopen('drc_test.pcm');
% nearEnds=fread(fid,'short');
% fclose(fid);
nearEnds=floor(32767*nearEnds);
N=64;
% L=4*N;
L=N;
% L=(N/2)*length(h0);
% wins=resample(h0',32,1)/32;
R=16;
Rs=60;
T=[30*ones(16,1);30*ones(17,1)];
Iters=R*floor(length(nearEnds)/R);
buffer=zeros(L,1);
putbuffer=zeros(L,1);
f=[0;1/N;1/R;1];
mm=[1;0.5^0.5;0;0];
% wins=fir2(L-1,f,mm)';
% wins=fir1(L-1,1/N,chebwin(L,Rs))';%2*pi/((N-1)*(2*B)+2*B),B=2*pi/(2*N)=pi/N
% fins=fir1(L-1,1/N,chebwin(L,Rs))';
wins=ones(N,1);
fins=wins;
gs=zeros(1+N/2,1);
tmpW=zeros(N,1);
bufferSB=zeros(N,L/R);
outputs=zeros(Iters+R-1,1);
outputsv=zeros(Iters+R-1,1);
out=zeros(R,1);
outv=zeros(R,1);
Sf=zeros(1+N/2,1);
alphaS=0.75;
Xsc=zeros(1+N/2,1);
dbg=[];
fid=fopen('dcr_real','rb');
Sr=fread(fid,'float');
fclose(fid);
sid=fopen('dcr_imag','rb');
Si=fread(sid,'float');
fclose(sid);
S=Sr+j*Si;
%[-pi/R,pi/R]
w=zeros(L,1);
W=zeros(L/R,1+N/2);
for m=1:1:Iters,
    buffer=[buffer(2:L);nearEnds(m,1)];
    out=w'*flipud(buffer);
    outputs(m)=out;
    if(mod(m-1,R)==0),
        tmpN=flipud(buffer.*wins);%xn[m]=x[n+m]*w[-m],xn[lN+r]
        tmpW(1,1)=sum(tmpN(1:N:end),1);
        for k=1:N-1,
            tmpW(1+k,1)=sum(tmpN(1+N-k:N:end),1);
        end
        tmpF=fft(tmpW,N);%BandPassed
        pF=diag(tmpF(1:1+N/2)*tmpF(1:1+N/2)');
        Sf=alphaS*Sf+(1-alphaS)*pF;
        En=log2(max(Sf,1));
        dbg=[dbg;En];
        Xsc=En;
        idx=find(En>T);
        Xsc(idx)=(En(idx)-T(idx))*0.5+T(idx);
        g=Xsc-En;
        alpha=ones(1+N/2,1).*447/32768;
        indx=find(gs>g);
        alpha(indx)=1751/32768;
        gs=(1-alpha).*gs+alpha.*g;
        gm=2.^(gs/2);
        W(1,:)=gm.';
        w=weight_transform(W,N,R,L/R);
        tmpF=[0;gm(2:N/2).*tmpF(2:N/2);0];
        tmpR=ifft([tmpF(1:1+N/2);conj(tmpF(N/2:-1:2))],N);
        %sum_{k}sum_{m}Vm[k]*f[n-m]*W(-(n-m)k,N),sum_{m=n:-1:n-L+1}f[n-m]sum_{k}(vm[n-m]=Vm[k]W(-(n-m)k,N))
        for k=0:1:R-1,%xn[m]=x[n+m]*w[-m],x[m]*w[n-m]*W(mk,N)=x[n-m]w[m]*W((n-m)k,N),W(nk,N)*x[n-m]w[m]W(-mk,N)
            gins=fins(1+k:R:end);%sum_{k}sum_{m=n:-1:n-L+1}Vm[k]*f[n-m]W(-(n-m)k,N)=sum_{m}vm[(n-m)]*f[n-m]
            acc=0;%y[n=rR]=sum_{m=n:-1:n-L+1}v[m=n:-R:end][(n-m)]*f[n-m]=f[0]*v(0)[0]+f[R]*v(1)[R]+f[2*R]*v()
            for p=0:1:L/R-1,
                acc=acc+bufferSB(1+mod(k+p*R,N),1+p)*gins(1+p,1);
            end
            outv(1+k,1)=N*R*acc;
        end
        outputsv(m:m+R-1)=outv;
    end
end

for m=1:R:Iters,
    buffer=[buffer(1+R:L);nearEnds(m:m+R-1)];%L-R
    tmpN=flipud(buffer.*wins);%xn[m]=x[m]*w[n-m],x'n[m]=x[n+m]*w[-m]
    tmpW(1)=sum(tmpN(1:N:end),1);%X'n[k]=sum_{m}(x'n[m]*W(mk,N))=sum_{r=0}^{N-1}(x(0)[(r)N]*W(k*r,N))
    for k=1:1:N-1,%x(0)[(r)N]=sum_{l}x'n[lN+r]
        tmpW(1+k,1)=sum(tmpN(N-k+1:N:end),1);%sum_{l}x'[lN+(r-n)]=x(0)[(r-n)N]
    end
    %tmpW=circshift(tmpW,m-1);
    %zn=tmpW(1:2:end)+j*tmpW(2:2:end);
    %tmp32=fft(zn,N/2);
    tmpF=fft(tmpW,N);%x(0)[(r-n)N]
    %dbg=[dbg;tmpF(1:1+N/2)];
    %index=1+((m-1)/R)*(1+N/2);
    %tmpF=S(index:index+N/2);
    pF=diag(tmpF(1:1+N/2)*tmpF(1:1+N/2)');
    Sf=alphaS*Sf+(1-alphaS)*pF;
    %En=10*log10(max(Sf,eps));%ref=1,q30
    En=log2(max(Sf,1));%10*log10(Sf)=-30dB,Sf=10^-3*2^30
    Xsc=En;
    idx=find(En>T);
    Xsc(idx)=(En(idx)-T(idx))*0.5+T(idx);
    g=Xsc-En;
    alpha=ones(1+N/2,1).*447/32768;
    indx=find(gs>g);
    alpha(indx)=1751/32768;
    gs=(1-alpha).*gs+alpha.*g;
    %dbg=[dbg;gs(2:N/2)];
    gm=2.^(gs/2);
    tmpF(1:1+N/2)=[0;gm(2:N/2).*tmpF(2:N/2);0];
    tmpR=real(ifft([tmpF(1:1+N/2);conj(tmpF(N/2:-1:2))],N));
    bufferSB=[circshift(tmpR,0),bufferSB(:,1:end-1)];%[new,...,old],Xm(e^jw)*W(n0k,N)<->um[(q-m-n0)N]
    %un[(q-n)N]=sum_{r}xn[(q-n)+rN]<->Xn(e^j(wk)),Xn(e^j(wk))*exp(-j*wkn0)<-
    %>un[(q-n+n0)N]
    %y[n]=Vk[n]*f[n]=sum_{k}(sum_{m}Vk[m]*f[n-m])W(-nk,N)=sum_{m}f[n-m]*vm[n]=sum_{m}f[-(m-n)]*vm[n]
    for k=0:1:R-1,%conv(Vm[k],g[m])=sum_{k}(sum_{m}(g[n-m]*Vm[k])W(-nk,N))=sum_{m}g[n-m]*vm[n]
        %out(1+k,1)=N*R*bufferSB(1+mod((m-1)+k,N),:)*fins(1+k:R:end);%sum_{k}Vm[k]*W(-nk,N)=vm[(n)N]
        gins=fins(1+k:R:end);
        acc=0;
        for p=0:1:L/R-1,%f[n-m]*(W(-(n-m)k,N)Vm[k]=vm[n-m]),sum_{m=n-L+1:n}f[n-m]*vm[n]
            acc=acc+bufferSB(1+mod(k+p*R,N),1+p)*gins(1+p,1);%vm[n]
        end
        out(1+k,1)=N*R*acc;
    end%f[-(m-n)]'*vm[n],f[-m]*[new,...,old]
    outputs(m:m+R-1)=out;%y[n]=sum_{m}g[n-m]vm[n]=sum_{m}g[-(m-n)]v[m]
    %um[q]=sum_{r}xm[q+rN],sum_{r}xm[(q-m+n0)+rN]==um[(q-m+n0)N]
    %     buffer=[buffer(end-N+R+1:end);nearEnds(m:m+R-1)];
    %     tmp=buffer.*wins;
    %     En=10*log10(max((tmp'*tmp)/N,eps));
    %     if(En>T),
    %         Xsc=(En-T)*0.5+T;
    %     else
    %         Xsc=En;
    %     end
    %     g=Xsc-En;
    %     gs=0.75*gs+0.25*g;
    %     gm=gs-T/2;
    %
    %     out=nearEnds(m:m+R-1)*10^(gs/20);
    %     putbuffer=[putbuffer(end-N+R+1:end);out];
    %     tmp=putbuffer.*wins;
    %     Et=10*log10(max(tmp'*tmp/N,eps));
    %     p=[p;out];
end
