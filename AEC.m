clc;
clear all;
close all;
addpath('./aec_record');
addpath('./Common');
addpath('./mfcc');

fid=fopen('TimeStretch.pcm','rb');
sid=fopen('ref004.pcm','rb');
near_end=fread(fid,'short');
far_end=fread(sid,'short');
fclose(fid);
fclose(sid);

%%
AEC_FLAG=0;
Fs=8000;
%K=64;
%R=K/2;
K=512;
R=K/2;
%R=K;%critical downsampling
Rs=60;
%L1=8*K;
L1=K;
L2=96/(Fs/1000)/R;
L3=48;%1.536s
buffer=zeros(R,L1/R);
buffern=zeros(L1,1);
bufferf=zeros(L1,1);
Iters=R*floor(length(near_end)/R);
%wins=ones(L1,1)/L1;
%wins=fir1(L1-1,1/K,chebwin(L1,Rs))';
dcgain=sum(hanning(L1))/R;
wins=sqrt(hanning(L1)/dcgain);
fins=wins;
tmpnW=zeros(K,1);
tmpfW=zeros(K,1);
bufferSB=zeros(K,L1/R);
bufferSS=zeros(K,L1/R);
outputs=zeros(R,1);
output=zeros(Iters+R-1,1);
P=reshape(wins,R,L1/R);
Gain=ones(1+K/2,2);
preSNR=zeros(1+K/2,2);

%%
%BSS
[s1,fs]=wavread('.\CHINESE(MANDARIN)\Ch_f1');
s1=resample(s1,Fs,fs);
s1=s1/max(abs(s1));
ISM_RIR_bank(my_ISM_setup_ii,'ISM_RIRs_ii.mat');
AuData_s1=ISM_AudioData('ISM_RIRs_ii.mat',s1);
[s2,fs]=wavread('.\CHINESE(MANDARIN)\Ch_m4');
s2=resample(s2,Fs,fs);
s2=s2/max(abs(s2));
ISM_RIR_bank(my_ISM_setup_iii,'ISM_RIRs_iii.mat');
AuData_s2=ISM_AudioData('ISM_RIRs_iii.mat',s2);
[s3,fs]=wavread('.\CHINESE(MANDARIN)\Ch_f2');
s3=resample(s3,Fs,fs);
s3=s3/max(abs(s3));
ISM_RIR_bank(my_ISM_setup_iv,'ISM_RIRs_iv.mat');
AuData_s3=ISM_AudioData('ISM_RIRs_iv.mat',s3);
[s4,fs]=wavread('.\CHINESE(MANDARIN)\Ch_f3');
s4=resample(s4,Fs,fs);
s4=s4/max(abs(s4));
ISM_RIR_bank(my_ISM_setup_v,'ISM_RIRs_v.mat');
AuData_s4=ISM_AudioData('ISM_RIRs_v.mat',s4);
idx=min([length(AuData_s1(:,1)),length(AuData_s2(:,1)),length(AuData_s3(:,1)),length(AuData_s4(:,1))]);
a1=1*AuData_s1(1:idx,1)+1*AuData_s2(1:idx,1)+0*AuData_s3(1:idx,1)+0*AuData_s4(1:idx,1);
a2=1*AuData_s1(1:idx,2)+1*AuData_s2(1:idx,2)+0*AuData_s3(1:idx,2)+0*AuData_s4(1:idx,2);
fid1=fopen('test_stream','rb');%near_end
fid2=fopen('test_stream1','rb');%far_end
s=fread(fid1,'short');
f=fread(fid2,'short');
% s=[s,s]/32768/4;
% f=f/32768;
s=[s,f]/32768;
fclose(fid1);
fclose(fid2);

% a1=wavread('x2-1');
% a2=wavread('x2-2');
% a1=wavread('mixture1');
% a2=wavread('mixture2');
%[s,fs]=wavread('mixed_babble');
%s=resample(s,Fs,fs);
% s=[[a1;a1],[a2;a2]];
% f=[s2;s2];

s=wavread('Capture.wav');
% s=[s;s;s;s;s;s];
% Tw=30;
% Ts=10;
% alpha=0.97;
% M=20;
% C=12;
% K=64;
% L=22;
% Fs=8000;
% LF=300;
% HF=3700;
% [MFCCs,FBEs,frames]=mfcc(s(:,1),Fs,Tw,Ts,alpha,@hamming,[LF,HF],M,C+1,L);
% u=mean(MFCCs,2);
% MFCCs=MFCCs-repmat(u,1,size(MFCCs,2));
aorig=[1,0;0,1];%(N*2)(2*2)
mixedsig=s*aorig';
%d=[0.02;-0.02];
d=[0.0275;-0.0275];%respect to origin
theta=[-pi/64;pi/64];
buffer_ch1=zeros(L1,1);
buffer_ch2=zeros(L1,1);
buffer_ch3=zeros(L1,1);
buffer_mxd=zeros(2,L3,1+K/2);
pca_inputs=zeros(2,L3,1+K/2);
ica_inputs=zeros(2,L3,1+K/2);
ym=zeros(2,L3,1+K/2);
zf=zeros(2,L3,1+K/2);
mf=ones(1+K/2,L3);
pref=zeros(2,L3,1+K/2);
nref=zeros(2,L3,1+K/2);%noise-ref
Iters=R*floor(length(s(:,1))/R);
dxh=zeros(2,2,1+K/2);
wxh=zeros(2,2,1+K/2);
pxh=zeros(2,2,1+K/2);
axh=zeros(2,2,1+K/2);
mxh=zeros(2,2,1+K/2);
nxh=zeros(2,2,1+K/2);
cxh=zeros(2,2,1+K/2);
fk=(0:1:K/2)'*Fs/K;
kc=2*pi*fk/346;%2*pi*f/c
for k=2:1:K/2,
    mxh(:,:,k)=exp(j*kc(k,1)*d*sin(theta'));
    wxh(:,:,k)=eye(2);%all
    pxh(:,:,k)=eye(2);%PCA
    dxh(:,:,k)=eye(2);%ICA
    axh(:,:,k)=eye(2);%beam
    cxh(:,:,k)=eye(2);
end
mu=0.01;
cntr=0;
Yk=zeros(1+K/2,2);
Pk=zeros(1+K/2,2);
Fk=zeros(1+K/2,1);
Output=zeros(2,Iters+R-1);
olv_buffer=zeros(L1-R,1);
ola_buffer=zeros(L1-R,1);
pm=zeros(1+K/2,1);
dbl=zeros(Iters/R,1);
rad=(-90:5:90)*pi/180;
DOA=zeros(2,1+K/2);
DOAs=zeros(2,1+K/2);%sorted DOA
res=size(rad,2);
ys2=zeros(2,1+K/2);
yc2=zeros(2,1+K/2);
ym2=zeros(2,1+K/2);
ym4=zeros(2,1+K/2);
%dbg=zeros(2,Iters/R+1);
%dbg=[];
dbm=zeros(2,Iters/R+1);
lambdaD=ones(1+K/2,2);
indices=zeros(2,1,1+K/2);
Vk=zeros(1+K/2,1);
Rnn=zeros(2,2,1+K/2);
T=zeros(2,2,1+K/2);
Atild=ones(2,1+K/2);%target-source
W=zeros(1+K/2,2);
for k=2:1:K/2,
    Rnn(1,:,k)=[1,sinc(kc(k)*d(1))];
    Rnn(2,:,k)=[sinc(kc(k)*d(2)),1];
    T(:,:,k)=pinv(Rnn(:,:,k)+eps*eye(2));
    W(k,:)=conj(Atild(:,k)'*T(:,:,k)/(Atild(:,k)'*T(:,:,k)*Atild(:,k)));
    %W(k,:)=exp(j*kc(k)*d*sin(pi/3))'/2;
end
egn=zeros(2,1+K/2);
hh=zeros(2,2*(1+K/2));
Sxx=ones(1+K/2,1);
Sxy=ones(1+K/2,1);
zc1=zeros(K,1);
zc2=zeros(K,1);
dbg_log=[];
fid=fopen('debug_log','wb');
G=@(W,z)...
    1./(0.1+abs(W'*z).^2);
F=@(W,z)...
    -1./(0.1+abs(W'*z).^2).^2;
%%
if 0,
    L=2048;
    K1=512;
    R1=128;
    state.cntr=0;
    state.wins=fir1(L-1,1/K1,'low',kaiser(L,2.5))';
    state.buffer1=zeros(L,1);
    state.buffer2=zeros(L,1);
    state.bufferSys1=zeros(K1,L/R1);
    state.bufferSys2=zeros(K1,L/R1);
    state.H=zeros(2,2,1+K1/2);
    for k=1:1:1+K1/2,
        state.H(:,:,k)=eye(2);
    end
    Iters=R1*floor(length(s(:,1))/R1);
    OutPuts=zeros(Iters,2);
    error=zeros(Iters/R1,1);
    for m=1:R1:Iters,
        [error(1+floor(m/R1)),OutPuts(m:m+R1-1,:),state]=pca_ica(mixedsig(m:m+R1-1,1),mixedsig(m:m+R1-1,2),K1,R1,state);
    end
end

for m=1:R:Iters,
    buffer_ch1=[buffer_ch1(1+R:L1);mixedsig(m:m+R-1,1)]+0*randn(L1,1)*0.0001;
    buffer_ch2=[buffer_ch2(1+R:L1);mixedsig(m:m+R-1,2)]+0*randn(L1,1)*0.0001;
    Zk1=fft(buffer_ch1.*wins(L1:-1:1),K);
    Zk2=fft(buffer_ch2.*wins(L1:-1:1),K);
    if(AEC_FLAG),
        buffer_ch3=[buffer_ch3(1+R:L1);f(m:m+R-1)];%s=[s[n];x[n]],[1,w;0,1]
        Zk3=fft(buffer_ch3.*wins(L1:-1:1),K);%x[n]
    end
    %dbg=[dbg;Zk1(2:K/2)];

    if(mod(cntr,L3)==0),
        if 0,
            for k=2:1:K/2,
                w=dxh(:,:,k)*pxh(:,:,k);%unmixing matrix
                u=pinv(w+eps*eye(2));
                h1=[abs(u(1,1));abs(u(2,1))*exp(j*angle(u(2,1)/u(1,1))/(kc(k)*2/pi)/abs(d(1)-d(2)))];
                h2=[abs(u(1,2));abs(u(2,2))*exp(j*angle(u(2,2)/u(1,2))/(kc(k)*2/pi)/abs(d(1)-d(2)))];
                h1=h1/norm(h1);
                h2=h2/norm(h2);
                wxh(:,:,k)=w;%not/corrected
                hh(:,1+2*(k-1))=h1;
                hh(:,2+2*(k-1))=h2;
            end
            Rh=(conj(hh)*hh.')/K/2;
            [v,D]=eig(Rh);
            dc=real(diag(D));
            [dummy,order]=sort(-dc);
            C=v(:,order);
            for k=2:1:K/2,
                p1=hh(:,2*(k-1)+1);
                p2=hh(:,2*(k-1)+2);
                j1=abs(C(:,1)'*p1);
                j2=abs(C(:,1)'*p2);
                [dummy,order]=sort([-j1;-j2]);
                wxh(:,:,k)=wxh(order,:,k);
            end
        else
            %p=sum(egn,2);
            for k=2:1:K/2,
                p=egn(:,k);
                if((p(1)>0.999*sum(p,1))),
                    w=wxh(:,:,k);
                else
                    w=dxh(:,:,k)*pxh(:,:,k);%(ICA/PCA),[x1;x2]=A*[s1;s2],w*[x1;x2]=(w*A)*[s1;s2]=(P*D)[s1;s2]=[d2s2;d1s1];
                end
                nxh(:,:,k)=pxh(:,:,k);
                %p=w*ica_inputs(:,:,k);
                %u=(ica_inputs(:,:,k)*p'/L3)*inv((p*p'+eps*eye(2))/L3);
                %u=pinv(w+eps*eye(2));%[a1,a2]
                %directive-patterns
                %amv=exp(j*kc(k,1)*d*sin(rad));%d=[d1;d2]
                %B1=w(1,:)*amv;%w1*exp(j*kc(k)*d*sin(*))
                %B2=w(2,:)*amv;%w2*exp(j*kc(k)*d*sin(*))
                %[f1,DOA1]=min(abs(B1),[],2);%beam-1
                %[f2,DOA2]=min(abs(B2),[],2);%beam-2
                %DOA(:,k)=-pi/2+[DOA1-1;DOA2-1]*(pi/(res-1));
                %[DOAs(:,k),idx]=sort(DOA(:,k));

                %DOA1=angle(u(1,1)/u(2,1))/(kc(k)*(d(1)-d(2)));%h(2,1)/h(1,1)
                %DOA2=angle(u(1,2)/u(2,2))/(kc(k)*(d(1)-d(2)));%h(2,2)/h(1,2)
                %DOA1=angle(-w(2,2)/w(2,1))/(kc(k)*(d(1)-d(2)));
                %DOA2=angle(-w(1,2)/w(1,1))/(kc(k)*(d(1)-d(2)));
                t1=w(2,2)*conj(w(2,1));
                DOA1=atan2(imag(-t1),real(-t1));
                t2=w(1,2)*conj(w(1,1));
                DOA2=atan2(imag(-t2),real(-t2));
                DOA(:,k)=[DOA1;DOA2]/kc(k)/(d(1)-d(2));%[exp(j*kc(k)*(d1*sin(theta)))][s1;s2]
                [dummy,indices(:,:,k)]=sort(abs(DOA(:,k)));
                DOAs(:,k)=DOA(indices(:,:,k),k);
                %u=u(:,indices(:,:,k));
                %wxh(:,:,k)=pinv(u+eps*eye(2));
                A=mxh(:,1,k)*mxh(:,1,k)'/4+ones(2)/2+mxh(:,2,k)*mxh(:,2,k)'/4;
                j1=real(w(1,:)*A*w(1,:)');
                j2=real(w(2,:)*A*w(2,:)');
                %[dummy,indices(:,:,k)]=sort([-j1;-j2]);
                wxh(:,:,k)=w(indices(:,:,k),:);
                %wxh(:,:,k)=w;
                axh(1,:,k)=0.5*[exp(-j*sign(d(1))/2*DOAs(1,k)),exp(-j*sign(d(2))/2*DOAs(1,k))];
                axh(2,:,k)=0.5*[exp(-j*sign(d(1))/2*DOAs(2,k)),exp(-j*sign(d(2))/2*DOAs(2,k))];
                tmp=wxh(:,:,k)*ica_inputs(:,:,k);
                u=pinv(wxh(:,:,k)+eps*eye(2));
                p1=u*[tmp(1,:);zeros(1,L3)];
                p2=u*[zeros(1,L3);tmp(2,:)];
                ym(:,:,k)=[sum(abs(p1),1);sum(abs(p2),1)]./2;
            end
            b1=abs(DOAs(1,:))<1;
            b2=abs(DOAs(2,:))<1;
            I=find(b1&b2);
            flag=[1;zeros(K/2-1,1);1];%[DC;Empty-Set;Nyquist]
            mDOA=mean(DOAs(:,I),2);
            vDOA=std(DOAs,1,2);
            [dummy,index]=sort(abs(mDOA));%[signal,noise]

            for k=2:1:K/2,%confidence/test,GetPermutations()
                if((norm(DOAs(:,k)-mDOA(:,1))<1.5*norm(vDOA(:,1)))),%deviates from mean DOAs
                    u1=wxh(1,:,k)*exp(j*kc(k)*d*DOAs(1,k));%null@w(2,:)
                    u2=wxh(2,:,k)*exp(j*kc(k)*d*DOAs(2,k));%null@w(1,:)
                    z1=wxh(1,:,k)*exp(j*kc(k)*d*DOAs(2,k));
                    z2=wxh(2,:,k)*exp(j*kc(k)*d*DOAs(1,k));
                    SIR1=10*log10(max(u1*u1',eps))-10*log10(max(z1*z1',eps));%SIR1=w1(DOA1)-w1(DOA2)
                    SIR2=10*log10(max(u2*u2',eps))-10*log10(max(z2*z2',eps));%SIR2=w2(DOA2)-w2(DOA1)
                    if(SIR1>10&&SIR2>10),
                        flag(k,1)=1;
                    end
                end
            end
            flag(1:1+K/2,1)=1;
            Corr=zeros(1+K/2,1)-512;
            Indices=cell(1+K/2,1);%recorded permutations
            I=find(flag(1:1+K/2)~=1);%I=freq/numbers
            while(~isempty(I)),
                for k=K/2:-1:2,
                    if(flag(k,1)~=1),
                        vf=ym(:,:,k);
                        vgm2=flag(max(k-2,2),1)*ym(:,:,max(k-2,2));
                        vgm1=flag(max(k-1,2),1)*ym(:,:,max(k-1,2));
                        vgp1=flag(min(k+1,K/2),1)*ym(:,:,min(k+1,K/2));
                        vgp2=flag(min(k+2,K/2),1)*ym(:,:,min(k+2,K/2));
                        cor=cov([vf',vgm2',vgm1',vgp1',vgp2'],1);
                        p=diag(sqrt(cor));
                        corrp1=cor(1,3)/max(eps,p(1)*p(3))+cor(1,5)/max(eps,p(1)*p(5))+cor(1,7)/max(eps,p(1)*p(7))+cor(1,9)/max(eps,p(1)*p(9));
                        corrp1=corrp1+cor(2,4)/max(eps,p(2)*p(4))+cor(2,6)/max(eps,p(2)*p(6))+cor(2,8)/max(eps,p(2)*p(8))+cor(2,10)/max(eps,p(2)*p(10));
                        corrp2=cor(2,3)/max(eps,p(2)*p(3))+cor(2,5)/max(eps,p(2)*p(5))+cor(2,7)/max(eps,p(2)*p(7))+cor(2,9)/max(eps,p(2)*p(9));
                        corrp2=corrp2+cor(1,4)/max(eps,p(1)*p(4))+cor(1,6)/max(eps,p(1)*p(6))+cor(1,8)/max(eps,p(1)*p(8))+cor(1,10)/max(eps,p(1)*p(10));
                        if(corrp1>=corrp2),
                            Corr(k,1)=corrp1;
                            Indices{k,1}=[1;2];
                        else
                            Corr(k,1)=corrp2;
                            Indices{k,1}=[2;1];
                        end
                    end
                end
                [cm,idx]=max(Corr(I));%may-track (flag==1)
                flag(I(idx),1)=1;
                Corr(I(idx),1)=-512;
                wxh(:,:,I(idx))=wxh(Indices{I(idx),1},:,I(idx));%optimum/permutate
                I=find(flag(1:1+K/2)~=1);
            end
            %% Harmonic/Tests
            if 0,
                for k=K/2:-1:2,
                    if(flag(k,1)~=1),
                        vf=abs(wxh(:,:,k)*ica_inputs(:,:,k));
                        vg2f=flag(2*k,1)*abs(wxh(:,:,2*k)*ica_inputs(:,:,2*k));
                        vg3f=flag(3*k,1)*abs(wxh(:,:,3*k)*ica_inputs(:,:,3*k));
                        vg2fm1=flag(2*k-1,1)*abs(wxh(:,:,2*k-1)*ica_inputs(:,:,2*k-1));
                        vg2fp1=flag(2*k+1,1)*abs(wxh(:,:,2*k+1)*ica_inputs(:,:,2*k+1));
                        vg3fm1=flag(3*k-1,1)*abs(wxh(:,:,3*k-1)*ica_inputs(:,:,3*k-1));
                        vg3fp1=flag(3*k+1,1)*abs(wxh(:,:,3*k+1)*ica_inputs(:,:,3*k+1));
                        cor=cov([vf',vg2fm1',vg2f',vg2fp1',vg3fm1',vg3f',vg3fp1'],1);
                        p=sqrt(diag(cor));
                        corrp1=cor(1,3)/max(eps,p(1)*p(3))+cor(2,4)/max(eps,p(2)*p(4));
                        corrp1=corrp1+cor(1,5)/max(eps,p(1)*p(5))+cor(2,6)/max(eps,p(2)*p(6));
                        corrp1=corrp1+cor(1,7)/max(eps,p(1)*p(7))+cor(2,8)/max(eps,p(2)*p(8));
                        corrp1=corrp1+cor(1,9)/max(eps,p(1)*p(9))+cor(2,10)/max(eps,p(2)*p(10));
                        corrp1=corrp1+cor(1,11)/max(eps,p(1)*p(11))+cor(2,12)/max(eps,p(2)*p(12));
                        corrp1=corrp1+cor(1,13)/max(eps,p(1)*p(13))+cor(2,14)/max(eps,p(2)*p(14));

                        corrp2=cor(2,3)/max(eps,p(2)*p(3))+cor(1,4)/max(eps,p(1)*p(4));
                        corrp2=corrp2+cor(2,5)/max(eps,p(2)*p(5))+cor(1,6)/max(eps,p(1)*p(6));
                        corrp2=corrp2+cor(2,7)/max(eps,p(2)*p(7))+cor(1,8)/max(eps,p(1)*p(8));
                        corrp2=corrp2+cor(2,9)/max(eps,p(2)*p(9))+cor(1,10)/max(eps,p(1)*p(10));
                        corrp2=corrp2+cor(2,11)/max(eps,p(2)*p(11))+cor(1,12)/max(eps,p(1)*p(12));
                        corrp2=corrp2+cor(2,13)/max(eps,p(2)*p(13))+cor(1,14)/max(eps,p(1)*p(14));

                        if(corrp1>corrp2),
                            if(corrp1>6*0.1*2),%|Ha|*2channels
                                flag(k,1)=1;
                                wxh(:,:,k)=wxh([1;2],:,k);
                            end
                        else
                            if(corrp2>6*0.1*2),%|Ha|*2channels
                                flag(k,1)=1;
                                wxh(:,:,k)=wxh([2;1],:,k);
                            end
                        end
                    end
                end
            end
        end

        for k=2:1:K/2,
            zf(:,:,k)=pxh(:,:,k)*ica_inputs(:,:,k);
            u=pinv(wxh(:,:,k)+eps*eye(2));
            u=pxh(:,:,k)*u;
            mf(k,:)=1./(1+exp(40*(angle(u(:,1)'*zf(:,:,k))-0.15*pi)));
            ym(:,:,k)=wxh(:,:,k)*ica_inputs(:,:,k);%org/sources
            ym4(:,k)=mean(abs(ym(:,:,k)).^4,2);%2*L3,E[|ym(:,1:L3)|.^4]
            ym2(:,k)=mean(abs(ym(:,:,k)).^2,2);%E[|ym(:,1:L3)|.^2]
            yc2(:,k)=mean(ym(:,:,k).^2,2);%E[ym(:,1:L3).^2]
            ys2(:,k)=std(ym(:,:,k),1,2);%E[|ym(:,1:L3)-uym|.^2]
        end
        Kurtosis=sum((ym4(:,2:K/2)-2*ym2(:,2:K/2).^2-abs(yc2(:,2:K/2)).^2)./max(ys2(:,2:K/2).^4,eps),2)/(K/2-1);
        for k=2:1:K/2,
            u=pinv(wxh(:,:,k)+eps*eye(2));%PB
            p=ym(:,:,k);%kurtosis/order,p=P*ym(:,:,k),p=(P*wxh(:,:,k))[x1;x2]
            pref(:,:,k)=u(:,1)*(mf(k,:).*p(1,:));%projection-back,inv(P*wxh(:,:,k))[p(1);0]=[s11;s21]
            %pref(:,:,k)=buffer_mxd(:,:,k);%post-process
            nref(:,:,k)=u(:,2)*p(2,:);%[s12;s22]
            ica_inputs(:,:,k)=buffer_mxd(:,:,k);%zero-means
            %Rc=cov(ica_inputs(:,:,k).',1);
            Rc=ica_inputs(:,:,k)*ica_inputs(:,:,k)'/L3;
            %dbg_log=[dbg_log;[Rc(1);Rc(3);Rc(2);Rc(4)]];
            %[v,D]=eig(Rc+trace(Rc)*eps*eye(2));
            %dc=real(diag(D));
            %[dummy,order]=sort(-dc);
            %pxh(:,:,k)=diag(dc(order,1)+eps)^-0.5*v(:,order)';
            b=zeros(2);
            for n=1:1:2,
                w=cxh(:,n,k);
                w=w/norm(w);
                for i=1:1:16,
                    w=Rc*w;
                    w=w-b*(b'*w);
                    w=w/max(eps,norm(w));
                end
                b(:,n)=w;
            end
            p=Rc*b;%[A*v1,A*v2]=[d1*v1,d2*v2]
            dc=sum(abs(p),1)./max(sum(abs(b),1),eps);
            [dummy,order]=sort(-dc);
            Rw=diag(max(dc(:,order),eps).^-0.5)*b(:,order)';
            %Rw=b(:,order)';
            %Rw=diag(max(dc,eps).^-0.5)*b';
            pxh(:,:,k)=Rw;
            %dbg_log=[dbg_log;[Rw(1);Rw(3);Rw(2);Rw(4)]];
            %dbg_log=[dbg_log;[Rc(1);Rc(3);Rc(2);Rc(4)]];
            %rDw=diag(max(dc(:,order),eps).^0.5)*b(:,order);
            rDw=b(:,order)'*diag(max(dc(:,order),eps).^-0.5);
            pca_inputs(:,:,k)=pxh(:,:,k)*ica_inputs(:,:,k);
            egn(:,k)=dc(order);

            dxh(:,:,k)=eye(2);%Initial/matrix
            %z=pca_inputs(:,:,k)';
            z=pca_inputs(:,:,k);
            %Phat=z.'*z/L3;
            Phat=z*z.'/L3;%diagonals
            if 1,
                if ~AEC_FLAG,
                    W=dxh(:,:,k)';%column-wise
                    Wold=zeros(2);%E[z*(conj(w'z)*g(*))]-E[g(*)+|w'z|^2*g'(*)]w
                    dbg_log=[];
                    if 0,
                        for n=1:1:L3,
                            Wold=W;
                            y=W'*z;
                            p=z*(conj(y).*G(W,z)).'/L3;%[nic,L]*[nic,L]'
                            q=((abs(y).^2).*F(W,z))*ones(L3,1)/L3;
                            r=G(w,z)*ones(L3,1)/L3;
                            W=p-W*diag(q)-W*diag(r);
                            [v,Lam]=eig(W'*W);
                            W=W*(v*diag(max(diag(real(Lam)),eps).^-0.5)*v');
                            error=norm(abs(Wold'*W)-eye(2),'fro');
                        end
                    else
                        dbg_log=[];
                        for n=1:1:L3,
                            Wold=W;%E[z*((w'*z)'*g(|w'*z|^2))]=E[z*(y'.*abs(y').^2)]
                            %y=z*W;
                            y=W'*z;
                            PhatW=Phat*W;%Pw=P*[w1,w2]=[Pw1,Pw2]
                            D=diag(diag(W.'*PhatW));
                            %W=(z'*(y.*abs(y).^2))/L3-2*W-conj(PhatW)*D;
                            W=(z*(y'.*abs(y').^2))/L3-2*W-conj(PhatW)*D;
                            [v,Lam]=eig(W'*W);%(W'*W)^-0.5*W'*W*(W'*W)^-0.5
                            W=W*(v*diag(max(diag(real(Lam)),eps).^-0.5)*v');
                            %                             for r=1:1:16,
                            %                                 W=W/max(sqrt(norm(W'*W,1)),eps);
                            %                                 W=1.5*W-0.5*W*W'*W;
                            %                             end

                            %                             Rc=W'*W;
                            %                             v=zeros(2);
                            %                             for q=1:1:2,
                            %                                 w=dxh(:,q,k);
                            %                                 for r=1:1:16,
                            %                                     w=Rc*w;
                            %                                     w=w-v*(v'*w);
                            %                                     w=w/max(norm(w),eps);
                            %                                 end
                            %                                 v(:,q)=w;
                            %                             end
                            %                             Lam=sum(abs(Rc*v),1)./max(sum(abs(v),1),eps);
                            %                             W=W*(v*diag((max(Lam,eps)).^-0.5)*v');
                            error=norm(abs(Wold'*W)-eye(2),'fro');%W'*W=v*D*v',v*D^-0.5*v'*(W'*W)*v*D^-0.5*v'=I
                            %fprintf(fid,'%f\n',error);
                        end
                    end
                    dxh(:,:,k)=W';
                else
                    t=rDw*[1;0];
                    b=t'*Rw;
                    b=b/(b(1)+eps);
                    pxh(:,:,k)=[[1,0];b];
                end
            else
                if 1,
                    t=zeros(2);
                    h1=[1;1];
                    cr(:,1)=pxh(:,:,k)*h1/norm(h1)^2;
                    cr(:,2)=pxh(:,:,k)*h1/norm(h1)^2;
                    A=mxh(:,1,k)*mxh(:,1,k)'/4+ones(2)/2+mxh(:,2,k)*mxh(:,2,k)'/4;
                    [v,D]=eig(A+trace(A)*eps*eye(2));
                    for q=1:1:2,
                        wold=zeros(2,1);
                        w=cr(:,q);
                        %w=dxh(q,:,k)';
                        for n=1:1:L3,
                            wold=w;
                            %PhatW=Phat*w;
                            zt=w'*z;
                            %gt=(abs(zt).^2).^3;%g=(z=|w'z|^2)^3
                            %gtt=3*(abs(zt).^2).^2;
                            gt=1./(0.1*ones(1,L3)+(abs(zt).^2));
                            gtt=-1./((0.1*ones(1,L3)+(abs(zt).^2)).^2);%E[z*conj(w'z)*g(*)]-E[g(*)+|w'z|^2*g'(*)]w
                            w=(z.*conj([1;1]*zt))*gt'/L3-(gt*ones(L3,1)/L3)*w-((abs(zt).^2)*gtt'/L3)*w;
                            a=(v*D^0.5*v')*(w-cr(:,q));
                            a=1e-4*a/max(eps,norm(a));
                            w=(v*D^-0.5*v')*a+cr(:,q);
                            %w=z*(zt'.*abs(zt').^2)/L3-2*w-conj(PhatW)*(w.'*PhatW);
                            %w=w/max(eps,abs(w'*cr(:,q)));
                            %w=w-t*(t'*w);
                            %w=w/max(eps,norm(w));
                        end
                        t(:,q)=w/max(eps,norm(w));
                    end
                    dxh(:,:,k)=t';
                else
                    t=zeros(2);
                    for q=1:1:2,
                        wold=zeros(2,1);
                        w=dxh(q,:,k)';
                        for n=1:1:L3,
                            wold=w;
                            PhatW=Phat*w;
                            %zt=z*w;
                            zt=w'*z;
                            %w=(z'*(zt.*abs(zt).^2))/L3-2*w-conj(PhatW)*(w.'*PhatW);
                            w=z*(zt'.*abs(zt').^2)/L3-2*w-conj(PhatW)*(w.'*PhatW);
                            w=w-t*t'*w;
                            w=w/norm(w);
                        end
                        t(:,q)=w;
                    end
                    dxh(:,:,k)=t';
                end
            end
        end
    end

    for k=2:1:K/2,
        %         for iter=1:1:384/L3,
        %             u=dxh(:,:,k)*pca_inputs(:,:,k);
        %             %nlr=-2*tanh(real(u))-j*2*tanh(imag(u));
        %             %nlr=-tanh(100*abs(u)).*exp(j*angle(u));
        %             nlr=-sign(real(u))-j*sign(imag(u));
        %             error=eye(2)+nlr*u'/L3;
        %             dxh(:,:,k)=dxh(:,:,k)+mu*error*dxh(:,:,k);
        %             [v,Lam]=eig(dxh(:,:,k)*dxh(:,:,k)');
        %             dxh(:,:,k)=(v*(diag(diag(real(Lam)).^-0.5)+eps*eye(2))*v')*dxh(:,:,k);
        %         end
        %         pm(k,1)=norm(error,'fro');
        if~AEC_FLAG,
            x=[Zk1(k,1);Zk2(k,1)];
        else
            x=[Zk1(k,1);Zk3(k,1)];
        end
        buffer_mxd(:,:,k)=[buffer_mxd(:,2:L3,k),x];
        y=wxh(:,:,k)*x;
        u=pinv(wxh(:,:,k)+eps*eye(2));
        nn=u(:,2)*y(2);
        ss=u(:,1)*y(1);
        p=nxh(:,:,k)*u(:,1);
        q=nxh(:,:,k)*x;
        mk=1./(1+exp(40*(angle(p'*q)-0.25*pi)));
        %Yk(k,:)=[x(1),nn(1)];
        %Yk(k,:)=[ss(1),nn(1)];
        %Yk(k,:)=reshape(pref(:,1+cntr,k),1,2);
        %Yk(k,:)=[y(1),y(2)];%[s1;s2]*[s1;s2]'=I
        Yk(k,:)=[ss(1),ss(2)]*mk;
        Pk(k,:)=[x(2),nn(2)];
    end
    %Sxy=(1-1/64)*Sxy+1/64*diag(Yk(:,2)*Yk(:,1)');
    %Sxx=(1-1/64)*Sxx+1/64*diag(Yk(:,2)*Yk(:,2)');
    lambdaD(:,1)=diag(Yk(:,2)*Yk(:,2)');
    postSNR=diag(Yk(:,1)*Yk(:,1)')./max(eps,lambdaD(:,1));
    %ksi=(1-1/64)*preSNR(:,1).*(Gain(:,1).^2)+1/64*max(postSNR-1,0);
    Gain(:,1)=(max(1-2*lambdaD(:,1)./max(eps,diag(Yk(:,1)*Yk(:,1)')),0)).^0.5;
    %Gain(:,1)=ksi./(1+ksi);
    %Gain(:,1)=Sxy./max(eps,Sxx);
    Vk(:,1)=Gain(:,1).*Yk(:,1);
    preSNR(:,1)=postSNR;
    lambdaD(:,2)=diag(Pk(:,2)*Pk(:,2)');
    postSNR=diag(Pk(:,1)*Pk(:,1)')./max(eps,lambdaD(:,2));
    %ksi=(1-1/64)*preSNR(:,2).*(Gain(:,2).^2)+1/64*max(postSNR-1,0);
    Gain(:,2)=max(1-2*lambdaD(:,2)./max(eps,diag(Pk(:,1)*Pk(:,1)')),0).^0.5;
    %Gain(:,2)=ksi./(1+ksi);
    Vk(:,2)=Gain(:,2).*Pk(:,1);
    preSNR(:,2)=postSNR;
    for k=2:K/2,
        Fk(k,1)=[0.5,0.5]*[Vk(k,1);Vk(k,2)];
    end
    zc1(1:1+K/2)=[zeros(20,1);Yk(21:K/2,1);0];
    zc1(2+K/2:K)=conj(zc1(K/2:-1:2));
    zc2(1:1+K/2)=[zeros(20,1);Yk(21:K/2,2);0];
    zc2(2+K/2:K)=conj(zc2(K/2:-1:2));
    z1=real(ifft(zc1,K));
    %dbg=[dbg;z1];
    z2=real(ifft(zc2,K));
    z1=z1.*wins;
    z2=z2.*wins;
    temp=ola_buffer+z1(1:L1-R);
    ola_buffer=[temp(1+R:end);z1(end-R+1:end)];
    Output(1,m:m+R-1)=temp(1:R)';
    temp=olv_buffer+z2(1:L1-R);
    olv_buffer=[temp(1+R:end);z2(end-R+1:end)];
    Output(2,m:m+R-1)=temp(1:R)';
    cntr=mod(cntr+1,L3);
end
fclose(fid);
rad=(-180:5:180)*pi/180;
doas=zeros(2,1+K/2);
bp1=zeros(1+K/2,length(rad));
bp2=zeros(1+K/2,length(rad));
for k=2:1:K/2,
    amv=exp(j*kc(k)*d*sin(rad));
    u=pinv(wxh(:,:,k)+eps*eye(2));
    doa1=angle(u(1,1)/u(2,1))/(kc(k)*(d(1)-d(2)));
    doa2=angle(u(1,2)/u(2,2))/(kc(k)*(d(1)-d(2)));
    doas(:,k)=[doa1;doa2];
    wxh(1,:,k)=wxh(1,:,k)/max(eps,norm(wxh(1,:,k)));
    wxh(2,:,k)=wxh(2,:,k)/max(eps,norm(wxh(2,:,k)));
    bp1(k,:)=wxh(1,:,k)*amv;
    bp2(k,:)=wxh(2,:,k)*amv;
end

subplot(1,2,1);contourf(rad,fk,abs(bp1),10);colorbar;xlabel('rad');ylabel('freq(Hz)');
subplot(1,2,2);contourf(rad,fk,abs(bp2),10);colorbar;xlabel('rad');ylabel('freq(Hz)');

if 0,%maximum-decimated/polyphase
    for m=1+R:R:Iters,
        buffer=[near_end(m:-1:m-R+1),buffer(:,1:L1/R-1)];%[new,...,old]
        for k=0:1:K-1,%polyphase_branches
            tmpnW(1+k,1)=buffer(1+k,:)*P(1+k,:)';%h[0,...,L1/R-1]
        end
        tmpnF=ifft(tmpnW,K);%sum_{r}z^(-r)sum_{l}hk[lN+r](z^N)^-l
        %hk[n]=h[n]*W(-nk,N),sum_{r}z^-r*W(-rk,N)*(H_{r}(z^N)=sum_{l}(h[lN+
        %r]*(z^N)^-l))
        tmpnR=real(fft([tmpnF(1:1+K/2);conj(tmpnF(K/2:-1:2))],K));%sum_{r}sum_{l}h[lN+r]z^(-lN-r),z^{-r}
        bufferSS=[tmpnR,bufferSS(:,1:L1/R-1)];%[new,...,old]
        for k=0:1:K-1,
            outputs(1+k,1)=K*R*bufferSS(K-k,:)*P(1+k,:)';%[0,1,...,R-1]
        end
        output(m-R:m-1)=outputs;%sum_{r}z^(-r)sum_{l}h[lN+r](z^N)^-l,[H0(z),...,HN-1(z)][1,z^-1,...,z^(-N+1)]'
    end%h[n]*W(-nk,N),sum_{r}z^-rsum_{l}h[lN+r]*W(-(lN+r)k,N)(z^N)^-l,sum_{r}z^-rW(-rk,N)sum_{l}h[lN+r](z^N)^-l
end

for m=1:R:Iters,%[1,...,1;1,W(-1*1,N),W(-2*1,N),...,W(-(N-1)*1,N);]
    buffern=[buffern(end-L1+R+1:end);near_end(m:m+R-1)];%nearEnds,[old,...,new]
    bufferf=[bufferf(end-L1+R+1:end);far_end(m:m+R-1)];%farEnds,[old,...,new]
    tmpn=flipud(buffern.*wins(L1:-1:1));%xn[m]=x[n+m]*w[-m]<-(I)DFS->Xn[k],un[0]=xn[0]+xn[0-N]+xn[0-2*N]+...
    tmpf=flipud(bufferf.*wins(L1:-1:1));%xn[m-n]:sum_{l}xn[lN+(r-n)],un[(r)N]=sum_{l}xn[lN+r]
    tmpnW(1,1)=sum(tmpn(1:K:end),1);%un[r=[0:1:N-1]]=sum_{l}xn[lN+r]
    tmpfW(2,1)=sum(tmpf(1:K:end),1);%un[(r)N]

    for k=1:1:K-1,
        tmpnW(1+k,1)=sum(tmpn(1+K-k:K:end),1);
        tmpfW(1+k,1)=sum(tmpf(1+K-k:K:end),1);
    end
    tmpnF=fft(tmpnW,K);
    tmpfF=fft(tmpfW,K);

    tmpnR=real(ifft([tmpnF(1:1+K/2);conj(tmpnF(K/2:-1:2))],K));

    bufferSB=[tmpnR(1:K),bufferSB(:,1:L1/R-1)];

    for k=0:1:R-1,%sum_{k}sum_{m=n:-1:n-L+1}Vm[k]*f[n-m]*W((n-m)k,N)
        acc=0;%sum_{m=n:-1:n-L+1}vm[n-m]*f[n-m]=f[0]*v[0][0]+f[R]*v[1][R]
        gins=fins(1+k:R:end);%f[1]*v[0][1]+f[1+R]*v[1][1+R],n-m=1,1+R,1+2R,...
        for p=0:1:L1/R-1,
            acc=acc+gins(1+p)*bufferSB(1+mod(k+p*R,K),1+p);
        end
        outputs(1+k,1)=K*R*acc;
    end
    output(m:m+R-1)=outputs;
end
