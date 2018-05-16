clc;
close all;
clear all;
addpath('./aec_record');
addpath('./Common');

%%
N = 256;
R = 64;
Fs = 8000;
hr = [zeros(32,1);(-0.875).^[0:1:63]'];
gr = (-0.75).^[0:1:63]';
[s1, fs] = wavread('.\CHINESE(MANDARIN)\Ch_f2');
[s2, fs] = wavread('.\CHINESE(MANDARIN)\Ch_m4');
s1 = [s1; zeros(16384,1); s2];
s1 = resample(s1, Fs, fs);
s2 = filter(hr,1,s1);
% [nn, fs] = wavread('.\CHINESE(MANDARIN)\Ch_f4');
[nn, fs] = wavread('.\FS16K\pink');
nn=resample(nn,Fs,fs);
% nn=[nn;nn;nn;nn];
% nn=awgn(nn,20,'measured');
u=nn(1:length(s1));
SNRc=10*log10(norm(s1)^2/norm(u)^2);
SNRt=55;
u=u*sqrt(10^(-(SNRt-SNRc)/10));
w=filter(gr,1,u);
s(:,1)=s1+u;
s(:,2)=s2+w;
% s=wavread('mixed_babble_i_8k');

%%
[s1, fs] = wavread('.\CHINESE(MANDARIN)\Ch_m2');
s1 = resample(s1, Fs, fs);
s1=s1/max(abs(s1));
ISM_RIR_bank(my_ISM_setup_ii, 'ISM_RIRs_ii.mat');
AuData_s1 = ISM_AudioData('ISM_RIRs_ii.mat', s1);
[s2, fs] = wavread('./FS16K/white');
s2 = resample(s2, Fs, fs);
s2=s2/max(abs(s2));
ISM_RIR_bank(my_ISM_setup_iii, 'ISM_RIRs_iii.mat');
AuData_s2 = ISM_AudioData('ISM_RIRs_iii.mat', s2);
idx = min(length(AuData_s1(:,1)), length(AuData_s2(:,1)));
AuData_s2 = AuData_s2(1:idx,:);
AuData_s1 = AuData_s1(1:idx,:);
a1=AuData_s1(:,1)+AuData_s2(:,1);
a2=AuData_s1(:,2)+AuData_s2(:,2);
s=[a1,a2];
s=[s;s];

%%
dcgain=R/sum(hanning(N),1);
wins = sqrt(dcgain*hanning(N));
phi_z2z1=zeros(1+N/2,1);
phi_z1z1=zeros(1+N/2,1);
Iter=R*floor(length(s(:,1))/R);
L=10;
M=10;
P=6;
buffer_near=zeros(N,1);
buffer_near1=zeros(N,1);
phiEst1=zeros(1+N/2,L);
phiEst2=zeros(1+N/2,L);
a=load('RTF.mat');
Atilt=[ones(1+N/2,1),ones(1+N/2,1)];%[1,A2(k,l)/A1(k,l)],[A1/A1,A2/A1],sqrt(A1^2+A2^2)/A1
% Atilt=a.Atilt;
%S[A1,A2][A1/A1,A2/A1]*(A1/sqrt(A1^2+A2^2))^2,FBF outputs, A1*S(jw)
B=[-Atilt(:,2).';ones(1,1+N/2)];
Vk=zeros(1+N/2,1);
Uk=zeros(1+N/2,1);
Yk=zeros(1+N/2,1);
Ek=zeros(1+N/2,1);
ola_buffer=zeros(N-R,1);
ovr_buffer=zeros(N-R,1);
output=zeros(Iter+R-1,1);
out=zeros(Iter+R-1,1);
vad=zeros(1+N/2,Iter/R);
nk=zeros(1+N/2,1);
S_ss=ones(1+N/2,P);
S_vs=zeros(1+N/2,P);
XFm=zeros(1+N/2,P);
S_cc=ones(1+N/2,P);
S_vc=zeros(1+N/2,P);
XFc=zeros(1+N/2,P);
GH1=ones(1+N/2,1);
ksi=zeros(1+N/2,1);
alphaS=0.98*ones(1+N/2,1);
alpha=0.9;
nD=0;
SY=zeros(1+N/2,1);
SU=zeros(1+N/2,1);
lambdaD=ones(1+N/2,1);
preSNR=ones(1+N/2,1);
beta=1.47;
Gmin=sqrt(10^(-30/10));
alphad=0.85;
alphas=0.92;
q=ones(1+N/2,1);
freq=[0:1:N/2]*(Fs/N);
kc=2*pi*(freq)/346;%wave number
Rnn=zeros(2,2,1+N/2);
rad=[-90:10:90]*pi/180;
bP=zeros(1+N/2,length(rad));
T=zeros(2,2,1+N/2);
e=0.01;
for k=1:1:1+N/2,
    Rnn(1,:,k)=[1,sinc(kc(k)*0.08)];
    Rnn(2,:,k)=[sinc(kc(k)*0.08),1];
    T(:,:,k)=inv(Rnn(:,:,k)+e*eye(2));
end

for m=1:R:Iter, %(1+N/2,2),MTF approximations
    buffer_near=[buffer_near(R+1:N);s(m:m+R-1,1)];
    Zk1=fft(buffer_near.*wins,N);
    buffer_near1=[buffer_near1(R+1:N);s(m:m+R-1,2)];
    Zk2=fft(buffer_near1.*wins,N);
    z=[Zk1(1:1+N/2,1),Zk2(1:1+N/2,1)];
    for k=1:1:1+N/2,
        % T=inv(Rnn(:,:,k)+e*eye(2));
        % W=Atilt(k,:)/(Atilt(k,:)*Atilt(k,:)');
        W=conj(Atilt(k,:))*T(:,:,k)/(conj(Atilt(k,:))*T(:,:,k)*Atilt(k,:).');
        amv=exp(j*kc(k)*[1;0.08]*sin(rad));
        bP(k,:)=W*amv;%beam-pattern
        Vk(k,1)=(W*z(k,:).')*exp(-j*2*pi*nD*(k-1)/N);%*F(jw)
        B(1,k)=-Atilt(k,2);%column space of B orthogonal to A
        H=B(:,k);
        Uk(k,1)=z(k,:)*H;
    end

    XFm=[XFm(:,2:end),Ek];
    S_ss(:,1:P-1)=S_ss(:,2:P);
    S_ss(:,P)=alphaS.*S_ss(:,P)+(1-alphaS).*diag(XFm(:,P)*XFm(:,P)');
    Thresh=abs(Vk);
    for k=P:-1:1,
        S_vs(:,k)=alphaS.*S_vs(:,k)+(1-alphaS).*diag(Vk*XFm(:,k)');
        h=S_vs(:,k)./max(eps,S_ss(:,k));
        Vk=Vk-XFm(:,k).*h;
    end
    absVf=max(abs(Vk),Thresh);
    absVf=Thresh./max(eps,absVf);
    Vk=Vk.*absVf;

    Ek=Uk;
    XFc=[XFc(:,2:end),Yk];
    Yk=Vk;
    S_cc(:,1:P-1)=S_cc(:,2:P);
    S_cc(:,P)=alphaS.*S_cc(:,P)+(1-alphaS).*diag(XFc(:,P)*XFc(:,P)');
    Thresh=abs(Ek);

    for k=P:-1:1,
        S_vc(:,k)=alphaS.*S_vc(:,k)+(1-alphaS).*diag(Ek*XFc(:,k)');
        hc=S_vc(:,k)./max(eps,S_cc(:,k));
        Ek=Ek-XFc(:,k).*hc;
    end
    absVf=max(abs(Ek),Thresh);
    absVf=Thresh./max(eps,absVf);
    Ek=Ek.*absVf;

    Pz=diag(Vk*Vk');
    Pu=diag(Ek*Ek');
    if m==1,
        parameters = initialise_parameters(Pz,Fs,'mcra');
        param = initialise_parameters(Pu,Fs,'mcra');
    else
        parameters = noise_estimation(Pz,'mcra',parameters);
        param = noise_estimation(Pu,'mcra',param);
    end
    v1=circshift(Vk,1);
    v2=circshift(Vk,-1);
    SY=alpha*SY+(1-alpha)*([diag(v1*v1'),diag(Vk*Vk'),diag(v2*v2')]*[0.25;0.5;0.25]);
    u1=circshift(Ek,1);
    u2=circshift(Ek,-1);
    SU=alpha*SU+(1-alpha)*([diag(u1*u1'),diag(Ek*Ek'),diag(u2*u2')]*[0.25;0.5;0.25]);
    r1=max(SY-parameters.noise_ps,0);
    r2=max(SU-param.noise_ps,e*parameters.noise_ps);
    tr=r1./r2;
    % p=parameters.pk;
    p1=SY./max(parameters.noise_ps,e);
    p2=SU./max(param.noise_ps,e);
    p3=diag(Vk*Vk')./max(parameters.noise_ps,e);
    c1=find((p1>1.67)&(p2<1.81));
    c2=find((p1>1.67)&(p2>=1.81)&(p3>4.6)&(tr>3));
    postSNR=diag(Vk*Vk')./max(lambdaD,e);
    ksi=alphas*GH1.^2.*preSNR+(1-alphas)*max(postSNR-1,0);
    preSNR=postSNR;
    nu=postSNR.*ksi./(1+ksi);
    b=(exp(nu)./(1+ksi))<1;
    % b(c1)=0;
    % b(c2)=0;
    q=0.1*b+(1-0.1)*q;
    p=(1-q)./(1-q+q.*(1+ksi).*exp(-nu));
    alphaD=alphad+(1-alphad)*p;
    lambdaD=alphaD.*lambdaD+(1-alphaD)*beta.*diag(Vk*Vk');
    expInt=exp(0.5*expint(max(nu,eps)));
    GH1=ksi./(1+ksi).*expInt;
    G=(GH1.^p).*(Gmin.^(1-p));
    Ok=G.*Vk;
    Idx=find(p>0.5);%replace via others VAD
    vad(:,1+floor(m/R))=p;
    %Idx=[];
    nk(Idx)=nk(Idx)+1;
    phi_z2z1(Idx)=phi_z2z1(Idx)+1/M*(Zk2(Idx).*conj(Zk1(Idx)));
    phi_z1z1(Idx)=phi_z1z1(Idx)+1/M*(Zk1(Idx).*conj(Zk1(Idx)));
    idx=find(nk==L);
    phiEst1(idx,1:L-1)=phiEst1(idx,2:L);
    phiEst1(idx,L)=phi_z1z1(idx);
    phiEst2(idx,1:L-1)=phiEst2(idx,2:L);
    phiEst2(idx,L)=phi_z2z1(idx);
    a1=sum(phiEst1(idx,1:L).*phiEst2(idx,1:L),2)/L;
    a2=sum(phiEst1(idx,1:L),2)/L;
    a3=sum(phiEst2(idx,:),2)/L;
    a4=sum(phiEst1(idx,:).^2,2)/L;
    a5=a2.^2;
    temp=(a1-a2.*a3);
    tmpe=(a4-a5);
    ndx=find(tmpe==0);
    tmpe(ndx)=eps;
    Atilt(idx,2)=temp./tmpe;
    phi_z2z1(idx)=0;
    phi_z1z1(idx)=0;
    nk(idx)=0;

    y=real(ifft([Vk;conj(Vk(N/2:-1:2))],N));
    y=y.*wins;
    temp=ovr_buffer+y(1:N-R);
    out(m:m+R-1)=temp(1:R);
    ovr_buffer=[temp(1+R:N-R);y(N-R+1:N)];

    y=real(ifft([Ok;conj(Ok(N/2:-1:2))],N));
    y=y.*wins;
    temp=ola_buffer+y(1:N-R);
    output(m:m+R-1)=temp(1:R);
    ola_buffer=[temp(1+R:N-R);y(N-R+1:N)];
end

savefile='RTF.mat';
save(savefile,'Atilt');
