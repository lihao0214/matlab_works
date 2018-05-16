clc;
clear all;
close all;
addpath('./aec_record');
addpath('./Common');
sim = 0;

%%
if ~sim,
    fid = fopen('echo004.raw', 'rb'); % near_end
    sid = fopen('ref004.raw', 'rb'); % far_end
    [nn, fs] = wavread('./FS16K/white');
    nn = resample(nn, 8000, fs);
    near_end = fread(fid, 'int16');
    far_end = fread(sid, 'int16');
    far_end = awgn(far_end,Inf,'measured');
    E_n = norm(near_end).^2;
    E_b = norm(nn(1:length(near_end))).^2;
    SNR = inf;
    SNR_c = 10*log10(E_n./E_b);
    SNR_d = SNR-SNR_c;
    nn = nn(1:length(near_end)).*10^(-SNR_d/20);
    near_end = near_end + nn;
    fclose(sid);
    fclose(fid);
else
    [near_end, fs] = wavread('.\CHINESE(MANDARIN)\Ch_f8');
    near_end = resample(near_end,8000,fs);
    ISM_RIR_bank(my_ISM_setup, 'ISM_RIRs.mat');
    AuData_near = ISM_AudioData('ISM_RIRs.mat', near_end);
    [far_end, fs] = wavread('.\CHINESE(MANDARIN)\Ch_m4');
    far_end = resample(far_end,8000,fs);
    ISM_RIR_bank(my_ISM_setup_i, 'ISM_RIRs_i.mat');
    AuData_far = ISM_AudioData('ISM_RIRs_i.mat', far_end);
    idx = min(length(AuData_near), length(AuData_far));
    near_end = [AuData_far(1:idx); AuData_near(1:idx); AuData_near(1:idx)+AuData_far(1:idx)];
    far_end = [far_end(1:idx); zeros(idx, 1); far_end(1:idx)];
end
%%
% near_end = 0.75*cos(2*pi*3000*[0:1:65535]'/16000);
% fid=fopen('echo00a.raw','rb');
% near_end=fread(fid,'short');
% fclose(fid);
% wins=hamming(65);
% b=firhalfband(64,wins);
% b0=b(1:2:end-1);%[h(0),...h(L-1)]
% b1=b(2:2:end);%
% ITER=2*floor(length(near_end)/2);
% xn=zeros(2,1);
% X=zeros(2,32);
% out=zeros(ITER/2,1);
% outH=zeros(ITER/2,1);
% for m=1:2:ITER,
%     xn=near_end(m:m+1);
%     X=[X(:,2:end),xn];
%     out(1+floor(m/2))=X(1,:)*fliplr(b1)'+X(2,:)*fliplr(b0)';
% end
% Y=zeros(2,32);
% for m=1:2:ITER,
%     xn=near_end(m:m+1).*[1,-1]';
%     Y=[Y(:,2:end),xn];
%     outH(1+floor(m/2))=Y(1,:)*fliplr(b1)'+Y(2,:)*fliplr(b0)';
% end
%
% Iter=length(out);
% Out=zeros(Iter*2,1);
% buffer=zeros(32,1);
% for m=1:1:Iter,
%     buffer=[buffer(2:end);out(m)];
%     Out(1+2*(m-1):2*m)=[fliplr(b0)*buffer;fliplr(b1)*buffer];
% end
%
% OutH=zeros(Iter*2,1);
% bufferH=zeros(32,1);
% for m=1:1:Iter,
%     bufferH=[bufferH(2:end);outH(m)];
%     OutH(1+2*(m-1):2*m)=[fliplr(b0)*bufferH;fliplr(b1)*bufferH].*[1,-1]';
% end
%%
% near_end=[near_end(1:23700);zeros(1024,1);near_end(23701:end-1024)];
% near_end=awgn(near_end,20,'measured');

Fs=8000;
alphaS=0.98;
N=256;
M=24;
MM=4;
R=32;
RR=128;
eps=1;
buffer_near=zeros(N,1);
buffer_far=zeros(N,1);
wins=hann(N,'period');
dcgain=RR/sum(wins);
wins=sqrt(wins*dcgain);
XFm=zeros(1+N/2,MM);
S_ss=zeros(1+N/2,MM);
S_vs=zeros(1+N/2,MM);
h=zeros(1+N/2,MM);
ola_buffer=zeros(N-RR,1);
iter = R*floor(length(near_end)/R);
Sdd=zeros(1+N/2,1);
GH1=ones(1+N/2,1);
alphaDD=0.6;
alphaQ=0.7;
Out=[];
gamma_prev=ones(1+N/2,1);
q=ones(1+N/2,1);
bandstart=5;
bandstop=129;
P=[];

% for m = 1:R:iter,
%     buffer_near = [buffer_near(1+R:end); near_end(m:m+R-1)];
%     yy = fft(buffer_near.*wins,N);
%     zr = ifft([yy(1:1+N/2);conj(yy(N/2:-1:2))],N);
%     zr = zr.*wins;
%     temp = OLA_Buf+zr;
%     out(m:m+R-1) = temp(1:R);
%     OLA_Buf = [temp(1+R:end);zeros(R,1)];
% end

mD=50+M;
DelaySxx=zeros(1+N/2,mD);
DelaySyy=zeros(1+N/2,mD);
S_xx=zeros(1+N/2,mD);
bxspectrum=zeros(1+N/2,mD);
byspectrum=zeros(1+N/2,1);
bcount=zeros(1+N/2,mD);
pcount=zeros(1+N/2,mD);
threshold=16;
pdelay=mD;
syy=zeros(1+N/2,1);
sxx=zeros(1+N/2,1);
D=[];
CC=[];
PP=[];
DD=[];
dbg=[];
dbl=[];
DBL=[];
fid = fopen('test_out_ii','rb');
ST=fread(fid,'int32');
fclose(fid);
sid = fopen('test_out_iii','rb');
SD=fread(sid,'int32');
fclose(sid);
Nr=ST(1:2:end)+j*ST(2:2:end);
Fr=SD(1:2:end)+j*SD(2:2:end);
Cntr=0;

for m = 1:R:iter,
    buffer_near = [buffer_near(end-N+R+1:end);near_end(m:m+R-1)];
    buffer_far = [buffer_far(end-N+R+1:end);far_end(m:m+R-1)];
    %     xx = fft(buffer_far.*wins,N);
    %     yy = fft(buffer_near.*wins,N);
    stamp=floor((m-1)/R);
    idx=[1+stamp*(1+N/2):(1+stamp)*(1+N/2)];
    xx=Fr(idx);
    yy=Nr(idx);
    Px = diag(xx(1:1+N/2)*xx(1:1+N/2)');
    Py = diag(yy(1:1+N/2)*yy(1:1+N/2)');
    %     if m == 1,
    %         parameters = initialise_parameters(Px./2^30,Fs,'imcra');
    %         parametersI=initialise_parameters(Py./2^30,Fs,'imcra');
    %     else
    %         parameters = noise_estimation(Px./2^30,'imcra',parameters);
    %         parametersI=noise_estimation(Py./2^30,'imcra',parametersI);
    %     end

    DelaySxx=[DelaySxx(:,2:end),xx(1:1+N/2)];
    DelaySyy=[DelaySyy(:,2:end),yy(1:1+N/2)];
    % xthreshold=median(abs(DelaySxx),2);
    % ythreshold=median(abs(DelaySyy),2);
    sxx=(1-1/64)*sxx+1/64*abs(xx(1:1+N/2));
    syy=(1-1/64)*syy+1/64*abs(yy(1:1+N/2));
    % xthreshold=mean(abs(DelaySxx(1:1+N/2,:)),2);
    % dbg=[dbg;sxx(1:1+N/2)];
    xthreshold=sxx;
    % ythreshold=mean(abs(DelaySyy(1:1+N/2,:)),2);
    ythreshold=syy;
    % dbg=[dbg;syy(1:1+N/2)];
    byspectrum=(abs(DelaySyy(1:1+N/2,mD))>=ythreshold);
    bxspectrum=[bxspectrum(:,2:end),(abs(DelaySxx(1:1+N/2,mD))>=xthreshold)];
    %     byspectrum=(parametersI.spp>0.75);
    %     bxspectrum=[bxspectrum(:,2:end),(parameters.spp>0.75)];

    for l=1:1:mD,
        bcount(:,l)=bitxor(bxspectrum(:,l),byspectrum);%s(q(x[n-D]),q(y[n]))
    end
    for l=1:1:mD,
        idx = find(bxspectrum(:,l))>0;
        pcount(idx,l)=(1-1/64)*pcount(idx,l)+(1/64)*bcount(idx,l);
    end
    Bcount=sum(pcount(bandstart:4:bandstop,:),1);%1^T
    % dbl=[dbl,Bcount];
    % dbl=[dbl;pcount(2:1+N/2,mD)];

    [Sxy,best_idx]=min(Bcount(M-3:mD),[],2);
    Idx=best_idx+M-4;
    if Sxy<=threshold,
        pdelay=Idx;
        threshold=Sxy;
    else
        threshold=min(max(threshold*(1+1/64),eps),8);
    end
    D=[D;pdelay];
    % P=[P;threshold];
    % XFm(:,1:MM)=DelaySxx(:,pdelay-M+4:4:pdelay);

    S_xx(:,1:mD-1)=S_xx(:,2:mD);
    S_xx(:,mD)=(1-1/256)*S_xx(:,mD)+(1/256)*diag(DelaySxx(:,mD)*DelaySxx(:,mD)');

    if mod(m-1,N/2)==96,
        XFm(:,1:MM-1)=XFm(:,2:MM);
        XFm(:,MM)=xx(1:1+N/2);
        S_ss(:,1:MM-1)=S_ss(:,2:MM);
        S_ss(:,MM)=alphaS*S_ss(:,MM)+(1-alphaS)*diag(XFm(:,MM)*XFm(:,MM)');
        % S_ss(:,1:MM)=S_xx(:,pdelay-M+4:4:pdelay);
        % PP=[PP;XFm(:,MM)];
        Vk=yy(1:1+N/2);
        Yk=yy(1:1+N/2);
        Ek=zeros(1+N/2,1);
        % DD=[DD;Vk];
        % PP=[PP;buffer_near(end-N/2+1:end)];
        for k=MM:-1:1,
            S_vs(:,k)=alphaS*S_vs(:,k)+(1-alphaS)*diag(Vk*XFm(:,k)');
            if(k==1&&Cntr==332),
                dbg=[dbg;diag(Vk*XFm(:,k)')];
            end
            h(:,k)=S_vs(:,k)./max(eps,S_ss(:,k));
            Vk=Vk-h(:,k).*XFm(:,k);
            Ek=Ek+h(:,k).*XFm(:,k);
        end
        PP=[PP;S_vs(2:N/2,MM-1)];
        CC=[CC;S_ss(2:N/2,MM-1)];
        DD=[DD;h(2:N/2,MM-1)];

        Thresh=abs(yy(1:1+N/2));
        absVf=max(Thresh,abs(Vk));
        absVf=Thresh./max(absVf,eps);
        Vk=Vk.*absVf;
        dbl=[dbl;Vk(2:N/2)];
        Cntr=Cntr+1;
    end
    %     gain=max(0,1-16*abs(Ek)./max(eps,abs(Vk)));
    %     Vk=Vk.*gain;
    %     Sdd=0.5*Sdd+0.5*diag(Ek*Ek');
    %     gamma=diag(Yk*Yk')./max(eps,Sdd);
    %     ksi=alphaDD*GH1.^2.*gamma_prev+(1-alphaDD)*max(gamma-1,0); %E[X(k,l)^2/lambdaD(k,l)|H1(k,l)]
    %     I=(gamma>3);
    %     q=alphaQ*q+(1-alphaQ)*I;
    %     gamma_prev=gamma;
    %     GH1=ksi./(1+ksi);%Sxx/(Sxx+Snn)
    %     Q=q./(1-q);%p(H1)/p(H0)
    %     nu=gamma.*ksi./(1+ksi);
    %     Delta=(1+ksi).^-1.*exp(nu);
    %     p=1./(1+Q.*Delta);
    % P=[P,parameters.spp];
    % Vk=(1-p).*GH1.*Yk;
    if mod(m-1,N/2)==96,
        temp=ifft([Vk;conj(Vk(N/2:-1:2))],N).*wins;
        out=temp(1:N-RR)+ola_buffer;
        Out=[Out;out(1:RR)];
        ola_buffer=[out(end-N+2*RR+1:end);temp(end-RR+1:end)];
    end
end
