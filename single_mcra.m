clc;
clear all;
close all;
addpath('./aec_record');

%%
Fs=8000;
[xn, fs] = wavread('./CHINESE(MANDARIN)/Ch_m1');
[yn, fs] = wavread('./CHINESE(MANDARIN)/Ch_f2');
[zn, fs] = wavread('./CHINESE(MANDARIN)/Ch_m1');
cn = [xn; zeros(16384, 1); yn; zeros(16384, 1); zn]; % 16384/16k = 1024ms
[nn, fs] = wavread('./FS16K/white');
idx = min(length(cn), length(nn));
cn = cn(1:idx)/norm(cn(1:idx));%0dB
nn = nn(1:idx)/norm(nn(1:idx));%0dB
SNR = 10;
gain = 10^(-SNR/10);
xn = cn + nn .* gain^0.5;
duration = 10;%seconds
simpar.fs = 16e3;%sampling frequency
simpar.L = idx;
simpar.type = 'gusts'; %type of generated wind noise: constant -> 'const' or gusty -> 'gusts'
n_sim=generate_wind_noise(simpar);
% [n_sim,fs] = wavread('wind_strong.wav');
% n_sim=resample(n_sim,16000,fs);
n_sim = n_sim(1:idx)'/norm(n_sim(1:idx));
swr = 0;
gain = 10^(-swr/10);
xn = xn+gain*n_sim(1:length(xn))*gain^0.5;
xn = resample(xn, Fs, fs);
% [xn,fs] = wavread('wind_test.wav');
% xn = resample(xn,Fs,fs);

%%
N = 256;
R = N/2;
Iter = floor(length(xn)/R)*R;
alpha_d = 0.85;
alpha_s = 0.9;
alpha = 0.9844; % ML or |X(k,l-1)|.^2./lambda(k,l-1)
D = 32; % R/8*32/1000=0.512s
beta = 1.47;
Bmin = 1.66;
gamma0 = 4.6;
zeta0 = 1.67;
gamma1 = 3;
wins=hann(N,'period');
dcgain=R/sum(wins);
wins=sqrt(wins*dcgain);

xi_min = 10^(-15/10);%-15dB
Gmin = sqrt(10^-(30/10));
buffer = zeros(N, 1);
ola_buffer = zeros(N-R, 1);
Sf = zeros(1+N/2, 1);
S = zeros(1+N/2, 1);
Sf_tilt = zeros(1+N/2, 1);
S_tilt = zeros(1+N/2, 1);
gamma_prev = ones(1+N/2, 1);
GH1 = zeros(1+N/2, 1); % conditional gain based on H_{1}
q = ones(1+N/2, 1); % a priori speech absence probability, page 264
S_min = zeros(1+N/2, 1);
S_min_sw = zeros(1+N/2, 1);
S_min_tilt = zeros(1+N/2, 1);
S_min_sw_tilt = zeros(1+N/2, 1);
noise_b = zeros(1+N/2, 1);%bark-Domain
U = 1;
V = D;
cntr = 0;
P = zeros(1+N/2, Iter/R);
out = zeros(Iter, 1);
store_min = zeros(1+N/2, U);
store_min_tilt = zeros(1+N/2, U);
DBL = zeros(Iter/R, 1);
DBG = zeros(Iter/R, 1);
dbg = [];
dbl = [];
dbn=[];
[CB_FREQ_INDICES] = find_CB_FREQ_INDICES(Fs,N,16,R);
fid = fopen('test_out_ii','rb');
ST=fread(fid,'int32');
fclose(fid);
S1=ST(1:2:end)/2^15;
S2=ST(2:2:end)/2^15;
ST=S1+j*S2;
eps=2^-30;
Sxx=zeros(1+N/2,1);
Nxx=zeros(1+N/2,1);
Sii=zeros(1+N/2,1);
alphaMin=0.025;
alphaMax=0.975;
v=1;
Sk=zeros(1+N/2,1);
mup=logical([1;zeros(N/2,1)]);
mon=false(1+N/2,1);
mI=false(1+N/2,1);
Pm=zeros(1+N/2,floor(length(xn)/R));
Pu=zeros(1+N/2,floor(length(xn)/R));
Pg=zeros(1+N/2,floor(length(xn)/R));
bufferG=zeros(1+N/2,3);
Sobel=zeros(1+N/2,1);

for m = 1:R:Iter,
    buffer = [buffer(end-N+R+1:end); xn(m:1:m+R-1)];%OLA,N-R,R-Input
    spectra = fft(buffer.*wins, N);
    %     idx = 1+floor(m/R)*(R);
    %     stp = ST(idx:idx+N/2);
    %     spectra = [stp;conj(stp(end-1:-1:2))];
    temp=diag(spectra(1:1+N/2)*spectra(1:1+N/2)'); % periodogram, |Y(k,l)|.^2
    if m==1,
        lambda=beta*temp;
    end
    gamma=temp/beta./max(lambda,eps);
    bufferG=[bufferG(:,end-1:end),10*log10(max(gamma,1))];
    padgamma=[[zeros(1,3);bufferG;zeros(1,3)]];
    for k=2:N/2+2,
        Sobel(k-1,1)=(padgamma(k-1,3)/4+padgamma(k,3)/2+padgamma(k+1,3)/4)-(padgamma(k-1,1)/4+padgamma(k,1)/2+padgamma(k+1,1)/4);
    end
    xi = max(alpha .* (GH1.^2 .* gamma_prev) + (1-alpha) .* max(gamma-1, 0), xi_min);
    Sf(2:N/2) = temp(1:N/2-1)./4 + temp(2:N/2)./2 + temp(3:1+N/2)./4; % Smooth over frequency
    if m == 1,
        S = Sf;
    else
        S = alpha_s * S + (1-alpha_s) * Sf;
    end
    if m == 1,
        S_min = Sf;
        S_min_sw = Sf;
    else
        S_min = min(S_min, S);
        S_min_sw = min(S_min_sw, S);
    end
    zeta = S/Bmin./max(S_min,eps);
    gamma_min = temp(1:1+N/2)/Bmin./max(S_min,eps); % posterior SNR
    I = (gamma_min < gamma0 & zeta < zeta0); % speech absence
    I(1)=0;
    I(1+N/2)=0;
    index = zeros(1+N/2, 1);
    index(2:N/2, 1) = (I(1:N/2-1) + I(2:N/2) + I(3:1+N/2));
    idx = find(index ~= 0);
    Sf_tilt = S_tilt;
    Sf_tilt(idx) = (I(idx-1).*temp(idx-1)./4 + I(idx).*temp(idx)./2 + I(idx+1).*temp(idx+1)./4)./(I(idx-1)./4 + I(idx)./2 + I(idx+1)./4);
    if m == 1,
        S_tilt = Sf;
    else
        S_tilt = alpha_s.*S_tilt + (1-alpha_s).*Sf_tilt; % alp * S + (1-alp) * S
    end
    if m == 1,
        S_min_tilt = Sf;
        S_min_sw_tilt = Sf;
    else
        S_min_tilt = min(S_min_tilt, S_tilt);
        S_min_sw_tilt = min(S_min_sw_tilt, S_tilt);
    end
    gamma_min_tilt = temp(1:1+N/2)/Bmin./max(S_min_tilt,eps);
    dbg=[dbg;gamma_min_tilt(2:N/2)];
    zeta_tilt = S/Bmin./max(S_min_tilt,eps);
    q = (gamma_min_tilt <= 1 & zeta_tilt < zeta0).*(ones(1+N/2, 1)) + (gamma_min_tilt > 1 & gamma_min_tilt < gamma1 & zeta_tilt < zeta0).*((gamma1-gamma_min_tilt)/(gamma1-1));
    gamma_prev = gamma;
    xi_tilt = xi;
    nu = gamma .* (xi_tilt./(1+xi_tilt));
    %     b = (exp(nu)./(1+xi_tilt))<1;
    %     q = 0.1*b+(1-0.1)*q;
    loglikelihoods=nu-log(1+xi_tilt);
    DBL(1+floor(m/R), 1) = sum(loglikelihoods);
    % Q=1-q;
    % p = Q./(Q+q.*(1+xi_tilt).*exp(-nu)); % p(H_{1}|Y_{k})
    p=1-q;
    % exp_int = exp(0.5*expint(max(nu,eps)));
    % GH1 = (xi_tilt./(1+xi_tilt)) .* exp_int; % E[X_{k}|Y_{k}] =
    % E[X_{k}|Y_{k}, H_{1}] * p(H_{1}|Y_{k}) + E[X_{k}|Y_{k}, H_{0}] * p(H_{0}|Y_{k})
    GH1 = xi_tilt./(1+xi_tilt);
    alpha_d_tilt = alpha_d + (1-alpha_d)*(p);
    % dbl=[dbl;alpha_d_tilt(2:N/2)];
    % G = Gmin.^(1-p) .* GH1.^p; % olsa
    G = p.*GH1;
    for i = 1:length(CB_FREQ_INDICES),
        noise_b(CB_FREQ_INDICES{i}, 1) = ...
            ones(size(CB_FREQ_INDICES{i},2), 1) * mean(lambda(CB_FREQ_INDICES{i}));
    end
    lambda=alpha_d_tilt.*lambda+(1-alpha_d_tilt).*temp(1:1+N/2);
    % dbl=[dbl;lambda(2:N/2)];
    Tk=G.*spectra(1:1+N/2);
    Sxx=diag(Tk*Tk');
    mup=(10*log10(max(Sxx,eps))-10*log10(max(beta*lambda,eps)))>0;
    % mup=Sobel>0;
    for k=2:1:1+N/2,
        mon(k,1)=mup(k,1)&((k<9)||mon(k-1,1)||(mon(k,1)||mI(k,1)));
    end
    Sii=(mon>0).*Sxx+(mon==0).*min(Sii,Sxx);

    mI=or(mI,mon);
    for k=2:1:1+N/2,
        tt=(Sii(k,1)>(beta*lambda(k,1)));
        if(mI(k,1)>0),
            mI(k,1)=tt&(mI(k-1,1)||(k<9));
        end
    end

    Nxx=(mI>0).*Sii;
    SSC = 8000/N*([1:96]*Sxx(2:97)/max(sum(Sxx(2:97)),eps));
    Sk(2:N/2)=abs(Tk(1:N/2-1))/4+abs(Tk(2:N/2))/2+abs(Tk(3:1+N/2))/4;
    SPL=20*sum(abs(log10(max(Sk(2:91),1))-log10(max(Sk(1:90),1))));
    if (SSC<200), % ((a|b)|(c|d))'=((a'*b')*(c'*d'))
        % Nxx=diag(Tk*Tk');
        alphaLm=alphaMin;
    elseif (SSC>600),
        % Nxx=zeros(1+N/2,1);
        alphaLm=alphaMax;
    else
        dn=zeros(1+N/2,1);
        for k=0:N/2,
            for t=2:1:N/2,
                if t*k<N/2,
                    dn(1+k)=dn(1+k)+Sxx(1+k*t);
                end
            end
        end
        MaxMag=-Inf;
        nIndex=1;
        for k=3:N/2,
            if Sxx(k)>Sxx(k-1)&&Sxx(k)>Sxx(k+1),
                if MaxMag<Sxx(k),
                    MaxMag=Sxx(k);
                    nIndex=k;
                end
            end
        end
        nIndex=nIndex-1;
        alphaLm=(alphaMax*(SSC-200)+alphaMin*(600-SSC))/(600-200);
        locs=[nIndex:nIndex:min(1+N/2,3*nIndex)]+1;

        if length(locs)>2,
            Idx=[];
            for k=locs(1)+1:locs(2)-1,
                if Sxx(k-1)>Sxx(k) && Sxx(k)<Sxx(k+1),
                    Idx=[Idx;k];
                    break;
                end
            end

            for k=locs(2)+1:locs(3)-1,
                if Sxx(k-1)>Sxx(k) && Sxx(k)<Sxx(k+1),
                    Idx=[Idx;k];
                    break;
                end
            end
            if length(Idx)>1,
                a=log2(Sxx(Idx(1))/max(Sxx(Idx(2)),eps));
                b=log2((Idx(2)-1)/(Idx(1)-1));
                if ((a/b>2)||(a/b<0.5)),
                    v=v;
                else
                    v=a/b;
                end
                dbn=[dbn;v];
                BETA=(Idx(1)-1)^v*Sxx(Idx(1));
                % Nxx=BETA*max((0:1:N/2)',eps).^-v;
            end
        end
    end
    % DBL(1+floor(m/R))=SSC;
    DBG(1+floor(m/R))=SPL;
    % Spectra=Tk;
    Spectra=(mI==0).*Tk;
    tmpno1 = [Spectra;conj(Spectra(N/2:-1:2))];
    output = (ifft(tmpno1, N).*wins);
    tmpno2 = ola_buffer+output(1:N-R);
    out(m:m+R-1)=tmpno2(1:R);
    ola_buffer=[tmpno2(1+R:N-R);output(N-R+1:N)];
    Pg(2:1+N/2,ceil(m/R))=mI(2:1+N/2,1);
    Pu(2:1+N/2,ceil(m/R))=mup(2:1+N/2);
    Pm(2:1+N/2,ceil(m/R))=mon(2:1+N/2);
    if mod(cntr, V) == 0,
        store_min = [store_min(:,2:end), S_min_sw];%U
        S_min = min(store_min, [], 2);
        S_min_sw = S;
        store_min_tilt = [store_min(:,2:end), S_min_sw_tilt];
        S_min_tilt = min(store_min_tilt, [], 2);
        S_min_sw_tilt = S_tilt;
        %         S_min=S_min_sw;
        %         S_min_sw=S;
        %         S_min_tilt=S_min_sw_tilt;
        %         S_min_sw_tilt=S_tilt;
    end
    cntr = mod(cntr + 1, V);
end
%%
