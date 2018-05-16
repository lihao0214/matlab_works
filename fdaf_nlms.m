% function [ouput] = fdaf_nlms(far_end, near_end),
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
    [nn, fs] = wavread('./FS16K/pink');
    nn = resample(nn, 8000, fs);
    near_end = fread(fid, 'int16');
    far_end = fread(sid, 'int16');
    near_end = floor((near_end+0)./1);
    E_n = norm(near_end).^2;
    E_b = norm(nn(1:length(near_end))).^2;
    SNR = 55;
    SNR_c = 10*log10(E_n./E_b);
    SNR_d = SNR-SNR_c;
    nn = nn(1:length(near_end)).*10^(-SNR_d/20);
    near_end = near_end + nn;
    fclose(sid);
    fclose(fid);
else
    [near_end, fs] = wavread('.\CHINESE(MANDARIN)\Ch_f1');
    near_end = resample(near_end, 8000, fs);
    ISM_RIR_bank(my_ISM_setup, 'ISM_RIRs.mat');
    AuData_near = ISM_AudioData('ISM_RIRs.mat', near_end);
    [far_end, fs] = wavread('.\CHINESE(MANDARIN)\Ch_m4');
    far_end = resample(far_end, 8000, fs);
    ISM_RIR_bank(my_ISM_setup_i, 'ISM_RIRs_i.mat');
    AuData_far = ISM_AudioData('ISM_RIRs_i.mat', far_end);
    idx = min(length(AuData_near), length(AuData_far));
    AuData_far = AuData_far(1:idx);
    AuData_near = AuData_near(1:idx);
    near_end = [AuData_far; AuData_near; AuData_near+AuData_far];
    near_end = awgn(near_end, 50, 'measured');
    near_end = near_end.*32768;
    far_end = [far_end(1:idx); zeros(idx, 1); far_end(1:idx)];
    far_end = far_end.*32768;
end

%%
M = 12; % M*N/2/Fs
N = 128;
eta = 4;
buffer = zeros(N, 1);
iter = floor(length(far_end) / (N/2)) * (N/2);
XFm = zeros(1+N/2, M);
WFb = zeros(1+N/2, M);
Bmin = 1.66;
alpha = 1/8;
beta = 0.98;
pn = zeros(1+N/2, 1);
out = zeros(iter+N/2-1, 1);
ee = zeros(N, 1);
yy = zeros(N, 1);
dd = zeros(N, 1);
xn = zeros(N, 1);
mufb = 0.5;
threshold = 3.8147e-6;
d = zeros(1+N/2, iter./(N./2));
gamma = 1-1./4;
mbuf = zeros(N, 1);
Sed = ones(1+N/2, 1);
Sxd = ones(1+N/2, 1);
Sdd = ones(1+N/2, 1);
Syy = ones(1+N/2, 1);
See = ones(1+N/2, 1);
Sxx = ones(1+N/2, 1);
Sex = ones(1+N/2, 1);
wins = [0; sqrt(hanning(N-1))];
window = hanning(N);
ppsnr = ones(1+N/2, 1);
ppsnr_p = ones(1+N/2, 1);
Gain = ones(1+N/2, 1);
gain_p = ones(1+N/2, 1);
ovrd = 2.*(1+sqrt(linspace(0,1,1+N/2)));
vctr = zeros(floor(length(far_end)./(N./2)), 1);
fL = 300;
fH = 3300;
BandRange = 1+[round(fL/(8000/N)):round(fH/(8000/N))];
Thresh = 0.4;
Thresh_hd = 0.5;
scale = 5;
timer = 4;
Hpf = ones(1+N/2, 1);
q = ones(1+N/2, 1);
alpha_q = 0.3;
lambda = ones(1+N/2, 1);
far_end = filter([1,-1],[1,-0.95],far_end);
near_end = filter([1,-1],[1,-0.95],near_end);
ramp = 1.0003; % Upward ramp
alp = 1/8;
Snn = ones(1+N/2, 1);
Sym = (2^31-1)*ones(N/2+1,1);
step = 1/4; % Downward step size
yout = zeros(iter+N/2-1, 1);
eout = zeros(iter+N/2-1, 1);
alpha_s = 0.75;
zc = [];
divergeState = 0;
for m = 1:N/2:iter,
    buffer = [buffer(1+N/2:end); far_end(m:m+N/2-1)];
    xx = fft(buffer, N);
    pn = (1-alpha) .* pn + alpha .* diag(xx(1:1+N/2) * xx(1:1+N/2)');
    d(:, 1+floor(m/(N/2))) = pn;
    XFm = [XFm(:,2:M), xx(1:1+N/2)];
    Yk = XFm .* WFb;
    Y = sum(Yk, 2);
    y = ifft([Y; conj(Y(N/2:-1:2))], N);
    yout(m:m+N/2-1) = y(1+N/2:end);
    en = near_end(m:m+N/2-1)-y(1+N/2:end);
    eout(m:m+N/2-1) = en;
    eN = [zeros(N/2, 1); en];
    Ek = fft(eN, N);
    Ek2 = Ek(1:1+N/2)./2.^floor(log2(1+M.*pn));
    absEf = max(abs(Ek2), threshold);
    absEf = ones(N/2+1,1) * threshold ./ absEf;
    Ek2 = Ek2 .* absEf;
    mEk = conj(XFm).*(ones(M, 1) * Ek2(1:1+N/2)')';
    Ekn = [mEk; conj(mEk(N/2:-1:2,:))];
    for k = 1:1:M,
        ekn = ifft([mEk(:, k); conj(mEk(N/2:-1:2, k))], N);
        Ekn(:,k) = [ekn(1:N/2,1); zeros(N/2, 1)];
    end
    PHI = fft(Ekn, N);
    % PHI = mEk;
    ee = [ee(1+N/2:end); en];
    yy = [yy(1+N/2:end); y(1+N/2:end)];
    dd = [dd(1+N/2:end); near_end(m:m+N/2-1)];
    xn = [xn(1+N/2:end); far_end(m:m+N/2-1)];
    ef = fft(ee.*wins, N);
    yf = fft(yy.*wins, N);
    df = fft(dd.*wins, N);
    xf = fft(xn.*wins, N);
    Yp = diag(df(1:1+N/2)*df(1:1+N/2)'); % Instantaneous power
    Snn = (1-alp)*Snn + alp*Yp; % Averaged power

    mm = min(Snn, Sym);
    Sym = ((1-step).*mm + step.*Sym).*ramp;

    Sed = gamma .* Sed + (1-gamma) .* real(diag(ef(1:1+N/2)*df(1:1+N/2)'));
    Syy = gamma .* Syy + (1-gamma) .* diag(yf(1:1+N/2)*yf(1:1+N/2)');
    Sdd = gamma .* Sdd + (1-gamma) .* diag(df(1:1+N/2)*df(1:1+N/2)');
    See = gamma .* See + (1-gamma) .* diag(ef(1:1+N/2)*ef(1:1+N/2)');
    Sxx = gamma .* Sxx + (1-gamma) .* diag(xf(1:1+N/2)*xf(1:1+N/2)');
    Sxd = gamma .* Sxd + (1-gamma) .* real(diag(xf(1:1+N/2)*df(1:1+N/2)'));
    Sex = gamma .* Sex + (1-gamma) .* diag(ef(1:1+N/2)*xf(1:1+N/2)');
    cohed = diag(Sed*Sed')./(See.*Sdd+1e-10);
    cohxd = diag(Sxd*Sxd')./(Sxx.*Sdd+1e-10);
    hnled = min(1-cohxd, cohed);

    psnr = Sdd./(eps+See); % See ~= near-end speech, gain = 0, if echo only, high SNR mode, gain ~= 1
    xi = beta .* Gain.^2 .* ppsnr + (1-beta) .* max(psnr-1, 0);
    ppsnr = psnr;
    Gain = xi./(1+xi);
    bf = df(1:1+N/2).*Gain - yf(1:1+N/2);
    lambda = gamma.*lambda+(1-gamma).*diag(bf*bf');
    psnr_p = See./(eps+lambda);
    xi_p = beta .* gain_p.^2 .* ppsnr_p + (1-beta) .* max(psnr_p-1, 0);
    gain_p = xi_p./(1+xi_p);
    nu = psnr_p.*(xi_p./(1+xi_p));
    ppsnr_p = psnr_p;
    WFb = WFb + ((mufb.*ones(1+N/2, 1))*ones(1,M)) .* PHI(1:1+N/2, :);
    gain = min(abs(Sed./(Sdd + eta.*Syy)), 1); % over-weight post-filter
    gain_2 = min(abs(Sed./(Sed + eta.*Syy)), 1);
    gain_all = gain.*gain_2;
    Fmix = gain_all.*ef(1:1+N/2);

    ekEn = sum(See);
    dkEn = sum(Sdd);
    if divergeState == 0
        if ekEn > dkEn
            ef = df;
            divergeState = 1;
        end
    else
        if ekEn*1.05 < dkEn,
            divergeState = 0;
        else
            ef = df;
        end
    end

    if ekEn > dkEn*19.95,
        WFb = zeros(N/2+1,M); % Block-based FD NLMS
    end

    T = mask(diag((Fmix*2.^-15)*(Fmix*2.^-15)'), N, 8000, 16);
    H = min(sqrt((2.^30.*T)./(eps+lambda)), 1);
    I = ones(1+N/2, 1).*(xi_p < 5);
    q = alpha_q.*q + (1-alpha_q).*I;

    p = (1-q)./(1-q+q.*(1+xi_p).*exp(-nu));
    H_ovrdsm = (p.*H).^2;
    Fmix = H_ovrdsm .* ef(1:1+N/2);
    zeta = Fmix(BandRange)' * Fmix(BandRange) / (1+ef(BandRange)'*ef(BandRange));
    if 0,
        if zeta > Thresh,
            zeta_T = 1;
        else
            zeta_T = zeta;
        end

        if zeta_T == 1,
            N_T = 1;
        else
            N_T = 2*round((1-zeta_T./Thresh)*scale)+1;
        end

        H_lambda = zeros(1+N/2, 1);
        I = [1:N_T];
        H_lambda(I) = 1./N_T;
        Hpf = filter(H_lambda, 1, H_ovrdsm);
    else
        if zeta > Thresh_hd,
            tgain = H_ovrdsm;
            timer = 4; % 32ms
        else
            if timer > 0,
                timer = timer - 1;
                tgain = H_ovrdsm;
            else
                tgain = 0;
            end
        end
        Hpf = alpha_s.*Hpf+(1-alpha_s).*tgain;
    end

    hpf = Hpf.*exp(-j.*pi.*[0:1:N/2]').*exp(-j.*pi.*[0:1:N/2]'./N);
    hn = real(ifft([hpf; conj(hpf(N/2:-1:2))]));
    gn = hn .* window;
    snn = sqrt(Bmin.*Sym);
    snn(1) = 0; % Reject LF noise
    Un = snn.*exp(-j*2*pi.*[0;rand(N/2,1)]);

    % Weight comfort noise by suppression
    Un = sqrt(1-Hpf.^2).*Un;
    Fmix = Hpf .* ef(1:1+N/2);
    Fmix = Fmix + Un;
    ff = real(ifft([Fmix; conj(Fmix(N/2:-1:2))], N));
    mixw = wins.*(1.*ff);
    mola = mbuf(end-N/2+1:end) + mixw(1:N/2);
    mbuf = mixw;
    % [out(m:m+N/2-1), zc] = filter(gn,1,en,zc);
    out(m:m+N/2-1) = mola;
    vctr(1+floor(m/(N/2)), 1) = sum(1-q);
end
