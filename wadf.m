clc;
close all;
clear all;
addpath('./aec_record');
addpath('./Common');
sim = 1;
Fs = 8000;

%%
if ~sim,
    fid = fopen('echo004.raw', 'rb'); % near_end
    sid = fopen('ref004.raw', 'rb'); % far_end
    [nn, fs] = wavread('./FS16K/white');
    nn = resample(nn, 8000, fs);
    near_end = fread(fid, 'int16');
    far_end = fread(sid, 'int16');
    near_end = floor(near_end);
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
    [near_end, fs] = wavread('.\CHINESE(MANDARIN)\Ch_f1');
    near_end = resample(near_end, Fs, fs);
    ISM_RIR_bank(my_ISM_setup, 'ISM_RIRs.mat');
    AuData_near = ISM_AudioData('ISM_RIRs.mat', near_end);
    [far_end, fs] = wavread('.\CHINESE(MANDARIN)\Ch_m4');
    far_end = resample(far_end, Fs, fs);
    ISM_RIR_bank(my_ISM_setup_i, 'ISM_RIRs_i.mat');
    AuData_far = ISM_AudioData('ISM_RIRs_i.mat', far_end);
    idx = min(length(AuData_near), length(AuData_far));
    AuData_far = AuData_far(1:idx);
    AuData_near = AuData_near(1:idx);
    near_end = [AuData_far; AuData_near; AuData_near+AuData_far];
    near_end = awgn(near_end, 50, 'measured');
    near_end = round(near_end.*2^15);
    far_end = [far_end(1:idx); zeros(idx, 1); far_end(1:idx)];
    far_end = round(far_end.*2^15);
end

%%
M = 6;
N = 256;
R = N/2;
alpha = 0.98;
H = zeros(1+N/2, M);
buffer_f = zeros(N,1);
buffer_n = zeros(N,1);
XFm = zeros(1+N/2, M+1);
wins = tukeywin(N);
% wins = hanning(N);
threshold = 2e-6;
pn = zeros(1+N/2,1);
ITER = R*floor(length(near_end)/R);
mu = 0.6;
out = zeros(ITER+R-1, 1);
for m = 1:R:ITER,
    buffer_f = [buffer_f(1+R:end); far_end(m:m+R-1)];
    buffer_n = [buffer_n(1+R:end); near_end(m:m+R-1)];
    xx = fft(wins.*buffer_f, N);
    XFm = [XFm(:,2:end), xx(1:1+N/2)];
    pn = alpha.*pn+(1-alpha).*diag(XFm(:,M+1)*XFm(:,M+1)');
    tr = (XFm(:,2:M+1)-j.*XFm(:,1:M)).*H(:,1:M); % X(n)-j*X(n-1)
    ti = (conj(XFm(N/2:-1:2,2:M+1))-j.*conj(XFm(N/2:-1:2,1:M))).*conj(H(N/2:-1:2,1:M));
    nn = sum(tr, 2);
    vv = sum(ti, 2);
    temp = ifft([nn; vv], N);
    tt = real([j*temp(1+N/2:end); temp(1+N/2:end)]);
    en = buffer_n-tt;
    out(m:m+R-1) = en(1+R:end);
    ek = fft(en.*wins, N);
    ek2 = ek(1:1+N/2)./(eps+M*pn);
    absEf = max(abs(ek2), threshold);
    absEf = ones(N/2+1,1)*threshold./absEf;
    ek2 = ek2.*absEf;
    mek = conj(XFm(:,2:M+1)).*(ones(M,1)*ek2')'; % [Mx1:1x1+N/2]
    for k = 1:1:M,
        tmp = ifft([mek(:,k); conj(mek(N/2:-1:2,k))], N);
        tmp = [tmp(1:N/2); zeros(N/2,1)];
        xxx = fft(tmp, N);
        H(:,k) = H(:,k)+mu.*xxx(1:1+N/2);
    end
end
