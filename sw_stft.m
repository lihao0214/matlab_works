clc;
clear all;
close all;
addpath('./aec_record');
addpath('./Common');
sim = 0;
Fs = 8000;
%%
if ~sim,
    fid = fopen('mic_signal.pcm', 'rb'); % near_end
    sid = fopen('speaker_signal.pcm', 'rb'); % far_end
    [nn, fs] = wavread('./FS16K/white');
    nn = resample(nn, 8000, fs);
    near_end = fread(fid, 'int16');
    far_end = fread(sid, 'int16');
    E_n = norm(near_end).^2;
    E_b = norm(nn(1:length(near_end))).^2;
    SNR = 10;
    SNR_c = 10*log10(E_n./E_b);
    SNR_d = SNR-SNR_c;
    nn = nn(1:length(near_end)).*10^(-SNR_d/20);
    near_end = near_end + nn;
    fclose(sid);
    fclose(fid);
else
    [near_end, fs] = wavread('.\CHINESE(MANDARIN)\Ch_f2');
    near_end = resample(near_end, Fs, fs);
    ISM_RIR_bank(my_ISM_setup, 'ISM_RIRs.mat');
    AuData_near = ISM_AudioData('ISM_RIRs.mat', near_end);
    [far_end, fs] = wavread('.\CHINESE(MANDARIN)\Ch_m5');
    far_end = resample(far_end, Fs, fs);
    ISM_RIR_bank(my_ISM_setup_i, 'ISM_RIRs_i.mat');
    AuData_far = ISM_AudioData('ISM_RIRs_i.mat', far_end);
    idx = min(length(AuData_near), length(AuData_far));
    AuData_far = AuData_far(1:idx);
    AuData_near = AuData_near(1:idx);
    near_end = [AuData_far; AuData_near; AuData_near+AuData_far];
    near_end = round(near_end.*2^15);
    far_end = [far_end(1:idx); zeros(idx, 1); far_end(1:idx)];
    far_end = round(far_end.*2^15);
end
%%
N = 256;
R = N/8; % filter length = 7L+1, overlap 7L, 8L-7L=L
K = 5;
M = 7*K-6;
B = 32; % # of mdf blocks
NN = 64;
alpha = 0.995;
as = 0.995;
fL = 300;
fH = 3300;
BandRange = 1+[round(fL/(8000/N)):round(fH/(8000/N))];
ITER = R*floor(length(near_end)/R);
wins = sqrt(hann(N));
buffer_far = zeros(N,1);
buffer_near = zeros(N,1);
buffer_fend = zeros(NN,1);
buffer_nend = zeros(NN,1);
Xfm = zeros(1+N/2,M);
Xf = zeros(1+N/2,M);
Xm = zeros(1+NN/2,B);
Phi_ss = ones(1+N/2,M);
Phi_vs = ones(1+N/2,K);
gain = ones(1+N/2, 1);
See = ones(1+N/2,1);
ppsnr = ones(1+N/2, 1);
lambda = ones(1+N/2, 1);
H = zeros(1+N/2, 1);
Wk = 0.05.*ones(1+N/2, K); % NLMS
Wb = zeros(1+NN/2, B);
ola_buffer = zeros(N,1);
out = zeros(ITER+R-1,1);
output = zeros(ITER+R-1,1);
zut = zeros(ITER+R-1,1);
tmp = fft(hann(N), N);
dc_gain = tmp(1)./R;
dbg = zeros(ITER/R, 1);
nk = ones(1+N/2,1);
ek = ones(1+N/2,1);
af = 0.5;
aa = 0.97;
gs = 0.8125;
mk = 0.6.*ones(K,1);
pn = zeros(1+N/2,1);
pp = zeros(1+NN/2,1);
temp = [];
ERLE = [];
ERLI = [];
threshold = 2e-6;
q = zeros(K,1);
s = 0;
g = 0.99;
alpha_d = 0.99;
for m = 1:R:ITER,
    buffer_far = [buffer_far(1+R:end); far_end(m:m+R-1)];
    buffer_near = [buffer_near(1+R:end); near_end(m:m+R-1)];
    buffer_fend = [buffer_fend(1+R:end); far_end(m:m+R-1)]; %mdf
    buffer_nend = [buffer_nend(1+R:end); near_end(m:m+R-1)]; %mdf
    ff = fft(buffer_fend, NN);
    nn = fft(buffer_nend, NN);
    ww = fft(wins.*buffer_far, N);
    xx = fft(buffer_far, N);
    yy = fft(wins.*buffer_near, N);
    Xm = [Xm(:,2:B), ff(1:1+NN/2)];
    Xfm = [Xfm(:,2:M), xx(1:1+N/2)];
    Xf = [Xf(:,2:M), ww(1:1+N/2)];
    Phi_ss(:,1:M-1) = Phi_ss(:,2:M);
    Phi_ss(:,M) = as.*Phi_ss(:,M)+(1-as).*diag(Xf(:,M)*Xf(:,M)');
    Dk = sum(Wb.*Xm(:,B:-1:1), 2);
    Yk = sum(Wk.*Xfm(:,M:-7:1), 2); %Wk[descending order]
    dn = real(ifft([Dk; conj(Dk(NN/2:-1:2))], NN));
    yn = real(ifft([Yk; conj(Yk(N/2:-1:2))], N));
    en = near_end(m:m+R-1)-yn(end-R+1:end);
    rn = near_end(m:m+R-1)-dn(end-R+1:end);
    ux = (1-aa).*K.*diag(Xfm(:,M)*Xfm(:,M)');
    id = find(ux>=pn);
    ii = find(ux<pn);
    pn(ii) = as.*pn(ii)+ux(ii);
    pn(id) = ux(id);
    % pn = as.*pn+(1-as).*K.*diag(Xfm(:,M)*Xfm(:,M)');
    pp = as.*pp+(1-as).*B.*diag(Xm(:,B)*Xm(:,B)');
    rk = fft([zeros(NN-R,1); rn], NN);
    rk2 = (rk(1:1+NN/2))./(1e-3+pp);
    absrf = max(abs(rk2), threshold);
    absrf = ones(1+NN/2,1)*threshold./absrf;
    rk2 = rk2.*absrf;
    Rk = fft([zeros(N-R,1); en], N);
    Rk2 = (Rk(1:1+N/2))./(1e-3+pn);
    absEf = max(abs(Rk2), threshold);
    absEf = ones(N/2+1,1)*threshold./absEf;
    Rk2 = Rk2.*absEf;

    %     for k = 1:1:K,
    %         q(k,1) = sum(abs(Wk(:,k)),1);
    %     end
    % s = g.*s+(1-g).*abs(en(R));
    % mk = (q+1/(eps+s))./(sum(q,1));
    for k = 1:1:K,
        temp = ifft(Rk2(1:1+N/2).*conj(Xfm(:,7*(K-k)+1)), N);
        % temp = [temp(1:N-R); zeros(R,1)];
        temp = fft(temp, N); % constrainted operations
        Wk(:,k) = Wk(:,k)+mk(k).*temp(1:1+N/2);
    end

    for k = 1:1:B,
        temp = ifft(rk2(1:1+NN/2).*conj(Xm(:,(B-k)+1)), NN);
        temp = [temp(1:NN-R); zeros(R,1)];
        temp = fft(temp, NN);
        Wb(:,k) = Wb(:,k)+0.6.*temp(1:1+NN/2);
    end

    Vk = yy(1:1+N/2);
    Ek = zeros(1+N/2,1);
    for k = K:-1:1,
        if 0,
            Phi_vs(:,k) = as.*Phi_vs(:,k)+(1-as).*diag(Vk*median(Xf(:,7*(k-1)+1:-1:7*(k-2)+2),2)');
        else
            Phi_vs(:,k) = as.*Phi_vs(:,k)+(1-as).*diag(Vk*Xf(:,7*(k-1)+1)');
        end
        H(:,k) = Phi_vs(:,k)./(1+Phi_ss(:,7*(k-1)+1));
        Vk = Vk-H(:,k).*Xf(:,7*(k-1)+1);
        Ek = Ek+H(:,k).*Xf(:,7*(k-1)+1);
    end
    % ek = af.*ek+(1-af).*abs(Ek);
    % nk = af.*nk+(1-af).*abs(Vk);
    % t = sum(sum(abs(H).^2, 2), 1);
    % temp = [temp; t];
    thresh = abs(yy(1:1+N/2));
    Vk = Vk./(eps+abs(Vk)).*min(thresh, abs(Vk));
    See = gs.*See+(1-gs).*diag(Vk*Vk');
    psnr = min(diag(yy(1:1+N/2)*yy(1:1+N/2)')./(eps+See), 8192);
    xi = alpha_d.*(gain.^2.*ppsnr)+(1-alpha_d).*max(psnr-1, 0);
    ppsnr = psnr;
    gain = xi./(1+xi);
    Bk = gain.*yy(1:1+N/2)-Ek;
    lambda = gs.*lambda+(1-gs).*diag(Bk*Bk');
    as = alpha.^(min(sum(lambda(BandRange))./(eps+sum(See(BandRange))), 1));
    dbg(1+floor(m/R)) = min(sum(lambda(BandRange))./(eps+sum(See(BandRange))), 1);
    v = wins.*ifft([Vk; conj(Vk(N/2:-1:2))], N);
    ola_buffer = [ola_buffer(end-(N-R)+1:end)+v(1:N-R); v(end-R+1:end)];
    out(m:m+R-1) = ola_buffer(1:R)./(dc_gain);
    zut(m:m+R-1) = rn;
    output(m:m+R-1) = en;
    ERLE = [ERLE; 10*log10((near_end(m:m+R-1)'*near_end(m:m+R-1))/(eps+en'*en))];
end

out = [out(225:end); zeros(224,1)];

for m = 1:R:ITER,
    ERLI = [ERLI; 10*log10((near_end(m:m+R-1)'*near_end(m:m+R-1))/(eps+out(m:m+R-1)'*out(m:m+R-1)))];
end
