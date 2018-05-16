clc;
clear all;
close all;
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
    E_n = norm(near_end).^2;
    E_b = norm(nn(1:length(near_end))).^2;
    SNR = inf;
    SNR_c = 10*log10(E_n/E_b);
    SNR_d = SNR-SNR_c;
    nn = nn(1:length(near_end))*10^(-SNR_d/20);
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
    near_end = awgn(near_end, Inf, 'measured');
    near_end = round(near_end*2^15);
    far_end = [far_end(1:idx); zeros(idx, 1); far_end(1:idx)];
    far_end = round(far_end*2^15);
end

%%
fid = fopen('test_out_i','rb');
vk_dbg = fread(fid,'int32');
fclose(fid);
vk_r = vk_dbg(1:2:end);
vk_i = vk_dbg(2:2:end);

fid = fopen('test_out_iv','rb');
yk_dbg = fread(fid,'int32');
fclose(fid);
yk_r = yk_dbg(1:2:end);
yk_i = yk_dbg(2:2:end);

M = 4;
Q = 2;
N = 256;
R = 128;
fL = 300;
fH = 3300;
BandRange = 1+[round(fL/(8000/N)):round(fH/(8000/N))];
alphaS = 0.98^((8000*(N/2))/128/Fs); % B = N/2
alphas = alphaS;
gamma = (0.25)^((8000*(N/2))/128/Fs);
eps = 1;
ITER = R*floor(length(near_end)/R);
buffer = zeros(N, 1);
buffer_y = zeros(N, 1);
wins=hann(N,'period');
dcgain=R/sum(wins,1);
wins = sqrt(wins*dcgain);
XFm = zeros(1+N/2, M);
XFm1=zeros(1+N/2,M);
XFm2=zeros(1+N/2,M);
XFp1=zeros(1+N/2,M);
XFp2=zeros(1+N/2,M);
S_ss = ones(1+N/2, M);
S_vs = ones(1+N/2, M);
S_vm = ones(1+N/2, M);
S_vm2 = ones(1+N/2, M);
S_vp = ones(1+N/2, M);
S_vp2 = ones(1+N/2, M);
S_xx = ones(1+N/2, Q);
S_xy = ones(1+N/2, Q);
lambda = ones(1+N/2, 1);
e_k = ones(1+N/2, 1);
v_k = ones(1+N/2, 1);
See = ones(1+N/2, 1);
ppsnr = ones(1+N/2, 1);
gain = ones(1+N/2, 1);
V_k = zeros(1+N/2, 1);
ola_buffer = zeros(N-R, 1);
olv_buffer = zeros(N-R, 1);
out = zeros(ITER+R-1, 1);
outv = zeros(ITER+R-1, 1);
alpha_d = 0.96.^((8000*(N/2))./128./Fs);
dbg = [];
h = zeros(1+N/2, M);
hr = zeros(1+N/2, M);
t = []; % debug
g = []; % debug
XFq = zeros(1+N/2,Q);
DBL=[];

for m = 1:R:ITER,
    idx = 1+(1+N/2)*floor((m-1)/R);
    buffer = [buffer(1+R:N); far_end(m:m+R-1)]; % time align within one block
    buffer_y = [buffer_y(1+R:N); near_end(m:m+R-1)];
    % DBL=[DBL;buffer_y];
    % xx = fft(buffer.*wins, N);
    xx = yk_r(idx:idx+N/2)+j.*yk_i(idx:idx+N/2);
    XFm = [XFm(:,2:end), xx(1:1+N/2)];
    % XFm1=[XFm1(:,2:end),circshift(xx(1:1+N/2),1)];
    % XFm2=[XFm2(:,2:end),circshift(xx(1:1+N/2),2)];
    % XFp1=[XFp1(:,2:end),circshift(xx(1:1+N/2),-1)];
    % XFp2=[XFp2(:,2:end),circshift(xx(1:1+N/2),-2)];
    % yy = fft(buffer_y.*wins, N);
    yy = vk_r(idx:idx+N/2)+j.*vk_i(idx:idx+N/2);
    V_k = yy(1:1+N/2);
    O_k = zeros(1+N/2, 1);
    S_ss(:,1:M-1) = S_ss(:,2:M);
    S_ss(:,M) = alphas.*S_ss(:,M)+(1-alphas).*diag(XFm(:,M)*XFm(:,M)');
    % t = [t;xx(1:1+N/2)];
    for k = M:-1:1,
        % Thresh = abs(V_k);
        S_vs(:,k) = alphas.*S_vs(:,k)+(1-alphas).*diag(V_k*XFm(:,k)');
        h(:,1) = S_vs(:,k)./max(eps,S_ss(:,k));
        V_k = V_k-XFm(:,k).*h(:,1);
        O_k = O_k+XFm(:,k).*h(:,1);
        % S_vm(:,k) = alphas.*S_vm(:,k)+(1-alphas).*diag(V_k*XFm1(:,k)');
        % h(:,2) = S_vm(:,k)./max(eps,circshift(S_ss(:,k),1));
        % V_k = V_k-XFm1(:,k).*h(:,2);
        % O_k = O_k+XFm1(:,k).*h(:,2);
        % S_vp(:,k) = alphas.*S_vp(:,k)+(1-alphas).*diag(V_k*XFp1(:,k)');
        % h(:,3) = S_vp(:,k)./max(eps,circshift(S_ss(:,k),-1));
        % V_k = V_k-XFp1(:,k).*h(:,3);
        % O_k = O_k+XFp1(:,k).*h(:,3);
        % absVf = max(abs(V_k), Thresh);
        % absVf = Thresh./max(eps,absVf);
        % V_k = V_k.*absVf;
    end
    Thresh = abs(yy(1:1+N/2));
    absVf = max(abs(V_k), Thresh);
    absVf = Thresh./max(eps,absVf);
    V_k = V_k.*absVf;

    temp = sum(sum(abs(h), 2), 1);
    t = [t; h(:,1)];
    g=[g;V_k(2:N/2)];
    % dbg = [dbg;V_k];
    Py=diag(yy(1:1+N/2)*yy(1:1+N/2)');
    if m==1,
        parameters = initialise_parameters(Py/2^30,Fs,'imcra');
    else
        parameters = noise_estimation(Py/2^30,'imcra',parameters);
    end

    XFq(:,1:Q-1)=XFq(:,2:Q);
    XFq(:,Q)=diag(xx(1:1+N/2)*xx(1:1+N/2)');
    S_xx(:,1:Q-1)=S_xx(:,2:Q);
    S_xx(:,Q)=alphas.*S_xx(:,Q)+(1-alphas).*diag(XFq(:,Q)*XFq(:,Q)');
    Yt = diag(V_k*V_k');
    eK = Yt;
    uK = zeros(1+N/2,1);
    for k = 1:1:Q,
        S_xy(:,k)=alphas.*S_xy(:,k)+(1-alphas).*diag(eK*XFq(:,1+Q-k)');
        hr(:,k)=S_xy(:,k)./max(1,S_xx(:,1+Q-k));
        eK = max(eK-hr(:,k).*XFq(:,1+Q-k),0);
        uK = uK+hr(:,k).*XFq(:,1+Q-k);
    end
    cK = max(0,1-1.2*sqrt(uK./max(eps,Yt))).*V_k;
    e_k = (1-50/128).*e_k+(50/128).*abs(O_k);
    v_k = (1-1/16).*v_k+(1/16).*abs(V_k);
    % dbg=[dbg;e_k(2:N/2)];
    See = gamma.*See + (1-gamma).*diag(V_k*V_k');
    psnr = diag(yy(1:1+N/2)*yy(1:1+N/2)')./max(eps,See);
    % dbg = [dbg; [0;diag(yy(2:N/2)*yy(2:N/2)');0]];
    % dbg = [dbg; [psnr(2:N/2)]];
    xi = alpha_d.*(gain.^2.*ppsnr)+(1-alpha_d).*max(psnr-1, 0);
    % dbg = [dbg;xi(2:N/2)];
    ppsnr = psnr;
    gain = xi./(1+xi);
    B_k = gain.*yy(1:1+N/2)-O_k;
    lambda = gamma.*lambda+(1-gamma).*diag(B_k*B_k');
    r=max(min(sum(lambda(BandRange),1)/max(sum(See(BandRange),1),eps),1),63/32768);
    % alphas = alphaS.^r;
    % dbg=[dbg;alphas];
    dbg=[dbg;lambda(2:N/2)];
    gainp = max(1-((16*r).*abs(e_k))./max(eps,abs(v_k)), 0);
    F_k = V_k.*gainp;

    v = real(ifft([F_k; conj(F_k(N/2:-1:2))]));
    e = real(ifft([V_k; conj(V_k(N/2:-1:2))]));
    v = v.*wins;
    e = e.*wins;
    ola_buffer=ola_buffer+v(1:N-R);
    olv_buffer=olv_buffer+e(1:N-R);
    out(m:m+R-1) = ola_buffer(1:R);
    outv(m:m+R-1) = olv_buffer(1:R);
    ola_buffer=[ola_buffer(1+R:N-R);v(N-R+1:N)];
    olv_buffer=[olv_buffer(1+R:N-R);e(N-R+1:N)];
end
