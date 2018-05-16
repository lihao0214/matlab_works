clc;
clear all;
close all;
%%

fid = fopen('mic_signal.pcm');
sid = fopen('speaker_signal.pcm');
near_end = fread(fid, 'short');
far_end = fread(sid, 'short');
fclose(sid);
fclose(fid);
%% prototype filter

c = -Table_analysis_window;
idx = [1:1:length(c)];
I = find(mod(floor((idx-1)/64), 2) == 1);
p = c;
p(I) = -p(I);

hann_size = 512; % partition modulated prototype filter into common item
K = 32;
L = 384;
R = 32; % cos((2*k+1)*(r-16)*pi/64)
M = zeros(32, 64);
for k = 0:1:31,
    for r = 0:1:63,
        M(k+1, r+1) = cos((2*k+1)*(r-16)*pi/64);
    end
end

%% synthesis Matrix
N = zeros(64, 32); % cos((2k+1)(n+16)*pi/64) = cos((2k+1)(64*q+r+16)*pi/64), cos((2k+1)(r+16)*pi/64)
for r = 0:1:2*R-1,
    for k = 0:1:31,
        N(r+1, k+1) = cos((2*k+1)*(16+r)*pi/64);
    end
end

%% analysis part
out = [];
weight = zeros(K, L/R);
epsi = 1e-5 * 2.^30;
mu = 0.075;
ITER = floor(length(near_end)/R) * R;
buffer_n = zeros(hann_size, 1);
buffer_f = zeros(hann_size, 1);
buffer_sub_filt = zeros(K, L/R);
% buffer_v = zeros(64, 16);
L1 = (512/R) * 64;
buffer_v = zeros(L1, 1);

for m = 1:R:ITER,
    buffer_n(end:-1:1+R) = buffer_n(end-R:-1:1); % conv, h_{k}[m] * x[n-m]
    buffer_n(R:-1:1) = near_end(m:1:m+R-1);
    buffer_f(end:-1:1+R) = buffer_f(end-R:-1:1);
    buffer_f(R:-1:1) = far_end(m:1:m+R-1);
    zn = c .* buffer_n;
    for r = 0:1:63,
        yn(r+1, 1) = sum(zn(r+1:64:end));
    end
    for k = 0:1:31,
        Sn(k+1, 1) = M(k+1, :) * yn;
    end
    zf = c .* buffer_f;
    for r = 0:1:63,
        yf(r+1, 1) = sum(zf(r+1:64:end));
    end
    for k = 0:1:31,
        Sf(k+1, 1) = M(k+1, :) * yf;
    end
    buffer_sub_filt = [buffer_sub_filt(:, 2:end), Sf];
    Sy = diag(buffer_sub_filt * weight');
    Se = Sn - Sy;
    for k = 1:1:32,
        weight(k, :) = weight(k, :) + mu .* Se(k, 1)' .* buffer_sub_filt(k, :) ./ (buffer_sub_filt(k,:) * buffer_sub_filt(k,:)' + epsi);
    end
    for r = 0:1:63,
        v(r+1, 1) = N(r+1,:) * Sn;
    end

    %  buffer_v(:, end:-1:2) = buffer_v(:, end-1:-1:1); % 16-depth fifo
    %  buffer_v(:, 1) = v;
    %  for n = 0:1:7,
    %    u(1+64*n:1+64*n+31, 1) = buffer_v(1:32, 1+2*n);
    %    u(1+32+64*n:1+32+64*n+31, 1) = buffer_v(33:64, 1+1+2*n);
    %    w(1:32, 1+2*n) = u(1+64*n:1+64*n+31) .* c(1+64*n:1+64*n+31);
    %    w(1:32, 1+1+2*n) = u(1+32+64*n:1+32+64*n+31) .* c(1+32+64*n:1+32+64*n+31);
    %  end
    buffer_v(end:-1:65, 1) = buffer_v(end-64:-1:1, 1);
    buffer_v(1:1:64, 1) = v;
    for n = 0:1:512/R/2-1,
        u(1+2*R*n:1+2*R*n+R-1, 1) = buffer_v(1+128*n:1+128*n+R-1);
        u(1+R+2*R*n:1+R+2*R*n+R-1, 1) = buffer_v(1+128*n+(128-R):1+128*n+(128-R)+R-1);
        w(1:R, 1+2*n) = u(1+2*R*n:1+2*R*n+R-1) .* c(1+2*R*n:1+2*R*n+R-1);
        w(1:R, 1+1+2*n) = u(1+R+2*R*n:1+R+2*R*n+R-1) .* c(1+R+2*R*n:1+R+2*R*n+R-1);
    end

    out = [out, sum(w')];
end
