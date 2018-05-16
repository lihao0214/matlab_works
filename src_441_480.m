clc;
close all;
clear all;

%%
L = 264;
Rs = 60;
a = fir1(L-1, 1/2, chebwin(L, Rs));
% a = firnyquist(L, 2);
% a = [a, 0];
h(:,1) = a(1:2:end);
h(:,2) = a(2:2:end);

% [near_end, fs] = wavread('.\CHINESE(MANDARIN)\Ch_f1');
% xn = resample(near_end, 44100, 16000);
t = [0:1/44100:1];
xn = chirp(t,0,1,22050,'linear');
% xn = sin(2*pi*17000*[0:1:16383]/44100); % sinusoidal

ITER = length(xn);
buffer = zeros(L/2, 1);
y = zeros(2, 1);
out = [];

for m = 1:1:ITER,
    buffer = [buffer(2:end); xn(m)];

    for k = 1:2,
        y(k) = buffer(end:-1:1)'*h(:,k);
    end
    out = [out; y];
end

iter = 147*floor(length(out)/147);
M = 192;
t = [0:1:M/2*80];
f = sinc(t/80)';
wins = hann(M);
buf = zeros(M, 1);
v = zeros(M+147, 1);
output = [];
for m = 1:147:iter,
    x = out(m:m+147-1);
    v = [v(end-M+1:end); x]; %[M-samples, new_frame]
    p = 0; % pos
    for k = 1:1:80,
        j = floor(p); % integer
        s = 1+(M/2+j);
        frac = p - j;
        fc = floor(frac * 80 + 0.5);
        fp = 80 - fc;
        idx = [1+fc:80:1+fc+80*(M/2-1)]; % M/2
        inx = [1+fp:80:1+fp+80*(M/2-1)]; % M/2
        buf = [v(s-M/2:s+M/2-1)];
        h = [flipud(f(idx)); f(inx)].*wins;
        g = sinc(frac+[M/2:-1:-M/2+1])'.*wins;
        output = [output; g'*buf];
        p = p + 147/80;
    end

    %     p = 1;
    %     for k = 0:1:79,
    %         ak = k*147/80-floor(k*147/80);
    %         sk = floor(k*147/80)-floor((k-1)*147/80);
    %         buf = [buf(end-(M-sk)+1:1:end); x(p:p+sk-1)];
    %         h = sinc(ak+[M/2:-1:-M/2+1]);
    %         nn = [nn; p+sk-1];
    %         output = [output; h*buf];
    %         p = p + sk;
    %     end
end
