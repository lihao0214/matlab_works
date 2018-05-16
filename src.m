%% resample from 32kHz to 96kHz
clc;
clear;
close all;
n = 264;
l = 3;
Rs = 60;
% b = firnyquist(n, l, 0.3);
b = fir1(n-1, 1/l, chebwin(n, Rs));
a = firnyquist(58, 2);
a = [a, 0];
out = [];

%% poly-phase filtering
h(:,1) = b(1:l:end)';
h(:,2) = b(2:l:end)';
h(:,3) = b(3:l:end)';
g(:,1) = a(1:2:end);
g(:,2) = a(2:2:end);

xn = sin(2*pi*15000*[0:1:8191]/32000);
ITER = length(xn);
buffer = zeros(floor(n/l), 1);
y = zeros(l,1);

for m = 1:ITER,
    buffer = [buffer(2:end); xn(m)];
    for k = 1:1:l,
        y(k, 1) = buffer(end:-1:1)'*h(:,k);
    end
    out = [out; y];
end

o = out(1:2:end);
