clc;
clear all;
close all;

%%
framelength = 80;
windowlength = 240;
P = 20;
[signal, fs] = wavread('./CHINESE(MANDARIN)/Ch_f2');
signal = resample(signal, 8000, fs);
signal = awgn(signal, 50, 'measured');
L = length(signal);
FrameNumber = floor(L/framelength)-2;

%% init values
exc = zeros(L,1);
zi_pre = zeros(P,1);
s_rec = zeros(L,1);
exc_syn_t = zeros(L,1);
s_syn_t = zeros(L,1);
last_syn_t = 0;
zi_syn_t = zeros(P,1);
hw = hamming(windowlength);

%%
for n = 3:FrameNumber,
    s_w = signal(n*framelength+1-windowlength:n*framelength).*hw;
    [A,E] = lpc(s_w, P); % A = [1,a2,a3,...aP]
    s_f = signal((n-1)*framelength+1:n*framelength); % \tilde{x[n]}=-a2*x[n-1]-a3*x[n-2]+...-aP*x[n-P+1]
    [exc1,zi_pre] = filter(A,1,s_f,zi_pre); % en[n]=x[n]+a2*x[n-1]+...+aP*x[n-P+1]
    exc((n-1)*framelength+1:n*framelength) = exc1; % E(z)=X(z)*A(z)
    s_Pitch = exc(n*framelength-222:n*framelength); % A(z)=a1*z^(P-1)+a2*z^(P-2)+...+aP*z^(0)
    PT = findpitch(s_Pitch);
    G = sqrt(E*PT);
    PT1 = floor(2*PT);
    poles = roots(A);
    deltaOMG = -150*2*pi/8000;
    for p = 1:P, % complex conj
        if imag(poles(p)) > 0,
            poles(p) = poles(p).*exp(j*deltaOMG);
        elseif imag(poles(p)) < 0,
            poles(p) = poles(p).*exp(-j*deltaOMG);
        end
    end

    A1 = poly(poles);
    tempn_syn_t = [1:n*framelength-last_syn_t]';
    exc_syn1_t = zeros(length(tempn_syn_t), 1);
    exc_syn1_t(mod(tempn_syn_t,PT1) == 0) = G; % impulse
    exc_syn1_t = exc_syn1_t((n-1)*framelength+1-last_syn_t:n*framelength-last_syn_t);
    [s_syn1_t,zi_syn_t] = filter(1,A1,exc_syn1_t,zi_syn_t);
    exc_syn_t((n-1)*framelength+1:n*framelength) = exc_syn1_t;
    s_syn_t((n-1)*framelength+1:n*framelength) = s_syn1_t;
    last_syn_t = last_syn_t + PT1*floor((n*framelength-last_syn_t)/PT1);
end
