% Example 5.1   Oversampled Filter Bank Using Relaxed Analysis/Synthesis
%               Filter Design Constraints(Compared With Critically-sample
%               Filter Banks)
%
% by Lee, Gan, and Kuo, 2008
% Subband Adaptive Filtering: Theory and Implementation
% Publisher: John Wiley and Sons, Ltd

addpath '..\Common';      % Functions in Common folder
clear all; close all;

% Prototype filter

Fs = 8000;                % Sampling frequency
Rp = 1;                   % Passband ripple = 1 dB
Rs = 60;                  % Stopband attenuation = -60dB
L = 143;                  % Length of FIR filter
we = 3500*2*pi/Fs;        % End of frequency range of interest
wc = 1.1/16;                 % Cutoff frequency (normalized freq wrt Fs/2)

b = fir1(L,wc, chebwin(L+1,Rs)); % Chebyshev window design method
figure; freqz(b); title('Frequency response of prototype filter (critically sampled)');

% Create Critically-sampled 8-band DFT Filter Banks

nbands = 16;               % Number of filter banks
D = 12;                    % Decimation factor (critically sampled)
len = length(b);
[H,G] = make_bank_DFT(b,nbands);
H = sqrt(D)*H;
G = sqrt(D)*G;
NFFT = 1024;              % FFT size
dist_alias(H',G',nbands); % Compute distortion and aliasing
[H,w,delta_w] = FreqResp(H,NFFT);
Hz = H';

% Plot frequency response

figure; subplot(3,1,[1 2]), hold on;
for k = 1:nbands
    if mod(k,2)==0
        plot(w/pi,20*log10(abs(Hz(k,:))+eps));
    else
        plot(w/pi,20*log10(abs(Hz(k,:))+eps),':');
    end
end
axis([0 2 -120 10]); box on; ylabel('Gain (dB)');
title('Frequency response of critically-sampled 8-band DFT filter bank');

% Distortion function

Tz = sum(abs(Hz).^2);
subplot(3,1,3); plot(w/pi,10*log10(Tz));
xlabel('Frequency,\omega (\pi)'); ylabel('Gain (dB)');
axis([0 2 0 20]);
title('Combined distortion response (critically sampled)');

% Create Oversampled 8-band DFT Filter Banks

% Define parameter values

Fs = 8000;                % Sampling frequency
Rp = 1;                   % Passband ripple = 1 dB
Rs = 60;                  % Stopband attenuation = -60dB
L = 255;                  % Order of FIR filter, length = 128
wc = 1/32;                 % Set cutoff frequency at 1/6

% Design a prototype lowpass filter using a Chebyshev window with Rs
% decibels of relative sidelobe attentuation

b_os = fir1(L,wc,chebwin(L+1,Rs));

figure;                   % Plot the prototype lowpass filter
freqz(b_os);              % Compute frequency response
title('Frequency response of prototype filter (oversampled)');

nbands = 64;               % Number of filter banks
D = 32;                    % Decimation factor (oversampled)

len = length(b_os);
[H,G] = make_bank_DFT(b_os,nbands);
% Design DFT filter banks
H = sqrt(D)*H;
G = sqrt(D)*G;
NFFT = 1024;              % FFT size
dist_alias(H',G',nbands); % Compute distortion and aliasing
[H,w,delta_w] = FreqResp(H,NFFT);
% Compute frequency response
Hz = H';
figure;                   % Plot results
subplot(3,1,[1 2]), hold on;
for k = 1:nbands
    if mod(k,2)==0
        plot(w/pi,20*log10(abs(Hz(k,:))+eps));
    else
        plot(w/pi,20*log10(abs(Hz(k,:))+eps),':');
    end
end
axis([0 2 -120 10]); box on; ylabel('Gain (dB)');
title('Frequency response of oversampled 8-band DFT filter bank');

% Distortion function

Tz = sum(abs(Hz).^2);
subplot(3,1,3); plot(w/pi,10*log10(Tz));
xlabel('Frequency,\omega (\pi)'); ylabel('Gain (dB)');
axis([0 2 0 20]);
title('Combined distortion response (oversampled)');


