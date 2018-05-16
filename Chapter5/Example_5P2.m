% Example 5.2    Cascade QMF Analysis and Synthesis Filters
%
% Reference:     R. E. Crochiere, L. R. Rabiner, Multirate Digital Signal Processing, 
%                Englewood Cliffs, NJ: Prentice Hall, 1983
%
% by Lee, Gan, and Kuo, 2008
% Subband Adaptive Filtering: Theory and Implementation
% Publisher: John Wiley and Sons, Ltd

addpath '..\Common';          % Functions in Common folder
clear all;close all; 

% QMF 32C
% No of taps = 32,
% Transition bandwidth = 0.0625
% Reconstruction error = 0.009 dB
% Stopband attenuation = 51 dB

% Define first half (16) of symmetric filter coefficients
h0_32C_half = [.69105790e-3 -0.14037930e-2 -0.12683030e-2 .4234195e-2 .1414246e-2 -0.9458318e-2 -0.1303859e-3 .1798145e-1 -0.4187483e-2 ...
    -0.3123862e-1 .1456844e-1 .5294745e-1 -0.3934878e-1 -0.9980243e-1 .1285579 .46640530];

% Create both lowpass and highpass symmetric FIR filter from 
%  the first half of coefficients
h0_32C = [h0_32C_half fliplr(h0_32C_half)];     % Lowpass filter
h1_32C = ((-1).^(0:length(h0_32C)-1)).*h0_32C;  % Highpass filter
figure;                                 % Plot frequency response
freqz(h0_32C); hold on;
freqz(h1_32C);
title('Frequency response of h0 in QMF 32C');

% Frequency response

DFTpoint = 4096;                        % FFT size
[Hz_lp,w] = FreqResp(h0_32C,DFTpoint);  % Frequency response of lowpass filter
[Hz_hp,w] = FreqResp(h1_32C,DFTpoint);  % Frequency response of hignpass filter
Hz = [Hz_lp; Hz_hp];
figure; subplot(3,1,[1 2]), hold on;    % Plot frequency responses
for k = 1:2
    if mod(k,2)==0
        plot(w/pi,20*log10(abs(Hz(k,:))+eps),':');
    else
        plot(w/pi,20*log10(abs(Hz(k,:))+eps));
    end
end
axis([0 1 -120 10]); box on; ylabel('Gain (dB)'); 
title('Magnitude response of QMF filter bank');

% Distortion function

Tz = sum(abs(Hz).^2);
subplot(3,1,3); plot(w/pi,10*log10(Tz));
xlabel('Frequency,\omega (\pi)'); ylabel('Gain (dB)');
axis([0 1 -0.01 0.01]);
title('Distortion of the QMF filter');

% (a) Plot a uniform filter bank
B1 = h0_32C;                           % Lowpass filter
B2 = h1_32C;                           % Highpass filter
B10 = zeros(1, 2*length(B1));
B10(1: 2: length(B10)) = B1;
B11 = zeros(1, 2*length(B2));
B11(1: 2: length(B11)) = B2;
C0 = conv(B1, B10);C1 = conv(B1, B11);
C2 = conv(B2, B10);C3 = conv(B2, B11);

% Compute the frequency responses of dgital filters
N=1024;                                % FFT size
[H00z, w] = freqz(C0, 1, N);
h00 = abs(H00z);M00 = 20*log10(h00);
[H01z, w] = freqz(C1, 1, N);
h01 = abs(H01z);M01 = 20*log10(h01);
[H10z, w] = freqz(C2, 1, N);
h10 = abs(H10z);M10 = 20*log10(h10);
[H11z, w] = freqz(C3, 1, N);
h11 = abs(H11z);M11 = 20*log10(h11);

figure;
plot(w/pi, M00,'r-','LineWidth',2); hold on;
plot(w/pi, M01,'g--','LineWidth',2);
plot(w/pi, M10,'b--','LineWidth',2);
plot(w/pi,M11,'c-','LineWidth',2);
legend('Band #1','Band #2','Band #3','Band #4','Location','Best');
xlabel('\omega/\pi'); ylabel('Gain, dB');grid
title('Uniform 4 subbands based on QMF');

% (b) Plot a 3-band octave filter bank
B10 = zeros(1, 2*length(B1));
B10(1: 2: length(B10)) = B1;
B11 = zeros(1, 2*length(B2));
B11(1: 2: length(B11)) = B2;
C0 = conv(B1, B10);
C00 = zeros(1,2*length(C0));
C00(1: 2: length(C00)) = C0;
C1 = conv(B1, B11);
C10 = zeros(1,2*length(C1));
C10(1: 2: length(C10)) = C1;
C20 = zeros(1,2*length(B11));
C20(1: 2: length(C20)) = B11;
C20 = [zeros(1,16) C20];

% Compute the frequency responses of digital filters
N=1024;                                 % FFT size
[H00z, w] = freqz(C00, 1, N);
h00 = abs(H00z);
M00 = 20*log10(h00);
[H01z, w] = freqz(C10, 1, N);
h01 = abs(H01z);
M01 = 20*log10(h01);
[H10z, w] = freqz(B11, 1, N);
h10 = abs(H10z);
M10 = 20*log10(h10);

figure;
plot(w/pi, M00,'r-','LineWidth',2); hold on;
plot(w/pi, M01,'g--','LineWidth',2);
plot(w/pi, M10,'b--','LineWidth',2);
legend('Band #1','Band #2','Band #3','Location','Best');
xlabel('\omega/\pi'); ylabel('Gain, dB');grid
title('Octave 3 subbands based on QMF');
