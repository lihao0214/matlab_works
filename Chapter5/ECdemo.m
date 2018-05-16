% ECdemo          Echo Cancellation Demonstration Using the 
%                 (i) NLMS, (ii) IPNLMS, and (iii) PMSAF Algorithms
%                 Under Different Sparse Impulse Responses
%
% Reference:      J. Homer, "Quantifying the convergence speed of LMS adaptive 
%                 FIR filter with autoregressive inputs," IEE Electronics Letters, 
%                 vol. 36, No. 6, pp. 585-586, March 2000
%
% Note:           This M-file uses a different impulse response as that in 
%                 Section 5.4.2. Simulation results shown in Section 5.4.2 are 
%                 based on ensemble averaging. 
%                 
% by Lee, Gan, and Kuo, 2008
% Subband Adaptive Filtering: Theory and Implementation
% Publisher: John Wiley and Sons, Ltd



addpath '..\Common';            % Functions in Common folder
clear all; close all;

% Filter and run parameters

M = 512;                        % Length of fullband adaptive filter
ITER = 8.0*80000;               % Total number of iterations                                                                  
truncate = M;                   % Number of samples to be truncated for running-in effect

% Load different room impulse responses

load('SparseImpulse_512_32.mat','ho_psi8','ho_psi20','ho_psi50','ho_psi200');
select_ho = input('Select impulse response, ho: (1) PSI8(sparse), (2) PSI20, (3) PSI50, (4) PSI200(dispersive) : ','s');

if select_ho == '1'
    ho = ho_psi8(1:M);          % Select the sparse impulse response (PSI8)
    psi = 8;                    %   with decay constant = 8 (highly sparse)
elseif select_ho == '2'
    ho = ho_psi20(1:M);         % Select the sparse impulse response (PSI20)
    psi = 20;                   %   with decay constant = 20
elseif select_ho == '3'
    ho = ho_psi50(1:M);         % Select the sparse impulse response (PSI50)
    psi = 50;                   %   with decay constant = 50
else 
    ho = ho_psi200(1:M);        % Select the sparse impulse response (PSI200)
    psi = 200;                  %   with decay constant = 200 (dispersive)
end

figure;
plot(ho); xlabel('Iteration'); ylabel('Amplitude'); grid on;
title(['Impulse response with sparse index, PSI=',num2str(psi)]);

% Generate input and desired signals

select=input('Select excitation signal, un: (1) White noise, (2) Color input, (3) Speech signal: ','s');
if select == '1'
   seed = sum(100*clock);
   randn('state',seed);         % Reset andom generator to diffrent state
   un = randn(1,ITER);
elseif select == '2'            % AR(10) signal (in J.Homer paper)
   seed = sum(100*clock);
   randn('state',seed);         % Reset random generator to diffrent state
   un = randn(1,ITER); 
   load('Homer.mat','a1','a2','a3'); a = a3;              
   un = filter(1,a,un);
else
   load('speech.mat','speech'); un = speech;
   load('speech.mat','speech'); un = [un; speech];
   load('speech.mat','speech'); un = [un; speech];
   load('speech.mat','speech'); un = [un; speech];
   un = un';
   ITER = length(un);
end
    
dn = filter(ho,1,un); 

% Normalization

un = un/std(dn); dn = dn/std(dn); 

% Add white noise to the desired signal

rand('twister',sum(100*clock));                     % Genearate additive noise
%vn = zeros(1,length(dn));                          % Noiseless
%vn = (rand(1,length(dn))-0.5)*sqrt(12*1e-10);      % 100 dB SNR
%vn = (rand(1,length(dn))-0.5)*sqrt(12*1e-6);       % 60 dB SNR
%vn = (rand(1,length(dn))-0.5)*sqrt(12*1e-5);       % 50 dB SNR
vn = (rand(1,length(dn))-0.5)*sqrt(12*1e-4);        % 40 dB SNR
%vn = (rand(1,length(dn))-0.5)*sqrt(12*3.162*1e-4); % 35 dB SNR
%vn = (rand(1,length(dn))-0.5)*sqrt(12*1e-3);       % 30 dB SNR
dn = dn + vn;

% Initialize variables for adaptive filtering

errnlms_sqr = zeros(1,ITER);
erripnlms_sqr = zeros(1,ITER);
errpmsaf_sqr = zeros(1,ITER);
norm_errnlms = zeros(1,ITER);
norm_erripnlms = zeros(1,ITER);
norm_pmsaf = zeros(1,ITER);                         % Length of adaptive filter 
w0 = zeros(M,1);                                    % Intialize filter coefficients to 0                              

% Perform adaptive filtering

% (1) NLMS algorithm

munlms = 0.1;                                       % Step size of NLMS algorithm
alpha_NLMS = 1e-4;                                  % Small constant used in NLMS algorithm     
tic;
disp(sprintf('NLMS, Step size = %.5f',munlms));
Snlms = NLMSinit(w0,munlms);
Snlms.unknownsys = ho;
[ynnlms,ennlms,Snlms] = NLMSadapt(un,dn,Snlms);
disp(sprintf('Total time = %.5f mins',toc/60));

err_sqrnlms = ennlms.^2;                            % Square error
EMLnlms = Snlms.eml.^2;                             % System error norm (normalized)

% (2) IPNLMS algorithm

muipnlms = 0.1;                                     % Step size for IPNLMS algorithm
alpha = 0;
leak = 0;
tic;
disp(sprintf('IPNLMS, Step size = %.5f',muipnlms));
Sipnlms = IPNLMSinit(w0,muipnlms,alpha,leak);
Sipnlms.unknownsys = ho;
[ynipnlms,enipnlms,Sipnlms] = IPNLMSadapt(un,dn,Sipnlms);
disp(sprintf('Total time = %.5f mins',toc/60));
    
err_sqripnlms = enipnlms.^2;                        % Square error
EMLipnlms = Sipnlms.eml.^2;                         % System error norm (normalized)

% (3) PMSAF algorithm

mupmsaf = 0.1;                                      % Step size of PMSAF algorithm
N = 4;                                              % Number of subbands, 4
L = 8*N;                                            % Length of analysis filters, L=2KN, 
                                                    %   overlapping factor K=4
tic;
disp(sprintf('PMSAF: Number of subbands, N = %d, Step size = %.5f',N,mupmsaf));
Spmsaf = PMSAFinit(w0,mupmsaf,N,L);
Spmsaf.unknownsys = ho;
[enpmsaf,Spmsaf] = PMSAFadapt(un,dn,Spmsaf);
disp(sprintf('Total time = %.5f mins',toc/60));

err_sqrpmsaf = enpmsaf.^2;                          % Square error
EMLpmsaf = Spmsaf.eml.^2;                           % System error norm (normalized)

   
% Performance Results

% (i) Final coefficients

figure;
plot(1:length(ho), ho,'r','LineWidth',3);hold on;
plot(1:length(ho), Snlms.coeffs,'b'); 
plot(1:length(ho), Sipnlms.coeffs,'g');
plot(1:length(ho), Spmsaf.coeffs,'m');
legend('Actual','Estimated NLMS','Estimated IPNLMS','Estimated PMSAF','Location','BEST');
title('System identification of FIR filter using LMS, NLMS and PMSAF');grid on;

% (ii) EML plot

figure; hold on; 
plot((0:ITER-1)/1024,10*log10(EMLnlms),'-r', (0:ITER-1)/1024,10*log10(EMLipnlms),'--b',(0:ITER-1)/1024,10*log10(EMLpmsaf),'--m');
legend('EML(NLMS)','EML(IPNLMS)','EML(PMSAF)');
xlabel('Number of iterations (\times 1024 input samples)'); ylabel('Misalignment (dB)');
title('Echo canceller demo in Section 5.4.2');
grid on;axis([0 ITER/1024 (min(10*log10(EMLpmsaf))-10) 10]);

% (iii) Square error signal

figure;
q = 0.99; 
MSEnlms = filter((1-q),[1 -q],err_sqrnlms);
MSEipnlms = filter((1-q),[1 -q],err_sqripnlms);
MSEpmsaf = filter((1-q),[1 -q],err_sqrpmsaf);
hold on; 
plot((0:ITER-1)/1024,10*log10(MSEnlms), (0:ITER-1)/1024,10*log10(MSEipnlms), (0:ITER-1)/1024,10*log10(MSEpmsaf));
legend('MSE(NLMS)','MSE(IPNLMS)','MSE(PMSAF)');
xlabel('Number of iterations (\times 1024 input samples)'); 
ylabel('Mean-square error (with delay)');
title('Echo canceller demo in Section 5.4.2');
grid on; axis([0 ITER/1024 -60 10]);






