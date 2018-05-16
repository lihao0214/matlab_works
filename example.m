%%
% Institute for Signal Processing (University of Luebeck, Germany)
% Copyright (c) 2011 by Radoslaw Mazur
%
% Permission to use, copy, modify, and distribute this software without
% fee is hereby granted FOR RESEARCH/EDUCATION PURPOSES only, provided
% that this copyright notice appears in all copies and in all supporting
% documentation, and that the software is not redistributed for any
% fee (except for a nominal shipping charge).
%
% For any other uses of this software, in original or modified form,
% including but not limited to consulting, production or distribution
% in whole or in part, specific prior permission must be obtained
% from the author.
% Signal processing methods and algorithms implemented by this
% software may be claimed by patents owned by others.
%
% The author makes no representation about the suitability of this
% software for any purpose. It is provided "as is" without warranty
% of any kind, either expressed or implied.
% Beware of the bugs.
%
%     Revision history
%
%     Ver     Date         Description
%     -------------------------------------------------------------------
%     0.5     21-08-2011   basic version
%


% read data

mix1 =wavread('mixture1.wav');
mix2 =wavread('mixture2.wav');


% Options for filterbank
fbank_opts.fft_len=1024*8;
fbank_opts.flen=fbank_opts.fft_len/2+1;
fbank_opts.data_len=1024*2;
fbank_opts.nshift=1024/8;
fbank_opts.window='han';


%calculate the Time-Frequency represen
spec1    = calc_spec(mix1,fbank_opts);
spec2    = calc_spec(mix1,fbank_opts);


iterations = 2000;
mue = 0.02;

% tic
% w_cuda_raw=do_gpu_frq_ica_2d(spec1, spec2, iterations, mue);
% toc

tic
w_matlab_raw=do_matlab_frq_ica_2d(spec1, spec2, iterations, mue);
toc


% compare results
% plot(abs(w_cuda_raw(:) - w_matlab_raw(:)))