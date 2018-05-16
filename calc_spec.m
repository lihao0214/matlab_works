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


function y = calc_spec(signal,fbankopts)
%%
fft_len = fbankopts.fft_len;
data_len = fbankopts.data_len;
nshift = fbankopts.nshift;
window = fbankopts.window;
%%


if strcmp(window,'han')
    window = hanning(data_len,'periodic');
end
if strcmp(window,'ham')
    window = hamming(data_len,'periodic');
end

if strcmp(window,'rect')
    window = ones(data_len,1);
end

%%

window = window(:);

signal = [signal(:); zeros(data_len,1)];
sig_length = length(signal);

count = fix((sig_length-data_len)/nshift);

y = zeros(data_len,count);

for ii=1:count
    current = signal((ii-1)*nshift+1:(ii-1)*nshift+data_len);
    y(:,ii) = (current.*window);
end
y  = fft(y,fft_len);
y = y(1:fft_len/2+1,:);
