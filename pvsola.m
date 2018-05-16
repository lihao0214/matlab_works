function [out,olag,oxcr] = pvsola(in,ratio,C,winsize,hopsize)
% usage : [out,olag,oxcr] = pvsola(in,ratio,C,winsize,hopsize)

% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%
% Copyright - University of Mons, written by Alexis Moinet - 2011-2012
%
% related publication is:
%
% A. MOINET, T. DUTOIT, "PVSOLA: a Phase Vocoder with Synchronized OverLap-Add", Proceedings of the 14th International Conference on Digital Audio Effects (DAFx-11), pp. 269-275, iRCAM, Paris, France, September 19-23, 2011.
%
% and can be downloaded at:
%
% http://recherche.ircam.fr/pub/dafx11/Papers/57_e.pdf
%
% or
%
% http://tcts.fpms.ac.be/publications/papers/2011/dafx2011_pvsola_amtd.pdf

% TODO rewrite this code in a clean and structured way.
%
% yes, this is development/test code that went through many iterations
% before being stable enough and is thus not very readeable nor clean.
% it's already a bit better than what was used for DAFx11 though.


if nargin < 1, help('pvsola'), return, end
if nargin < 2, ratio = 1; end
if nargin < 3, C = 3; end
if nargin < 4, winsize = 512; end % 512 is good for Fs = 16000. adapt as needed.
if nargin < 5, hopsize = round(winsize/4); end

LOCK = 1;
PVSOLA = 0; % if 0 --> normal phase vocoder, if 1 --> PVSOLA
flagreset = 1;% don't touch this

in = in(:);

out = zeros(ceil(length(in)*ratio),1); %could be round instead of ceil, but "better safe than sorry"
olag = zeros(ceil(length(in)/hopsize),1);
oxcr = zeros(ceil(length(in)/hopsize),1);

win = pvsola_hann(winsize);
nrfft = winsize/2 + 1;

K = 1/1.5;

resetwin = zeros(winsize,1);
resetwin(1:hopsize) = win(end-hopsize+1:end).^2;
resetwin(1:2*hopsize) = resetwin(1:2*hopsize) + win(end-2*hopsize+1:end).^2;
resetwin(1:3*hopsize) = resetwin(1:3*hopsize) + win(end-3*hopsize+1:end).^2;

resetwin = resetwin * K;

pos = 1;
start = 1;
stop = winsize;
outstart = 1;
outstop = winsize;
finalstop = length(in)-2*winsize-hopsize;

%NB : changed this
%x1 = 2*winsize+hopsize; %x1 = length(out(outstart-2*hopsize:outstart+winsize+2*hopsize))
x1 = 2*winsize; %x1 = length(out(outstart-2*hopsize:outstart+winsize+2*hopsize))
x2 = 2*winsize+4*hopsize;%+3*hopsize;

lf = in(start:stop);
rf = in(start+hopsize:stop+hopsize);

fftlf = fft(lf.*win);
fftrf = fft(rf.*win);
alf = abs(fftlf);
arf = abs(fftrf);
plf = angle(fftlf);
prf = angle(fftrf);

phi = plf;
dphi = prf-plf;

accerr = 0;
k = 0;
lfsign = 1;
nframes = 0;
framecount = 0;
resetcount = 0;
count = 0;
count2 = 0;
flagplot = 0;

llag = 0;

while pos < finalstop
    count2 = count2+1;

    %     if mod(count2,100)==0
    %         disp(pos/finalstop)
    %     end
    %     disp('loop')

    %     disp(['pos : ' num2str(pos)])

    % the first 'if PVSOLA' is to be used when we want to take into account the
    % shift/lag (so that constraint C is globally respected over the whole sentence)
    % the second 'if' does not take into account the shift/lag and thus the algo.
    % resets a bit more often than asked by C (because shifts are generally
    % backward, even though we try to even things out)
    % in the end I think that the first 'if' (+/- controlled) is of slightly lesser quality than the
    % uncontrolled one but this has to be double-checked
    % IMPORTANT NOTE: this doesn't mean that the second 'if' causes a deviation
    % from the asked speed ratio !!! it just means that overall it resets a bit more
    % often (e.g. every 2.75 frames instead of every 3 frames)
    % IMPORTANT NOTE2: for the moment there is a logical bug with the two 'if'
    % and the line 'pos = outstart/ratio;' (it's working but it's not exactly
    % the intended algorithm, there seems to be some mismatch)
    % NB: '&& outstart > 2*hopsize' is just for the case C==1 which has some problem
    % with negative boundary limits for the 1st iteration of the algo.
    %
    %     if PVSOLA == 1 && framecount >= (C - fix(accerr)) && outstart > 2*hopsize
    %        accerr = accerr - fix(accerr);
    if PVSOLA == 1 && framecount >= C && outstart > 2*hopsize
        framecount = 0;
        resetcount = resetcount + 1;
        % --- COMMENT/REMOVE THIS BLOCK to have a 'normal' pvoc (Dan Ellis')
        % --- or set variable PVSOLA to 0
        % 1. generate samples for next positionS (so that xcorr is not
        % biaised toward negative lags)
        rf = in(start+hopsize:stop+hopsize);
        fftrf = fft(rf.*win);
        arf = abs(fftrf);
        prf = angle(fftrf);

        dphi = prf-plf;
        tmpoutstart = outstart;
        tmpoutstop = outstop;

        for kk=1:5%tF
            fftf = alf .* exp(1i*phi);
            frame = real(ifft(fftf));
            if LOCK == 0
                phi = phi + dphi;
            else
                [peaks, rois] = pvsola_peakpicking(alf,2);
                phi(peaks) = phi(peaks) + dphi(peaks);

                for pk=1:length(peaks)
                    start_roi = rois(pk);
                    if pk == length(peaks)
                        stop_roi = rois(pk+1);
                    else
                        stop_roi = rois(pk+1)-1;
                    end
                    bin_range=[start_roi:stop_roi];
                    index_peaks = bin_range==peaks(pk);
                    bin_range(index_peaks) = [];
                    phi(bin_range) = phi(peaks(pk)) - ( plf(peaks(pk)) - plf(bin_range) );
                end

                phi = [phi(1:nrfft);-phi(nrfft-1:-1:2)];
            end

            out(tmpoutstart:tmpoutstop) = out(tmpoutstart:tmpoutstop) + K*1*frame .* win;
            tmpoutstart = tmpoutstart + hopsize;
            tmpoutstop = tmpoutstop + hopsize;
        end

        % 1bis. move positionning to the frame (i.e. change lf and rf and all related positions)
        rpos = round(pos);
        start = rpos;
        stop = rpos + winsize-1;

        lf = in(start:stop);
        fftlf = fft(lf.*win);
        alf = abs(fftlf);
        plf = angle(fftlf);

        rf = in(start+hopsize:stop+hopsize);
        fftrf = fft(rf.*win);
        arf = abs(fftrf);
        prf = angle(fftrf);

        dphi = prf-plf;


        % 2. compute xcorr to get optimal positionning

        [xcr,lags] = xcorr(out(outstart-2*hopsize:outstart+winsize+2*hopsize),lf(1:3*hopsize).*(win(1:3*hopsize).^2));

        % TODO : change this to look for a max-peak instead of just a max
        % DONE (pvsola_xcorrpeaks)
        [xcval,ind] = max(abs(xcr(x1:x2)));% TODO : max(abs()) + flip
        %[xcval,ind] = max(pvsola_xcorrpeaks(abs(xcr(x1:x2)),2));% TODO : max(abs()) + flip
        % [ind] = pvsola_maxxcorrpeaks(abs(xcr(x1:x2)),2);% optimized function --> much faster (~2x fater)
        % [ind] = mex_getmaxpeak(abs(xcr(x1:x2)),2);% optimized mex function --> a bit faster than pvsola_maxxcorrpeaks (but can be optimized)

        llag = ind-2*hopsize-2;%+hopsize;%-hopsize; %si x1 = winsize+hopsize

        k = k+1;
        olag(k) = llag;

        % 3. shift positions
        accerr = accerr + llag/hopsize;
        outstart = round(outstart + llag);
        outstop = outstart + winsize - 1;

        %take the shift into account --> this implicitely forces the algorithm to
        %cancel the accumulated error in positionning (accerr) due to llag
        %since the algorithm doesn't stop as long as pos < finalstop
        pos = outstart/ratio;

        lfsign = sign(max(xcr(ind+x1-1),eps));

        phi = plf; % reset phase
        flagreset = 1;
        % --- END OF BLOCK

    elseif pos >= start + hopsize
        count = count + 1;
        start = start + hopsize;
        stop = stop + hopsize;

        fftlf = fftrf;
        alf = arf;
        plf = prf;
        lf = rf;

        rf = in(start+hopsize:stop+hopsize);
        fftrf = fft(rf.*win);
        arf = abs(fftrf);
        prf = angle(fftrf);

        dphi = prf-plf;
    else
        %         count = count + 1 ;
        %         disp(count)
    end

    if flagreset == 0
        % normal vocoder
        rpos = pos-floor(pos);
        %fftf = (alf*(1-rpos)+arf*(rpos)) .* exp(1i*phi);
        fftf = alf .* exp(1i*phi);
        frame = real(ifft(fftf));
    else
        % we are going to insert a frame directly from input --> window properly
        flagreset = 0;
        frame = lf.*win;%ones(winsize,1);%
        out(outstart:outstop) = out(outstart:outstop) .* resetwin;
        out(outstop:outstop+3*winsize) = 0;% adjust margin to sample-precision, 3*winsize is just a dirty hack
    end

    if LOCK == 0
        phi = phi + dphi;
    else
        [peaks, rois] = pvsola_peakpicking(alf,2);
        phi(peaks) = phi(peaks) + dphi(peaks);

        for pk=1:length(peaks)
            start_roi = rois(pk);
            if pk == length(peaks)
                stop_roi = rois(pk+1);
            else
                stop_roi = rois(pk+1)-1;
            end
            bin_range=[start_roi:stop_roi];
            index_peaks = bin_range==peaks(pk);
            bin_range(index_peaks) = [];
            phi(bin_range) = phi(peaks(pk)) - ( plf(peaks(pk)) - plf(bin_range) );
        end

        phi = [phi(1:nrfft);-phi(nrfft-1:-1:2)];
    end

    out(outstart:outstop) = out(outstart:outstop) + K*lfsign*frame .* win;
    framecount = framecount + 1;
    nframes = nframes + 1;

    outstart = outstart + hopsize;
    outstop = outstop + hopsize;

    pos = pos + hopsize/ratio;
end


%%% some useful functions %%%

function [peaks,rois] = pvsola_peakpicking(frame,N)
% how many neighbors do we look on both sides ?
if nargin < 2, N = 2; end
% must be > 0
if N <= 0, N = 1; end

L = length(frame);
peaks = zeros(L,1);

for k=1:L
    leftStop = k-N;
    rightStop = k+N;
    if leftStop < 1, leftStop = 1;end
    if rightStop > L, rightStop = L;end

    if frame(k) == max(frame(leftStop:rightStop))
        peaks(k) = 1;
    end
end

peaks = find(peaks==1);

if nargout > 1
    rois = zeros(length(peaks)+1,1);

    % start of first ROI
    rois(1) = 1;
    % end of last ROI
    rois(end) = length(frame);

    for k = 1:length(peaks)-1
        subFrame = frame(peaks(k)+1:peaks(k+1)-1);
        [m,ind] = min(subFrame);
        rois(k+1) = ind + peaks(k);
    end
end

function [peaks,rois] = pvsola_xcorrpeaks(frame,N)

peakspos = pvsola_peakpicking(frame,N);
peaks = zeros(length(frame),1);
peaks(peakspos) = frame(peakspos);
peaks(1:N) = 0;
peaks(end-N+1:end) = 0;

function [pos,val] = pvsola_maxxcorrpeaks(frame,N)

% this can be further optimized but it's already 2x faster than
% pvsola_xcorrpeaks or an exhaustive search

if N<0, N=0; end

L = length(frame);
maxpk = 0;
pos = 1;

k=1+N;
while k <= L-N
    if frame(k) > maxpk % if frame(k) is smaller than current max peak, no need to compare with neighbours
        subframe = frame(k-N:k+N);
        if frame(k) == max(subframe)
            maxpk = frame(k);
            pos = k;
            k = k+N+1;%if it's a peak, no need to search for a peak in the next N samples
        else
            k=k+1;
        end
    else
        k=k+1;
    end
end

function [out] = pvsola_hann(n)
% usage : [out] = pvsola_hann(winsize)
%
% central peak value is 1 whether n is even or odd

if mod(n,2)==0
    out = hann(n+1);
    out = out(1:end-1);
else
    out = hann(n);
end
