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




function [WW px] = do_matlab_frq_ica_2d(spec1, spec2, iterations, mue)
%%

[fbins, mlen] = size(spec1);
observations = 2;

WW = zeros(observations,observations,fbins);

II = eye(observations);

t = clock;

for ii=1:fbins
%%
    if (true)
        if mod(ii,100)==1;
            elap_time = etime(clock,t);
            eta_time = elap_time/ii*(fbins-ii);
            fprintf('separating %3.0f%%    Time: %4.0fs    ETA: %4.0fs \n', ii/fbins*100,elap_time,eta_time)
        end
    end

    
    X=[spec1(ii,:);spec2(ii,:)];

    W = eye(observations);

    for n=1:iterations
        Y=W*X;
        nonlinearF = sign(Y);
        d_tmp = (II - nonlinearF * Y'/mlen)*W;
        W = W + mue*d_tmp;
    end

    WW(:,:,ii) = W;
    %WW(:,:,ii) = diag(diag(W_raw^-1))*W_raw;
end
