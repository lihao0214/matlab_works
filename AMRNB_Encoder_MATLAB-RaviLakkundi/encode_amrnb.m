function out_bytes = amrnb_encode(in_speech, out_bytes)
format long e
load  lag_wind.mat
% /* LSF prediction factors (not in MR122) */
load pred_fac.mat
load interpol_window.mat
load qua_gain_code.mat
% /* LSF means ->normalize frequency domain */
load table_gain_lowbitrates.mat
load mean_lsf_3.mat
load dico1_lsf_3.mat
qua_gain_code = qua_gain_code';
table_gain_lowbitrates = table_gain_lowbitrates;
% Main Encode function
global Final_St_80Hz
global working_frame;
global lsp_old;
global past_rq;
global lsp_old_q;
global mem_w_prebig;
global past_in;
global past_out;
global past_frame_160_sample;
global mem_err;
global mem_w0;
global final_exc;
global h1;
global past_qua_en;
global qua_ener_index;
global synth;
% filter the input speech with HighPass filter

%  Down-scaling and high-pass filtering are combined by dividing
%  the coefficients at the numerator by 2.
working_frame(1:80) = working_frame(161:240); % 160
[out_filt_80Hz Final_St_80Hz]=filter([0.927246093 -1.8544941 0.927246903]./2,[1 -1.906005859 0.911376953],in_speech, Final_St_80Hz);

working_frame(81:end) = out_filt_80Hz; % 160
load window_200_40.mat

% Multiply the working_frame with the windows

windowd_input = window_200_40.*working_frame';

% Take autocorr of these values
windowd_input = windowd_input';
windowd_input = [windowd_input zeros(1,11)];
autocorr_r = xcorr(windowd_input,windowd_input, 10);
autocorr_r = autocorr_r(11:end);

autocorr_r = autocorr_r.*lag_wind;
A_unquant = levinson(autocorr_r);
A_unquant_lsp = real(exp(i*poly2lsf(A_unquant)))';

% Interpolate all the lpc's based on old lpc's

A_unquant_subframes_lsp(1,:) = lsp_old.*0.75 + A_unquant_lsp.*0.25;
A_unquant_subframes_lsp(2,:) = lsp_old.*0.5 + A_unquant_lsp.*0.5;
A_unquant_subframes_lsp(3,:) = lsp_old.*0.25 + A_unquant_lsp.*0.75;;
A_unquant_subframes_lsp(4,:) = A_unquant_lsp;


A_unquant_subframes_lpc(1,:) =  lsf2poly(real(acos(A_unquant_subframes_lsp(1,:))));
A_unquant_subframes_lpc(2,:) =  lsf2poly(real(acos(A_unquant_subframes_lsp(2,:))));
A_unquant_subframes_lpc(3,:) =  lsf2poly(real(acos(A_unquant_subframes_lsp(3,:))));
A_unquant_subframes_lpc(4,:) =  A_unquant;


lsp_old = A_unquant_lsp;%(2:end);

% Compute Normalized Lsf

NormLsf = real(acos(A_unquant_subframes_lsp(4,:))) .*(4000.0/pi);

% Compute LSF weighting factors

lsf_wf = filter([1 0 -1],1,NormLsf);
lsf_wf = [lsf_wf(2:end) 4000.0 - NormLsf(9)]; % di

temp_1 =1.8 - (((1.8-1.0)/(1500.0-450.0)).*(lsf_wf - 450.0));
temp_2 =1.8 - (((3.347-1.8)/(450.0-0.0)).*(lsf_wf - 450.0));

temp_1 = temp_1.*temp_1;
temp_2 = temp_2.*temp_2;

ind = find(lsf_wf < 450.0 ); temp_1(ind) = 0;
ind = find(lsf_wf >= 450.0 ); temp_2(ind) = 0;
lsf_wf1 =   temp_1 + temp_2;

% Compute predicted LSF and prediction error

lsf_p = mean_lsf_3 + past_rq.*pred_fac;
lsf_r1 = NormLsf - lsf_p;

% After VQ we will get past_rq
% 3, 3, 4 Split-VQ

dico1_lsf_3=reshape(dico1_lsf_3,3,256);
dist_lsf_3_part1 = sum(((((repmat(lsf_r1(1:3),256,1) - dico1_lsf_3').*(repmat(lsf_wf1(1:3),256,1)))).^2)');
lIndex1 = find(dist_lsf_3_part1 == min(dist_lsf_3_part1));

out_bytes(1) = lIndex1 - 1;

load dico2_lsf_3;

dico2_lsf_3=reshape(dico2_lsf_3,3,512);
dist_lsf_3_part2 = sum(((((repmat(lsf_r1(4:6),512,1) - dico2_lsf_3').*(repmat(lsf_wf1(4:6),512,1)))).^2)');
lIndex2 = find(dist_lsf_3_part2 == min(dist_lsf_3_part2));

out_bytes(2) = lIndex2 - 1;

load dico3_lsf_3;

dico3_lsf_3=reshape(dico3_lsf_3,4,512);
dist_lsf_3_part3 = sum(((((repmat(lsf_r1(7:10),512,1) - dico3_lsf_3').*(repmat(lsf_wf1(7:10),512,1)))).^2)');
lIndex3 = find(dist_lsf_3_part3 == min(dist_lsf_3_part3));

out_bytes(3) = lIndex3 - 1;

past_rq = [dico1_lsf_3(:,lIndex1)' dico2_lsf_3(:,lIndex2)' dico3_lsf_3(:,lIndex3)' ] ; % quantified

lsf1_q = past_rq + lsf_p;

% Check for minimum distance of lsf's to be 50
lsf1_q(find(lsf1_q < 50.0)) = 50.0;

% convert LSFs to the cosine domain
lsp1_q = cos(lsf1_q*(pi/4000.0));

% Now interpolate these quantized lsp's for each subframe
% and get AzQ values

A_quant_subframes_lsp(1,:) = lsp_old_q.*0.75 + lsp1_q.*0.25;
A_quant_subframes_lsp(2,:) = lsp_old_q.*0.5 + lsp1_q.*0.5;
A_quant_subframes_lsp(3,:) = lsp_old_q.*0.25 + lsp1_q.*0.75;;
A_quant_subframes_lsp(4,:) = lsp1_q;


A_quant_subframes_lpc(1,:) =  lsf2poly(real(acos(A_quant_subframes_lsp(1,:))));
A_quant_subframes_lpc(2,:) =  lsf2poly(real(acos(A_quant_subframes_lsp(2,:))));
A_quant_subframes_lpc(3,:) =  lsf2poly(real(acos(A_quant_subframes_lsp(3,:))));
A_quant_subframes_lpc(4,:) =  lsf2poly(real(acos(lsp1_q)));

lsp_old_q = lsp1_q;%(2:end);

% dtx_buffer and check_lsp not implemented

gamma1 = 0.9400024414063;
gamma2 = 0.6000061035156;

for h = 1:2
    wsp_work = [];
    for k = 2*h-1:2*h % 1:2, 3:4
        % process two subframes (which form the "big" subframe) 80 each * 2
        Ap1 = A_unquant_subframes_lpc(k,:).*[1 gamma1.^[1:10]];
        Ap2 = A_unquant_subframes_lpc(k,:).*[1 gamma2.^[1:10]];

        past_in = working_frame(41+((k-1)*40) - 10:41+((k-1)*40) - 1);
        past_in = flipud(past_in')';
        yy = filtic(Ap1,Ap2,flipud(past_out')',past_in);
        wsp = filter(Ap1,Ap2,working_frame(41+((k-1)*40):40 + 40*k), yy);
        past_out = wsp(31:end);
        wsp_work = [wsp_work wsp];

    end

    % Find open loop pitch lag for two subframes
    % [pitch_min - 20, pitch_max - 143]
    % We get only 80 new samples here
    temp1 = wsp_work;
    temp2 = past_frame_160_sample(18:end);
    temp22 = [temp2  temp1];

    corr_wsp = xcorr(temp1,temp22,143);
    corr_wsp = flipud(corr_wsp(21:144)')';
    past_frame_160_sample = [past_frame_160_sample(81:160) wsp_work];

    corr_wsp = [corr_wsp zeros(1,20)]; corr_wsp(find(corr_wsp < 0.0005)) = 0.0000001;

    max1 = max(corr_wsp(1:64));
    lag1_Index = find(corr_wsp(1:64) == max1);
    max1_c = 1/(sqrt(sum(temp22(lag1_Index(1): lag1_Index(1)+80).*temp22(lag1_Index(1): lag1_Index(1)+80))) + 0.0000001) * max1;
    lag1_Index = 144 - lag1_Index(end) + 1;

    max2 = max(corr_wsp(65:104));
    lag2_Index = find(corr_wsp(65:104) == max2);
    max2_c = 1/(sqrt(sum(temp22(65+lag2_Index(1): lag2_Index(1)+80 + 65).*temp22(lag2_Index(1) + 65: lag2_Index(1)+80 + 65))) + 0.0000001) * max2;
    lag2_Index = 144 - 64 - lag2_Index(end) + 1;

    max3 = max(corr_wsp(105:144));
    lag3_Index = find(corr_wsp(105:144) == max3);
    max3_c = 1/(sqrt(sum(temp22(105+lag3_Index(1): 105+lag3_Index(1)+80).*temp22(105+lag3_Index(1): 105+lag3_Index(1)+80))) + 0.0000001) * max3;
    lag3_Index = 144 - 104 - lag3_Index(end) + 1;

    if(max1_c *.85 < max2_c)
        max1_c = max2_c;
        lag1_Index = lag2_Index;
    end
    if(max1_c *.85 < max3_c)
        lag1_Index = lag3_Index;
    end

    %     if(lag1_Index <=40)
    %         lag1_Index = 42;
    %     end
    T_OP(h) = lag1_Index;
end

% Loop for subframe ABS


for index_subframe=1:4
    % Compute the Impulse response of weigted synthesis filter
    % H(z) * W(z) = A(z/g1) / ( A'(z) * A(z/g2) ) [ Every sub-frame]

    Ap1 = A_unquant_subframes_lpc(index_subframe,:).*[1 gamma1.^[1:10]];
    Ap2 = A_unquant_subframes_lpc(index_subframe,:).*[1 gamma2.^[1:10]];

    h1=filter(1,A_quant_subframes_lpc(index_subframe,:),[Ap1 zeros(1,29)]);
    h1=filter(1,Ap2,h1); % Compute the Impulse response of weigted synthesis filter

    work_absframe = working_frame((40*index_subframe + 1):40*index_subframe+40);

    xxx=working_frame((40*index_subframe + 1) - 10:(40*index_subframe));
    xxy=filtic(A_quant_subframes_lpc(index_subframe,:),1,1,flipud(xxx')');
    pred_residual = filter(A_quant_subframes_lpc(index_subframe,:),1,work_absframe,xxy);

    z=filtic(1,A_quant_subframes_lpc(index_subframe,:),mem_err); % Check if you need flipud to mem_err
    error = filter(1,A_quant_subframes_lpc(index_subframe,:),pred_residual,z);

    target = filter(Ap1,1,error, filtic(Ap1,1,1,mem_err));

    z=filtic(1,Ap2,mem_w0); % Check if you need flipud to mem_err
    target = filter(1,Ap2,target,z);

    pred_residual2 = pred_residual;

    % Checking for cl_ltp, closed lop LTP

    % Currently Fractional pitch search is not performed
    % Later will be adding fractional pitch search
    % Assume for a while that the open loop pitch estimation is correct
    % i.e. work on T_OP(1) and T_OP(2)

    % Since we are working only for MR59 mode, Pitch Index is given by
    % For [19 1/3, 84 2/3] the pitch lag is encoded as: [for 1 and 3rd frames]
    % pitch_index = pitch_lag_int*3 + pitch_lag_frac - (19 1/3)*3;
    % pitch_lag_frac = 0, since we are not calculating it
    % for 2 and 4th in range
    % pitch_lag_int  = pitch_index - (85 - 19 1/3)*3 + 85 = pitch_index - 112

    % Ref: http://wiki.multimedia.cx/index.php?title=AMR-NB

    if(index_subframe == 1 || index_subframe == 3)
        if(index_subframe == 3)
            index_temp = 2;
            lag = T_OP(2);
        else
            index_temp = 1;
            lag = T_OP(1);
        end

        if(T_OP(index_temp) <=85)
            pitch_index = (T_OP(index_temp) -1)*3 - 58; % -1 is to match C language indexing, MATLAB Starts from 1
        else
            pitch_index = (T_OP(index_temp) -1) + 112;
        end



    else (index_subframe == 2 || index_subframe == 4)
        % Code using 4 bit for the bitrate we chose, MR59
        pitch_index = 9; % T0*3 - (T0 - 2)*3 + 3; Ref: Standard Code or wiki link, http://wiki.multimedia.cx/index.php?title=AMR-NB
        lag = T_OP(index_subframe/2);
    end

    out_bytes(4 + (index_subframe -1)*4) = pitch_index;

    final_exc(160 + index_subframe*40 - 39:160 + index_subframe*40) = pred_residual;

    % if lag < 40 there is a problem, thus a special case, as exc[-T0]
    % cannot  cover entire excitation use the computed one's

    for tj = 1: 40
        val = ...
            sum(final_exc(160 + index_subframe*40 - 39 - lag + tj - 10:160 + index_subframe*40 - 39 - lag + tj - 10 + 20).*interpol_window);
        val = val + 16384;
        val = floor(val / 2^15);
        final_exc(160 + index_subframe*40 - 40 + tj) = val;
    end

    tempexc =  final_exc(160 + index_subframe*40 - 39:160 + index_subframe*40);%final_exc(160 + index_subframe*40 - 39 - lag + 1 : 160 + index_subframe*40 - 39 - lag  + 40);%

    y1 = conv(tempexc,h1);
    y1 = y1(1:40);

    % Compute pitch gain

    gain_pitch = sum(target.*y1)/(sum(y1.*y1)+0.01) + 0.00001;

    gcoeff = [0 0];
    gcoeff(1) = sum(y1.*y1)+0.01;
    gcoeff(2) = sum(target.*y1); % xn -> target

    if(gain_pitch < 0)
        gain_pitch = 0;
    end
    if(gain_pitch > 1.2)
        gain_pitch = 1.20;
    end

    %out_bytes(5 + (index_subframe -1)*4) = gain_pitch;

    gp_limit = 2.0;

    target2 = target - gain_pitch*y1; % xn2 -> target2
    pred_residual2 = pred_residual2 - gain_pitch*tempexc;

    % Compute Codebook Gain and Vector
    % Already mentioned, this works only for MR59, bitrate 5.9 kbps
    % This code is to be only used as a study material
    % 5.90 kbps mode

    %     * 2 pulse positions coded using 4 and 5 bits
    %     * signs coded using 1 bit for each pulse
    %       Pulse 	Positions
    % i0 	1, 6, 11, 16, 21, 26, 31, 36
    %       3, 8, 13, 18, 23, 28, 33, 38
    % i1 	0, 5, 10, 15, 20, 25, 30, 35
    %       1, 6, 11, 16, 21, 26, 31, 36
    %       2, 7, 12, 17, 22, 27, 32, 37
    %       4, 9, 14, 19, 24, 29, 34, 39

    dn = xcorr(target2,h1);
    dn = dn(40:end);

    % get a sign aray of dn

    sign_arr = zeros(1,40);
    sign_arr(find(dn <0)) = -1;
    sign_arr(find(sign_arr ==0)) = 1;

    if(length(sign_arr) ==0)
        sign_arr = ones(1,40);
    end

    %    *  PURPOSE:  Computes correlations of h[] needed for the codebook search;
    %    *            and includes the sign information into the correlations.
    %    *
    %    *  DESCRIPTION: The correlations are given by
    %    *         rr[i][j] = sum_{n=i}^{L-1} h[n-i] h[n-j];   i>=j; i,j=0,...,L-1
    %    *
    %    *  and the sign information is included by
    %    *         rr[i][j] = rr[i][j]*sign[i]*sign[j]

    xx = h1;
    cor_codemat = zeros(40,40); % phi = H'H
    for ind=1:39
        tc = xcorr(xx(1:40-ind+1),xx(1:40-ind+1));
        cor_codemat(ind,:) = [zeros(1,ind-1) tc(40-ind+1:end)];
        cor_codemat(ind,ind) = cor_codemat(ind,ind)/2;
    end
    cor_codemat(40,:) = [zeros(1,39) xx(1)/2];
    cor_codemat = cor_codemat + cor_codemat';

    % Compute sign array, temporarily in C fashion

    for t1=1:40
        for t2=t1:40
            sign_s(t1,t2)=sign_arr(t1)*sign_arr(t2);
            sign_s(t2,t1) = sign_s(t1,t2);
        end
    end

    cor_codemat = cor_codemat.*sign_s;

    % search_2i40_11bits, StartPos1 and StartPos2: Check the Pulses table
    % first values

    % we can do the below in matrix way

    %     * 2 pulse positions coded using 4 and 5 bits
    %     * signs coded using 1 bit for each pulse
    %       Pulse 	Positions
    % i0 	1, 6, 11, 16, 21, 26, 31, 36  Track_1[1]
    %       3, 8, 13, 18, 23, 28, 33, 38  Track_1[2]
    %
    % i1 	0, 5, 10, 15, 20, 25, 30, 35  Track_2[1]
    %       1, 6, 11, 16, 21, 26, 31, 36  Track_2[2]
    %       2, 7, 12, 17, 22, 27, 32, 37  Track_2[3]
    %       4, 9, 14, 19, 24, 29, 34, 39  Track_2[4]

    codevec=[1 2];
    % We can actually do it by index
    Track_1 = [2 4];    % +1 for Matlab Indexing, as it starts from 1
    Track_2 = [1 2 3 5]; % +1 for Matlab Indexing, as it starts from 1
    psk = -1;
    alpk = 1;
    dn = abs(dn);
    for track1_ind = 1:2
        for track2_ind = 1:4
            % Construct one Algebraic code vector
            ipos = [Track_1(track1_ind) Track_2(track2_ind)];

            for i0 = ipos(1):5:40
                ps0 = dn(i0);
                sq = -1;alp = 1;
                ix = ipos(2);

                for i1 = ipos(2):5:40
                    ps1 = ps0 + dn(i1);
                    alp1= cor_codemat(i0,i0)*0.25 + cor_codemat(i1,i1)*0.25 + cor_codemat(i0,i1)*0.5;
                    sq1 = ps1*ps1;
                    if( sq1/alp1 > sq/ alp)
                        sq= sq1;
                        alp = alp1;
                        ix = i1;
                    end
                end

                if(alpk*sq > psk * alp)
                    psk = sq;
                    alpk = alp;
                    codvec(1) = i0;
                    codvec(2) = ix;
                end

            end
        end
    end

    % Update the target vector and build
    % based on codvec, allocate bits, 4 bits for i0 and 5 bits for i1
    % for example for codevec(1) = 13 (5th Index in table, vertical),
    % codevec(2) = 22 (18th index in table of i1)
    % 13 in binary -> 1101, 18 in binary -> 10010
    % thus the 9 bit number becomes -> 10010 1101 = 301
    % plus two bits for sign, total 11 bits for Algebraic codevector

    f_ind = 0;
    c_sign=[0 0];
    rsign = 0;
    cod = zeros(1,40);
    for k=1:2
        ki = codvec(k) - 1;
        ss = sign_arr(codvec(k));
        track = mod(ki,5);
        index = fix(ki/5);

        if(track == 0)
            track = 1;
            index = index * 2^6;
        elseif(track == 1)
            if(k == 1)
                track = 0;
                index = index * 2;
            else
                track = 1;
                index = (index * 2^6) + 16;
            end
        elseif (track == 2)
            track = 1;
            index = index * 2^6 + 32;
        elseif (track == 3)
            track = 0;
            index = index * 2 + 1;
        elseif (track ==4)
            track = 1;
            index = index *2^6 + 48;
        end

        if(ss >0)
            cod(ki+1) = 0.9998779296875;
            rsign = rsign + 2^track;
        else
            cod(ki+1) = -1;
        end

        f_ind = f_ind + index;

    end

    % Filter the Algebraic codevector with Weighted Synthesis Speech
    % filters impulse response
    y = conv(h1,cod);
    y = y(1:40);
    y2 = y;
    % Quantizing, big job ahead

    ener_code = sum(cod.*cod);
    ener = ener_code * 134217728;

    [s e] =log2(ener); s = fix(s*(2^15));e = e - 1;
    tmp = fix(e*-49320 + ((s * -24660)/2^15)*2) + 2134784; % For MR59 mode
    tmp = tmp * 2^9;
    tmp = tmp + 5571 * qua_gain_code(past_qua_en(1) + 1);
    tmp = tmp + 4751 * qua_gain_code(past_qua_en(2) + 1);
    tmp = tmp + 2785 * qua_gain_code(past_qua_en(3) + 1);
    tmp = tmp + 1556 * qua_gain_code(past_qua_en(4) + 1);
    tmp = fix((fix(tmp/2^15) * 10886)/2^9);

    gcode0_exp = fix(tmp / 2^15);
    gcode0_fra = fix(tmp - gcode0_exp * 32768);


    % y1, y2, target, target2

    ener_init = 0.01;
    coeff=zeros(1,5);
    coeff(1) = gcoeff(1);
    coeff(2) = -2.0 * gcoeff(2);
    coeff(3) = sum(y2.*y2) + ener_init;
    coeff(4) = sum(target.*y2)*-2.0 + ener_init;
    coeff(5) = sum(y1.*y2)*2.0 + ener_init;

    gcode0 = round(2^(gcode0_exp + gcode0_fra/32767));

    table_len = 64;
    dist_min = Inf;
    for cd = 1:table_len
        g_pitch = table_gain_lowbitrates(cd,1);
        g_code = table_gain_lowbitrates(cd,2);
        if(g_pitch < gp_limit)
            g_code = g_code * gcode0;
            g2_pitch = g_pitch^2;
            g2_code = g_code * g_code;
            g_pit_cod = g_code * g_pitch;
            tmp = coeff(1) * g2_pitch;
            tmp  = tmp + coeff(2) * g_pitch;
            tmp = tmp + coeff(3) * g2_code;
            tmp = tmp + coeff(4) * g_code;
            tmp = tmp + coeff(5) * g_pit_cod;

            if(tmp < dist_min)
                dist_min = tmp;
                t_index = cd;
            end
        end
    end

    out_bytes(5 + (index_subframe -1)*4) = f_ind;
    out_bytes(6 + (index_subframe -1)*4) = rsign;
    out_bytes(7 + (index_subframe -1)*4) = t_index - 1;

    g_pitch = table_gain_lowbitrates(t_index,1);
    g_code_tmp = fix(table_gain_lowbitrates(t_index, 2) * 4096);
    gcode_0 = round(2^14 * 2^(gcode0_fra/32767));
    if(gcode0_exp < 11)
        gain_cod = fix((g_code_tmp * gcode_0) / (2^(25 - gcode0_exp)));
    else
        t = fix((g_code_tmp * gcode_0) * (2^(gcode0_exp - 9)));
        if(fix(t / 2^(gcode0_exp - 9)) ~= ( g_code_tmp * gcode_0))
            gain_cod = 32767;
        else
            gain_cod = fix(t / 2^16);
        end
    end


    gain_cod = gain_cod * 0.5;

    for tt = 4:-1:2
        past_qua_en(tt) = past_qua_en(tt-1);
    end
    past_qua_en(1) = 160 + t_index - 1;
    % Subframe post processing

    final_exc(160 + index_subframe*40 - 39:160 + index_subframe*40) = floor((g_pitch * final_exc(160 + index_subframe*40 - 39:160 + index_subframe*40) ...
        + gain_cod * cod)+ 0.5); % to be completed

    par_in = final_exc(160 + index_subframe*40 - 39:160 + index_subframe*40);
    mem_synth = flipud(synth(31:40)')';
    %mem_synth=[100:-10:10];
    z=filtic(1,A_quant_subframes_lpc(index_subframe,:),(mem_synth)); % Check if you need flipud to mem_err
    [synth t_zf] = filter(1,A_quant_subframes_lpc(index_subframe,:),par_in,z);

    mem_err =  work_absframe(31:40) - synth(31:40);%[100:-10:10];
    mem_err = flipud(mem_err')';
    mem_w0 = target(31:40) - y1(31:40)*g_pitch - y2(31:40)*gain_cod;
    mem_w0 = flipud(mem_w0')';

end

final_exc(1:160)=final_exc(161:320);
final_exc(161:320) = zeros(1,160);

return