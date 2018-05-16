% AMR-NB Encoder
% Started on 2nd October, 2009
% Ravi Lakkundi

% Encodes only in MR59 mode
% Can be decoded by standard AMR NB Decoder provied by 3gpp
% The dcoder exe is attached with this file or can be downloaded
% from mathworks

% Limitations:
% 1. Supports Only MR59
% 2. Supports only basic encoding, Not supporting VAD, SID, DTX etc
% 3. Fractional pitch not supported


% If you want to hire me please send a mail at
% ravilakkundi@gmail.com
% For my profile please log-on to linkedin

clear all
clc


lsp_old=[ 0.9595
    0.8413
    0.6549
    0.4154
    0.1423
    -0.1423
    -0.4154
    -0.6549
    -0.8413
    -0.9595];

load order_MR59.mat
lsp_old= lsp_old';
lsp_old_q = lsp_old;

format long e
speech = wavread('speaker.wav');
fp = fopen('test.out','w');
NumOfFrames = floor(length(speech)/160);

% orig_speech is a matrix of frames of speech with 160 cols and NumOfFrames
% Rows

orig_speech = reshape(speech(1:160*NumOfFrames), 160, NumOfFrames); % Will help to plot
orig_speech = orig_speech';

% Make the file multiple of 160 samples
% Ignore the last frame no problem

Encoded_bytes = zeros(1,19);

% Some initial states
Final_St_80Hz=[0 0];
past_rq = zeros(1,10);

global Final_St_80Hz
global lsp_old;
global lsp_old_q;
global past_rq;
global mem_w_prebig;
global past_in;
global past_out;
mem_err = zeros(1,10);
mem_w0 = mem_err;
global past_frame_160_sample;
global mem_err;
global mem_w0;
global final_exc;
global h1;
global past_qua_en;
global qua_ener_index;
global synth;
save;



global window_200_40;
working_frame = [];
global working_frame;
working_frame = zeros(1,240);
mem_w_prebig=zeros(1,10);
past_in = zeros(1,10);
past_out = zeros(1,10);
past_frame_160_sample = zeros(1,160);
final_exc = zeros(1,320); % Two past excitation frames
h1 = zeros(1,320);
past_qua_en = [783 783 783 783];
synth = zeros(1,40);
Frame_Indx = 1;
qua_ener_index = 160;
fwrite(fp,[35 33 65 77 82 10],'uint8'); % #AMR!
while NumOfFrames ~= 0
    Encoded_bytes = encode_amrnb(orig_speech(Frame_Indx, :).*32767, Encoded_bytes);

    %     Bit ordering according to Std
    %     Table B.3: Ordering of the speech encoder bits for the 5.9 kbit/s mode: table2(j)
    %      0	   1	   4	   5	   3	   6	   7	   2	  13	  15
    %      8	   9	  11	  12	  14	  10	  16	  28	  74	  29
    %     75	  27	  73	  26	  72	  30	  76	  51	  97	  50
    %     71	  96	 117	  31	  77	  52	  98	  49	  70	  95
    %    116	  53	  99	  32	  78	  33	  79	  48	  69	  94
    %    115	  47	  68	  93	 114	  46	  67	  92	 113	  19
    %     21	  23	  22	  18	  17	  20	  24	 111	  43	  89
    %    110	  64	  65	  44	  90	  25	  45	  66	  91	 112
    %     54	 100	  40	  61	  86	 107	  39	  60	  85	 106
    %     36	  57	  82	 103	  35	  56	  81	 102	  34	  55
    %     80	 101	  42	  63	  88	 109	  41	  62	  87	 108
    %     38	  59	  84	 105	  37	  58	  83	 104

    temp = 0;
    fin_out = zeros(1,14);
    k = 1;
    ft = 1;
    for p = 1:118
        if(bitand(Encoded_bytes(order_MR59(k) + 1) , order_MR59(k+1)))
            temp = temp + 1;
        end
        if(mod(p,8))
            temp = temp *2;
        else
            fin_out(ft) = temp;
            temp = 0;
            ft = ft + 1;
        end
        k = k + 2;
    end

    fin_out = [20 fin_out temp*2]; % Mode appended
    fwrite(fp,fin_out,'uint8');
    NumOfFrames = NumOfFrames - 1;
    Frame_Indx = Frame_Indx + 1;
    disp('...');
end

fclose(fp);

disp('Iam done finally');