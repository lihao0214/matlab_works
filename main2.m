clc;
close all;
clear all;
fid=fopen('speaker_signal.pcm'); far=fread(fid,'short'); fclose(fid);
fid=fopen('mic_signal.pcm'); mic=fread(fid,'short'); fclose(fid);
far = resample(far,2,1);
mic = resample(mic,2,1);

%fid=fopen('QA1far.pcm'); far=fread(fid,'short'); fclose(fid);
%fid=fopen('QA1near.pcm'); mic=fread(fid,'short'); fclose(fid);

start=0 * 8000+1;
stop= 30 * 8000;
microphone=mic(start:stop);
TheFarEnd=far(start:stop);
avtime=1;

% 16000 to make it compatible with the C-version
[emicrophone,tdel]=compsup(microphone,TheFarEnd,avtime,16000);

spclab(8000,TheFarEnd,microphone,emicrophone);


