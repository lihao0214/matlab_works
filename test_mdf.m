%%
clc;
close all;
clear all;
addpath('./aec_record');
fid = fopen('echo004.raw', 'rb');
d = fread(fid, 'int16');
d = d./32768;
fclose(fid);
sid = fopen('ref004.raw', 'rb');
u = fread(sid, 'int16');
u = u./32768;
fclose(sid);
u = u(1:floor(length(u)/64)*64);
d = d(1:floor(length(d)/64)*64);

%%

speex_mdf_out = speex_mdf(8000, u, d, 1024, 64);