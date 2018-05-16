%FASTICA算法分离混合语音
clear
load('Speech4.mat');         %加载语音
s1=Speech4(1,:);                    %第一个语音
s2=Speech4(3,:);                    %第二个语音
s1=s1/max(abs(s1));
s2=s2/max(abs(s2));
figure                              %第一个图形
subplot(121)                        %建立子图
plot(s1)                            %第一个语音波形图
xlabel('样点数')
ylabel('幅度')
subplot(122)                        %建立子图
plot(s2)                            %第二个语音波形图
xlabel('样点数')
ylabel('幅度')
s=[s1;s2];                          %将两个源信号放到一个矩阵中
aorig=rand(size(s,1));              %产生混合矩阵
mixedsig=aorig*s;                   %混合观测信号
ms1=mixedsig(1,:);
ms2=mixedsig(2,:);
figure                             %第二个图形
subplot(121)
plot(ms1)                           %第一个混合语音信号波形图
xlabel('样点数')
ylabel('幅度')
subplot(122)
plot(ms2)                           %第二个混合语音信号波形图
xlabel('样点数')
ylabel('幅度')
icasig=FASTICA(mixedsig)            %调用FASTICA函数
is1=icasig(1,:);
is2=icasig(2,:);
figure                             %第三个图形
subplot(121)
plot(is1)                          %第一个分离语音信号波形图
xlabel('样点数')
ylabel('幅度')
subplot(122)
plot(is2)                          %第二个分离语音信号波形图
xlabel('样点数')
ylabel('幅度')
