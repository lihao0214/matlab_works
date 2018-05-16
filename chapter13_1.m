%均匀分布混合信号的三种基本变换
clear
t=0:0.05:100;                          %产生时间信号
s1=sin(0.5*t);                          %第一个源信号
s2=2*rand(size(s1))-1;                   %第二个源信号
s=[s1;s2];                              %将两个源信号放到一个矩阵中
aorig=rand(size(s,1));                     %产生混合矩阵
mixedsig=aorig*s;                       %混合观测信号
ms1=mixedsig(1,:);
ms2=mixedsig(2,:);
whitesig=whitening(mixedsig);             %白化处理
figure                                 %第一个图形
plot(s1,s2,'k.')                           %源信号的变量密度分布图
axis([-1.5,1.5,-1.5,1.5])                   %处理坐标轴
figure                                 %第二个图形
plot(ms1,ms2,'k.')                      %混合后信号的变量密度分布图
axis([-1.5,1.5,-1.5,1.5])
figure                                %第三个图形
wis1= whitesig(1,:)
wis2= whitesig(2,:)
plot(wis1, wis2,'k.')                     %白化后信号的变量密度分布图
axis([-3,3,-3,3])
% icasig=FASTICA(mixedsig)             %调用FSATICA函数
icasig=INFOMAXICA(mixedsig,0.025,0.0001,1000)
is1=icasig(1,:);
is2=icasig(2,:);
figure                                %第四个图形
plot(is1,is2,'k.')                        % ICA分离后信号的变量密度分布图
axis([-2,2,-2,2])
