function [icasig,W]=FASTICA(mixedsig)
% 输入mixedsig是混合信号
% 输出icasig 是恢复信号， W 是解混矩阵
mixedmean=mean(mixedsig')';                       %混合信号的均值
[x, whiteningmatrix ]=whitening(mixedsig)             %白化处理
[vectorsize,numsamples]=size(x);                     %矩阵的大小
b=zeros(vectorsize);
numofic=vectorsize
for r=1:numofic
    i=1;
    maxnumiterations=5000;                         %最大迭代次数
    w=[1;0];                        %迭代矩阵初值
    w=w/norm(w);
    while i<=maxnumiterations+1
        w=(x*((x'*w).^3))/numsamples-3*w;
        w=w/norm(w);
        w=w-b*b'*w;
        w=w/norm(w);
        i=i+1;
    end
W(r,:)=w'*whiteningmatrix;                            %解混矩阵
b(:,r)=w;
   end
icasig=W*mixedsig+(W*mixedmean)*ones(1,numsamples);  %恢复源信号
