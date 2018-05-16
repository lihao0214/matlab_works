function [unmixedsig,W,WQ]=INFOMAXICA(mixedsig,miu,error,Maxdiedaicishu)
%???mixedsig??????miu?????error?W????Maxdiedaicishu???????
%???W??????unmixedsig????????WQ?????
[mixedsig_white,Q]=whitening(mixedsig);%??whitening??
% X = mixedsig_white;                     %?????????
X=mixedsig;
Q=eye(2);
daxiao=size(X);                         %???????
W=rand(daxiao(1),daxiao(1));              %????
W(1,:)=W(1,:)/norm(W(1,:));
W(2,:)=W(2,:)/norm(W(2,:));
%????
maxdiedaicishu=Maxdiedaicishu;           %??????
diedaicishu=0;                %???????
Wold=W;
tmp=[];
for i=1:maxdiedaicishu
    diedaicishu=diedaicishu+1;
    if  diedaicishu>maxdiedaicishu
        break;
    end
    y=Wold*X;
    gy=1-2./(1+exp(-y));
    % gy=-2*tanh(y);
    % gy=-sign(y);
    gyy=gy*y';
    t=eye(size(gyy))+gyy/daxiao(2);
    Wzeng=miu*t*Wold;
    W=W+Wzeng;
    if sum(sum(abs(Wold-W)))<error
        Wold=W;
        break;
    end
    Wold=W;
    tmp=[tmp;norm(t,'fro')];
end
W(1,:)=Wold(1,:)/norm(Wold(1,:),2);
W(2,:)=Wold(2,:)/norm(Wold(2,:),2);
unmixedsig=W*X;                %????
WQ=W*Q;                             %????
