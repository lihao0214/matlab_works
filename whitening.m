function  [whitesig,whiteningmatrix]=whitening(mixedsig)
% 输入mixedsig 是混合信号
% 输出whitsig白化信号，whiteningmatrix是白化矩阵
omixedsig=zeros(size(mixedsig));             %产生初值
mixedmean=mean(mixedsig')';                  %混合信号的均值
omixedsig=mixedsig-mixedmean*ones(1,length(mixedsig)); %去均值处理
covariancematrix=cov(omixedsig',1);
[e,d]=eig(covariancematrix);
eigenvalues=flipud(sort(diag(d)));
whiteningmatrix=inv(sqrt(d)+0.01*eye(2))*e';             %求得白化矩阵
dewhiteningmatrix=e*sqrt(d);
whitesig=whiteningmatrix*mixedsig;            %白化处理
