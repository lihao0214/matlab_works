clc;
clear;
close all;
%%
g=@(w,z) ((w'*z).*exp(-(w'*z).^2./2));
f=@(w,z) (1-(w'*z).^2).*exp(-(w'*z).^2./2);
r=@(w,z,L) ((w'*z).*g(w,z)*ones(L,1))./L;
G=@(W,z) (W'*z).^3;
F=@(W,z) 3*(W'*z).^2;
%G=@(W,z) (W'*z).*exp(-(W'*z).^2/2);
%F=@(W,z) (1-(W'*z).^2).*exp(-(W'*z).^2/2);
%%
load('Speech4.mat');
s1=Speech4(1,:);
s2=Speech4(4,:);
s3=Speech4(2,:);
s1=s1/max(abs(s1));
s2=s2/max(abs(s2));
s3=s3/max(abs(s3));
s=[[s1;s2],[s1;s2]];
[nic,L]=size(s);
aorig=rand(nic);
%aorig=eye(nic);
x=aorig*s;%A*[s1;s2]
u=mean(x,2);
x=x-repmat(u,1,L);
p=ica_signal(x,2,1);
cor=x*x'/L;
[v,DD]=eig(cor);
W=diag(max(diag(DD),eps).^-0.5)*v';%[d1,0;0,d2]*[v1';v2']
z=W*x;
Phat=z*z.'/L;
d=eye(nic);
pv=zeros(nic);
for m=1:1:nic,
    w=d(:,m);
    for k=1:1:32,
        w=cor*w;
        w=w-pv*(pv'*w);
        w=w/max(eps,norm(w));
    end
    pv(:,m)=w;
end
C=cor*pv;
D=diag((sum(abs(C),1)./max(sum(abs(pv),1),eps)).^-0.5);
Rw=D*pv';%D^-0.5*[v1';v2']
z=Rw*x;
% z=pv'*x;
Iters=64;
Errs=zeros(Iters,1);
b=eye(nic);
wold=zeros(nic);
for k=1:1:Iters,
    %     for n=1:1:nic,
    %         w=b(:,n);
    %         p=z*g(w,z)'/L;
    %         q=f(w,z)*ones(L,1)/L;
    %         w=p-q*w;
    %         b(:,n)=w;
    %     end
    p=z*G(b,z)'/L;%E[z*g(w'*z)]-E[g'(w'*z)]*w,W'=[w1,w2]'=[w1';w2']*z=[w1'*z;w2'*z]'=z*[w1'*z;w2'*z]'
    q=F(b,z)*ones(L,1)/L;%W*diag([w1'*z;w2'*z]*ones(L,1)/L)
    b=p-b*diag(q);%[v1,v2]*[d1,0;0,d2]=[d1*v1,d2*v2]
    [v,D]=eig(b'*b);%eig(W'*W)=v*D*v',rW'*(b'*b)*rW=I
    rW=v*diag(max(diag(D),eps).^-0.5)*v';%(b'*b)^-1/2*(b'*b)*(b'*b)^-1/2=I
    b=b*rW;
    Errs(k,1)=norm(abs(b'*wold)-eye(nic),'fro');
    wold=b;
end
W=b';

% b=zeros(nic);%(R=WV)A=PD,R=P*D*inv(A),inv(R)=A*D^-1*P=[h1,h2]*[1/d1,0;0,1/d2]*P=[h1/d1,h2/d2]*P,y=P*D*[s1;s2],inv(R)*(y=P*[d1*s1;d2*s2])=A*s
% for m=1:1:nic,
%     wold=zeros(nic,1);
%     w=d(:,m);
%     for k=1:1:32,
%         p=z*g(w,z).'/L-r(w,z,L)*w;
%         q=f(w,z)*ones(L,1)/L-r(w,z,L);
%         %p=(w'*z).^3;
%         %w=z*p'/L-3*w;
%         %p=z*g(w,z).'/L;
%         %q=f(w,z)*ones(L,1)/L;
%         %w=p-q*w;
%         w=w-p/max(eps,q);
%         w=w-b*(b'*w);
%         w=w/max(eps,norm(w));
%     end
%     b(:,m)=w;
% end
% W=b';
V=W*Rw;
P=inv(V);

% for m=1:1:nic,
%     %     wold=zeros(nic,1);
%     %     w=d(:,m);
%     %     w=w/norm(w);
%     %     for k=1:1:64,
%     %         PhatW=Phat*w;
%     %         w=(z*((z'*w).^3))/L-2*w-conj(PhatW)*(w.'*PhatW);
%     %         w=w-b*b'*w;
%     %         w=w/norm(w);
%     %         error=norm(w-wold);
%     %         wold=w;
%     %         dbg=[dbg;error];
%     %     end
%     %     b(:,m)=w;
%     W=eye(2);
%     for k=1:1:32,
%         W=(z*(z'*W).^3)/L-3*W;
%         [v,Q]=eig(W'*W);
%         W=W*v*Q^-0.5*v';
%         %         for n=1:1:16,
%         %             W=W/norm(W,2);
%         %             W=1.5*W-0.5*W*W'*W;
%         %         end
%     end
% end
%a=b'*z;
a=W*z;
u=P*a;
p1=P(:,1)*a(1,:);
p2=P(:,2)*a(2,:);
plot(x(1,:),x(2,:),'.green');
hold on;
plot(z(1,:),z(2,:),'.blue');
hold on;
plot(a(1,:),a(2,:),'.red');
hold on;
plot(p1(1,:),p2(1,:),'.black');
