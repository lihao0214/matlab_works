function [ERR,outputs,state]=pca_ica(ch1_data,ch2_data,K,R,state),
lambda=0.05;
state.buffer1=[state.buffer1(1+R:end);ch1_data];
state.buffer2=[state.buffer2(1+R:end);ch2_data];
tmpD1=flipud(state.wins.*state.buffer1);
tmpD2=flipud(state.wins.*state.buffer2);
tmpW1(1,1)=sum(tmpD1(1:K:end),1);
tmpW2(1,1)=sum(tmpD2(1:K:end),1);
for k=1:1:K-1,
    tmpW1(1+k,1)=sum(tmpD1(1+K-k:K:end),1);
    tmpW2(1+k,1)=sum(tmpD2(1+K-k:K:end),1);
end

tmpF1=fft(circshift(tmpW1,state.cntr),K);%circshift/n=rR
tmpF2=fft(circshift(tmpW2,state.cntr),K);%circshift/n=rR
pm=zeros(1+K/2,1);

for k=2:1:K/2,
    x=[tmpF1(k,1);tmpF2(k,1)];
    h=state.H(:,:,k);
    z=h*x;
    error=eye(2)-z*z';
    h=h+0.001*error*h;
    state.H(:,:,k)=h;
    tmpF1(k,1)=z(1);
    tmpF2(k,1)=z(2);
    pm(k,1)=norm(error,'fro');
end
ERR=sum(pm,1);

tmpn1=ifft([tmpF1(1:1+K/2);conj(tmpF1(K/2:-1:2))],K);
tmpn2=ifft([tmpF2(1:1+K/2);conj(tmpF2(K/2:-1:2))],K);

state.bufferSys1=[tmpn1,state.bufferSys1(:,1:end-1)];%[new,...,old]
state.bufferSys2=[tmpn2,state.bufferSys2(:,1:end-1)];%[new,...,old]
outputs=zeros(R,2);

for k=0:1:R-1,
    outputs(1+k,1)=K*R*state.bufferSys1(1+mod(state.cntr+k,K),:)*state.wins(1+k:R:end);
    outputs(1+k,2)=K*R*state.bufferSys2(1+mod(state.cntr+k,K),:)*state.wins(1+k:R:end);
end
state.cntr=mod(state.cntr+R,K);%K/bands/bounds
end
