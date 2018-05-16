function [W] = weight_transform(w, K, R, Ms)
if 0,
    % Ms = Ms*2; % resolutions
    Ms = Ms*4;
    G = fft(w(:,1:1+K/2),Ms,1);
    X = zeros(Ms/2,K/2-1);

    for m = 1:2:K/2-1, % odd bands
        x = reshape(G(:,m+1),Ms/4,4);
        x = [x(:,2); x(:,3)];
        X(:,m) = x;
    end

    for m = 2:2:K/2-1, % even bands
        x = reshape(G(:,m+1),Ms/4,4);
        x = [x(:,4); x(:,1)];
        X(:,m) = x;
    end
    X = X(:);

    m = 0;
    x = reshape(G(:,m+1),Ms/4,4); X = [x(:,1); X];

    m = K/2;
    if mod(m,2)==0                    % For N/2 is even, eg. N = 8, 16 subbands
        x = reshape(G(:,m+1),Ms/4,4); X = [X; x(:,4)];
    else                              % For N/2 is odd, eg. N = 10 subbands
        x = reshape(G(:,m+1),Ms/4,4); X = [X; x(:,2)];
    end

    X = [X; 0; conj(X(end:-1:2))];
    x = ifft(X);
    W = x(1:Ms*R/4);

else
    Ms=2*Ms;%double-resolution
    G=fft(w(:,1:1+K/2),Ms,1);
    X=zeros(Ms/4,K/2-1);

    for m=4:4:K/2-1,
        x=reshape(G(:,1+m),Ms/8,8);
        x=[x(:,8);x(:,1)];
        X(:,m)=x;
    end

    for m=1:4:K/2-1,
        x=reshape(G(:,1+m),Ms/8,8);
        x=[x(:,2);x(:,3)];
        X(:,m)=x;
    end

    for m=2:4:K/2-1,
        x=reshape(G(:,1+m),Ms/8,8);
        x=[x(:,4);x(:,5)];
        X(:,m)=x;
    end

    for m=3:4:K/2-1,
        x=reshape(G(:,1+m),Ms/8,8);
        x=[x(:,6);x(:,7)];
        X(:,m)=x;
    end

    X=X(:);
    x=reshape(G(:,1),Ms/8,8);
    X=[x(:,1);X];
    x=reshape(G(:,1+K/2),Ms/8,8);
    switch mod(K/2,4),
        case 0
            x=x(:,8);
        case 1
            x=x(:,2);
        case 2
            x=x(:,4);
        case 3
            x=x(:,6);
    end
    X=[X;x];
    X=[X;0;conj(X(end:-1:2))];
    x=ifft(X);
    Ms=Ms/2;
    L=Ms*R;
    W=x(1:L);
end

end
