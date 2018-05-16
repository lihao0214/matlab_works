% function [Y] = radix_2_2(x, N),
x = floor(128 * cos(2*pi*1023*[0:1:255]./8000));
N = 256;
stage = ceil((log2(N))/2);
group = 1;
n_tilt = floor(N/4);

for m = 0:1:stage-1,
    t = [];
    for n = 1:1:group,
        a = 1 + (n-1)*(4*n_tilt) + 0 * n_tilt;
        b = 1 + (n-1)*(4*n_tilt) + 1 * n_tilt;
        c = 1 + (n-1)*(4*n_tilt) + 2 * n_tilt;
        d = 1 + (n-1)*(4*n_tilt) + 3 * n_tilt;
        for k = 0:1:n_tilt-1,
            t = [t; x(a+k); x(b+k); x(c+k); x(d+k)];
            s1 = x(a+k) + x(c+k);
            s2 = x(a+k) - x(c+k);
            x(a+k) = s1;
            x(c+k) = s2;
            s3 = x(b+k) + x(d+k);
            s4 = x(b+k) - x(d+k);
            x(b+k) = s3;
            x(d+k) = s4 * (-j);
            s1 = x(a+k) + x(b+k);
            s2 = x(a+k) - x(b+k);
            x(a+k) = s1;
            x(b+k) = s2;
            s3 = x(c+k) + x(d+k);
            s4 = x(c+k) - x(d+k);
            x(c+k) = s3;
            x(d+k) = s4;
            x(b+k) = x(b+k) * exp(-j*2*pi*k*2/(4*n_tilt));
            x(c+k) = x(c+k) * exp(-j*2*pi*k*1/(4*n_tilt));
            x(d+k) = x(d+k) * exp(-j*2*pi*k*3/(4*n_tilt));
        end
    end
    group = group * 4;
    n_tilt = n_tilt/4;
end
index = [0:1:length(x)-1];
ii = bitrevorder(index);
Y = x(1+ii);
