function PT = findpitch(s),
[B,A] = butter(5,700/4000);
s = filter(B,A,s);
R = zeros(143,1);

for k = 1:143,
    R(k) = s(144:223)'*s(144-k:223-k);
end

[R1,T1] = max(R(80:143));
R1 = R1/(norm(s(144-T1:223-T1))+1);
T1 = T1+79;
[R2,T2] = max(R(40:79));
R2 = R2/(norm(s(144-T2:223-T2))+1);
T2 = T2+39;
[R3,T3] = max(R(20:39));
R3 = R3/(norm(s(144-T3:223-T3))+1);
T3 = T3+19;
Top = T1;
Rop = R1;

if R2 >= 0.85 * Rop,
    Rop = R2;
    Top = T2;
end

if R3 > 0.85 * Rop,
    Rop = R3;
    Top = T3;
end

PT = Top;

return;


