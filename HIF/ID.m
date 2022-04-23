function [Ahat,T,p,k] = ID(A)
%ID Compute the ID decomposition of a low-rank matrix A.

[Q,R,p] = qr(A,'vector');

k = 0;
while R(k+1,k+1) ~=0
    k = k+1;
end

R1 = R(:,1:k);
R2 = R(:,k+1:end);
Ahat = Q*R1;
T = R1\R2;

p1 = p(1:k);
p2 = p(k+1:end);
[p1,Ip1] = sort(p1);
[p2,Ip2] = sort(p2);
Ahat = Ahat(:,Ip1);
T = T(Ip1,Ip2);
p = [p1,p2];
% Finally, we have A(:,p) = [Ahat,Ahat*T] and p1, p2 is sorted

end

