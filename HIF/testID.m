m = 6;
[A,xy] = grid5(m);
n = m^2;
n = n/2;
A = A(1:n,n+1:end);
A = full(A);
[Ahat,T,p,k] = ID(A);
norm(A(:,p) - [Ahat, Ahat*T])