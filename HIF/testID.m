%%
m = 6;
[A,xy] = grid5(m);
n = m^2;
n = n/2;
A = A(1:n,n+1:end);
A = full(A);
%% 
A = eye(8)
%% 
[T,p1,p2] = ID(A,0.01);
p = [p1,p2];
disp(norm(A(:,p) - [A(:,p1), A(:,p1)*T]))
