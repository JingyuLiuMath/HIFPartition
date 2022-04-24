function [T,p1,p2] = ID(A,tol)
%ID Compute the ID decomposition of a low-rank matrix A.

assert(tol < 1, "the tol must be less than 1!")

[p1,p2,T] = klhoID(A,tol);
[p1,Ip1] = sort(p1);
[p2,Ip2] = sort(p2);
T = T(Ip1,Ip2);

% Finally, we have A(:,p2) \approx A(:,p1) * T and p1, p2 is sorted

end

