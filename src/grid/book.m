function A = book(k,n)
% BOOK k-page graph
%
%   A = book(k,n) Generate a book of k pages, each with n*2 points.
%   A is the Laplacian; xy is coordinates for a planar drawing.

% Yingzhou Li, 2017
% John Gilbert, 1994.
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

A = grid5(n);
A = blockdiags(A,0,k,k);
blocks = 1:(k*n^2);
for it = 1:k-1
    idx = n^2*it;
    blocks(idx+1 : idx+n) = blocks(1:n);
end
A = contract(A,blocks);

end