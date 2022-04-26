% Recording some aborted code.
%% ID decomposition.
% [~,R,p] = qr(A,0);
% (m * n) = (m * n) (n * n)
% sizeR = size(R,1);
% k = 0;
% while k + 1 <= sizeR && R(k+1,k+1) ~=0
%     k = k + 1;
% end
% 
% R1 = R(1:k,1:k);
% R2 = R(1:k,k+1:end);
% Ahat = Q(:,1:k)*R1;
% Ahat = A(:,p(1:k));
% T = R1\R2;
% 
% p1 = p(1:k);
% p2 = p(k+1:end);
% [p1,Ip1] = sort(p1);
% [p2,Ip2] = sort(p2);
% Ahat = Ahat(:,Ip1);
% T = T(Ip1,Ip2);
% p = [p1,p2];