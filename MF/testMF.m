%%
% 1 - 5
% |   |
% 2 - 6
% |   |
% 3 - 7
% |   |
% 4 - 8

A1 = [4,-1,0,0;
    -1,4,-1,0;
    0,-1,4,-1;
    0,0,-1,4];
A2 = [0,-1;
    -1,0];
I2 = eye(2);
I4 = eye(4);
A = kron(I2,A1) + kron(A2,I4);
clear A1 A2 I2 I4;
xy = [1,4;
    1,3;
    1,2;
    1,1;
    2,4;
    2,3;
    2,2;
    2,1];
Axy.A = A;
Axy.xy = xy;

MF = MFGraph(Axy);
MF = BuildTree(MF);
MF = SetNbNodes(MF);
% DemoPart(testHIF)
MF = FillTree(MF);
MF = Factorization(MF);
x = ones(size(A,1),1);
b = A*x;
MF = MFSolve(MF,b);
disp('Relative error:')
disp(norm(MF.solution - x)/norm(x))

%%
[A,xy] = grid5(16);
A = full(A);
Axy.A = A;
Axy.xy = xy;
% Axy.xy = []

MF = MFGraph(Axy);
MF = BuildTree(MF);
MF = SetNbNodes(MF);
DemoPart(MF)
MF = FillTree(MF);
MF = Factorization(MF);
x = rand(size(A,1),1);
b = A*x;
MF = MFSolve(MF,b);
disp('Relative error:')
disp(norm(MF.solution - x)/norm(x))
