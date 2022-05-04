%% 2D example
[A,xy] = grid5(16);
A = full(A);
Axy.A = A;
Axy.xy = xy;

method = "Specpart";
% method = "Geopart";

MF = MFGraph(Axy);
MF = BuildTree(MF,method);
MF = SetNbNode(MF);
% DemoPart(MF)
DemoFinalPart(MF);
MF = FillTree(MF);
MF = Factorization(MF,0);

x = rand(size(A,1),1);
b = A*x;
MF = MFSolve(MF,b);
disp('Relative error:')
disp(norm(MF.solution - x)/norm(x))

%% 3D example.
[A,xy] = grid3d(8);
A = full(A);
Axy.A = A;
Axy.xy = xy;

method = "Specpart";
% method = "Geopart";

MF = MFGraph(Axy);
MF = BuildTree(MF,method);
MF = SetNbNode(MF);
% DemoPart(MF)
DemoFinalPart(MF);
MF = FillTree(MF);
MF = Factorization(MF,0);

x = rand(size(A,1),1);
b = A*x;
MF = MFSolve(MF,b);
disp('Relative error:')
disp(norm(MF.solution - x)/norm(x))

%% Triangular example
[A,xy] = gridt(32);
A = full(A);
Axy.A = A;
Axy.xy = xy;
% Axy.xy = [];

method = "Specpart";
% method = "Geopart";

MF = MFGraph(Axy);
MF = BuildTree(MF,method);
MF = SetNbNode(MF);
% DemoPart(MF)
DemoFinalPart(MF);
MF = FillTree(MF);
MF = Factorization(MF,0);

x = rand(size(A,1),1);
b = A*x;
MF = MFSolve(MF,b);
disp('Relative error:')
disp(norm(MF.solution - x)/norm(x))

%% Solve Ainv
Ainv = zeros(size(A));
for i = 1:size(A,1)
    ei = zeros(size(A,1),1);
    ei(i) = 1;
    MF = MFSolve(MF,ei);
    Ainv(:,i) = MF.solution;
end
dA = Ainv - inv(A);
disp('Norm(dA):')
disp(norm(dA))
