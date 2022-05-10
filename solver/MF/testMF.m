%% 2D example.
[A,xy] = grid5(32);

%% 3D example.
[A,xy] = grid3d(16);

%% Triangular example.
[A,xy] = gridt(32);

%% Partition Methods.
A = full(A);
Axy.A = A;
if exist("xy")
    Axy.xy = xy;
else
    Axy.xy = [];
end
method = "Specpart";
% method = "Geopart";

%% MF process.
MF = MFGraph(Axy);
MF = BuildTree(MF,method);
MF = SetNbNode(MF);
% DemoPart(MF)
% DemoFinalPart(MF);
MF = FillTree(MF);
MF = Factorization(MF,0);

%% Solve linear systems.
x = rand(size(A,1),1);
b = A*x;
MF = MFSolve(MF,b);
disp(" Relative error:")
disp(norm(MF.solution - x)/norm(x))

%% Solve Ainv.
I = eye(size(A,1));
HIF = HIFSolve(HIF,I);
Ainv = HIF.solution;
dA = Ainv - inv(A);
disp(" Norm(dA):")
disp(norm(dA,"Inf"))