%% 2D example.
[A,xy] = grid5(32);

%% 3D example.
[A,xy] = grid3d(16);

%% Triangular example.
[A,xy] = gridt(32);

%% Basic settings.
A = full(A);
Axy.A = A;
if exist("xy")
    Axy.xy = xy;
else
    Axy.xy = [];
end
method = "metis";
% method = "meshpart_specpart";
% method = "meshpart_geopart";
demoMF = 0;

%% MF process.
% profile on
MF = MFGraph(Axy);
MF = BuildTree(MF,method);
MF = SetNbNode(MF);
% DemoPart(MF)
% DemoFinalPart(MF);
MF = FillTree(MF);
MF = Factorization(MF,demoMF);
% profile viewer
% profsave(profile('info'),'profile_MF')

%% Solve linear systems.
x = rand(size(A,1),1);
b = A*x;
MF = MFSolve(MF,b);
disp(" Relative error:")
disp(norm(MF.solution - x)/norm(x))

%% Solve Ainv.
I = eye(size(A,1));
MF = MFSolve(MF,I);
Ainv = MF.solution;
dA = Ainv - inv(A);
disp(" Norm(dA):")
disp(norm(dA,"Inf"))
