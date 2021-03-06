%% 2D example.
[A,xy] = grid5(32);

%% 3D example.
[A,xy] = grid3d(16);

%% Triangular example.
[A,xy] = gridt(32);

%% Basic settings.
Axy.A = A;
if exist("xy")
    Axy.xy = xy;
else
    Axy.xy = [];
end
minvtx = 64;
method = "metis";
% method = "meshpart_specpart";
% method = "meshpart_geopart";
demoMF = 0;

%% MF process.
% profile on
MF = MFGraph(Axy,minvtx,method);
% DemoPart(MF)
% DemoFinalPart(MF);
MF = Factorization(MF,demoMF);
% profile viewer
% profsave(profile('info'),'profile_MF')

%% Solve linear systems.
x = rand(size(A,1),1);
b = A*x;
xsol = MFSolve(MF,b);
disp(" Relative error:")
disp(norm(xsol - x)/norm(x))

%% Solve Ainv.
I = eye(size(A,1));
Ainv = MFSolve(MF,I);
dA = Ainv - inv(A);
disp(" norm(dA) / Norm(A):")
disp(norm(dA,"Inf") / norm(A, "Inf"))
