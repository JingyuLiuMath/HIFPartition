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
tol = 1e-3;
demoHIF = 0;

%% HIF process.
% profile on
HIF = HIFGraph(Axy,minvtx,method);
% DemoPart(HIF)
% DemoFinalPart(HIF);
HIF = Factorization(HIF,tol,demoHIF);
% profile viewer
% profsave(profile('info'),'profile_HIF')

%% Solve linear systems.
x = ones(size(A,1),1);
b = A*x;
xsol = HIFSolve(HIF,b);
disp(" Relative error:")
disp(norm(xsol - x)/norm(x))

%% Solve Ainv.
I = eye(size(A,1));
Ainv = HIFSolve(HIF,I);
dA = Ainv - inv(A);
disp(" norm(dA) / Norm(A):")
disp(norm(dA,"Inf") / norm(A, "Inf"))

 imagesc(abs(dA))

 % 2986: 10, 186, 0
 % 1706: 10, 106, 0