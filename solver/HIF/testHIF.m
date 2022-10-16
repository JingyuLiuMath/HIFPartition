%% 2D example.
[A, xy] = grid5(128);

%% 3D example.
[A, xy] = grid3d(32);

%% Triangular example.
[A, xy] = gridt(32);

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
HIFalg = 1;
demoHIF = 0;

%% HIF process.
% profile off
% profile on
HIF = HIFGraph(Axy, minvtx, method);
% DemoPart(HIF)
% DemoFinalPart(HIF);
HIF = Factorization(HIF, tol, HIFalg, demoHIF);
% profile viewer

%% Solve linear systems.
x = randn(size(A, 1), 1);
b = A * x;
xsol = HIFSolve(HIF, b);
disp("Relative error:");
disp(norm(xsol - x) / norm(x));

%% Solve Ainv.
I = eye(size(A, 1));
Ainv = HIFSolve(HIF, I);
dA = Ainv - inv(A);
disp("norm(dA) / Norm(A):");
disp(norm(dA, "Inf") / norm(A, "Inf"));

%% All one vector.
x = ones(size(A, 2), 1);
b = A * x;
xsol = HIFSolve(HIF, b);
disp("Relative error:");
disp(norm(xsol - x) / norm(x));
