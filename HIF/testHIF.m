%% 2D example.
[A,xy] = grid5(32);
A = full(A);
Axy.A = A;
Axy.xy = xy;

%% 3D example.
[A,xy] = grid3d(16);
A = full(A);
Axy.A = A;
Axy.xy = xy;

%% Triangular example.
[A,xy] = gridt(32);
A = full(A);
Axy.A = A;
Axy.xy = xy;

%% Basic settings.
method = "Specpart";
% method = "Geopart";
tol = 1e-3;
demoHIF = 0;

%% HIF process.
HIF = HIFGraph(Axy);
HIF = BuildTree(HIF,method);
HIF = SetNbNode(HIF);
% DemoPart(HIF)
% DemoFinalPart(HIF);
HIF = FillTree(HIF);
HIF = Factorization(HIF,tol,demoHIF);

%% Solve linear systems.
x = rand(size(A,1),1);
b = A*x;
HIF = HIFSolve(HIF,b);
disp(" Relative error:")
disp(norm(HIF.solution - x)/norm(x))

%% Solve Ainv.
I = eye(size(A,1));
HIF = HIFSolve(HIF,I);
Ainv = HIF.solution;
dA = Ainv - inv(A);
disp(" Norm(dA):")
disp(norm(dA,"Inf"))
