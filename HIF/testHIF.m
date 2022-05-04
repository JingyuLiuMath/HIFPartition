%% 2D example
[A,xy] = grid5(16);
A = full(A);
Axy.A = A;
Axy.xy = xy;

method = "Specpart";
% method = "Geopart";

HIF = HIFGraph(Axy);
HIF = BuildTree(HIF,method);
HIF = SetNbNode(HIF);
% DemoPart(HIF)
DemoFinalPart(HIF);
HIF = FillTree(HIF);
HIF = Factorization(HIF,1e-3,0);

x = ones(size(A,1),1);
b = A*x;
HIF = HIFSolve(HIF,b);
disp('Relative error:')
disp(norm(HIF.solution - x)/norm(x))

%% 3D example.
[A,xy] = grid3d(8);
A = full(A);
Axy.A = A;
Axy.xy = xy;

method = "Specpart";
% method = "Geopart";

HIF = HIFGraph(Axy);
HIF = BuildTree(HIF,method);
HIF = SetNbNode(HIF);
% DemoPart(HIF)
% DemoFinalPart(HIF);
HIF = FillTree(HIF);
HIF = Factorization(HIF,1e-3,1);

x = rand(size(A,1),1);
b = A*x;
HIF = HIFSolve(HIF,b);
disp('Relative error:')
disp(norm(HIF.solution - x)/norm(x))

%% Triangular example
[A,xy] = gridt(32);
A = full(A);
Axy.A = A;
Axy.xy = xy;

method = "Specpart";
% method = "Geopart";

HIF = HIFGraph(Axy);
HIF = BuildTree(HIF,method);
HIF = SetNbNode(HIF);
% DemoPart(HIF)
DemoFinalPart(HIF);
HIF = FillTree(HIF);
HIF = Factorization(HIF,1e-3,0);

x = rand(size(A,1),1);
b = A*x;
HIF = HIFSolve(HIF,b);
disp('Relative error:')
disp(norm(HIF.solution - x)/norm(x))

%% Solve Ainv
Ainv = zeros(size(A));
for i = 1:size(A,1)
    ei = zeros(size(A,1),1);
    ei(i) = 1;
    HIF = HIFSolve(HIF,ei);
    Ainv(:,i) = HIF.solution;
end
dA = Ainv - inv(A);
disp('Norm(dA):')
disp(norm(dA))
