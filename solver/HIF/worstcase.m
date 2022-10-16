% run("/home/jyliu/GenerateMatrices/meshpart/meshpart_startup.m")

nvec = [16];
minvtx = 64;
method = "metis";
tol = 1e-3;
demoHIF = 0;
for i = 1:length(nvec)
    n = nvec(i);
    [A,xy] = grid3d(n);
    % Basic settings.
    Axy.A = A;
    if exist("xy")
        Axy.xy = xy;
    else
        Axy.xy = [];
    end
    HIF = HIFGraph(Axy,minvtx,method);
    HIF = Factorization(HIF,tol,demoHIF);
    
    % Worst case.
    AHIFinvAop = @(x,str) (str == "notransp") * HIFSolve(HIF,A*x) + (str == "transp") * A*HIFSolve(HIF,x);     
%     AHIFinvA = A;
%     AHIFinvA = HIFSolve(HIF, AHIFinvA);
%     B = AHIFinvA - speye(size(AHIFinvA));
%     [u,s,v] = svds(B,1);
    [u,s,v] = svds(AHIFinvAop, [size(A,1),size(A,2)], 1);
    x = v;
    b = A*x;     
    xsol = HIFSolve(HIF,b);
    disp(" Worst case")
    disp(" Relative error:")
    disp(norm(xsol - x)/norm(x))
    
    % All one case.
    x = ones(size(A,2),1);
    b = A*x;
    xsol = HIFSolve(HIF,b);
    disp(" All one case")
    disp(" Relative error:")
    disp(norm(xsol - x)/norm(x))

end