function Solver_Startup()
% Solver_Startup  Startup file for HIF.

% meshpart
run("..\meshpart\meshpart_startup.m")
% metis
run("..\metismex\METIS_startup.m")
% part
addpath("./Part")
% MF
addpath("./MF")
% HIF
addpath("./HIF")

end