function Solver_Startup()
% Solver_Startup  Startup file for partitions.

% meshpart
run("..\meshpart\meshpart_startup.m")
% metis
run("..\metismex\METIS_startup.m")
% part
addpath("./Part")

end