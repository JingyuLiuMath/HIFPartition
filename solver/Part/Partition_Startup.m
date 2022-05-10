function Partition_Startup()
% Partition_Startup  Startup file for partitions.

% meshpart
run("..\meshpart\meshpart_startup.m")
% metis
run("..\metismex\METIS_startup.m")

end