function [p1,p2,sep1,sep2] = GraphPart(Axy,method)
% GraphPart Graph partition.

if nargin == 1
    method = "metis";
end

switch method
    case "meshpart_geopart"
        numtries = 30;
        [p1,p2,sep1,sep2] = geopart(Axy.A,Axy.xy,numtries);
    case "meshpart_specpart"
        [p1,p2,sep1,sep2] = specpart(Axy.A);
    case "metis"
        [p1,p2,sep1,sep2] = metispart(Axy.A);
end

end

