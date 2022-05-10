function [p1,p2,sep1,sep2] = metispart(A)
% METISPART Graph parition using metis. 

[p1,p2,sep] = METIS_SepPartition(A);
% p1 = p1 + sep, p2 = p2, sep1 = sep, sep2 need to be assigned.
p1 = sort([p1,sep]);
sep1 = sep;
sep2 = [];
for i = 1:length(sep1)
    sep1i = sep1(i);
    sep2i = find(A(sep1i,:)~=0);
    sep2i = intersect(sep2i,p2);
    sep2 = sort(unique([sep2,sep2i]));
end

end

