%% Skel 1 
%             mysepx = mysep;
%             for i = 1:length(mysep)
%                 mysepi = mysep(i);
%                 % NEED TO BE CHANGED!
%                 % interact_mysepi = find(obj.root.inputAxy.A(:,mysepi) ~= 0);
%                 index_mysepi = find(obj.sep==mysepi,1);
%                 index_mysepi = find(obj.ANS(:,index_mysepi) ~= 0);
%                 interact_mysepi = obj.nb(index_mysepi);
%                 interact_mysepi = setdiff(interact_mysepi,nodek.sep);
%                 if ~isempty(interact_mysepi)
%                     mysep(i) = -mysepi;
%                 end
%             end
%             s12 = -mysep(mysep < 0);
%             mysep = mysep(mysep > 0);   

%             for i = 1:length(mynb)
%                 mynbi = mynb(i);
%                 % NEED TO BE CHANGED!
%                 % interact_mysepi = find(obj.root.inputAxy.A(:,mysepi) ~= 0);
%                 nodekindex_mynbi = find(nodek.sep==mynbi,1);
%                 index_mynbi = find(nodek.ANS(:,nodekindex_mynbi) ~= 0);
%                 interact_mynbi = nodek.nb(index_mynbi);
%                 interact_mynbi = setdiff(interact_mynbi,obj.sep);
%                 if ~isempty(interact_mynbi)
%                     mynb(i) = -mynbi;
%                 end
%             end
%             s21 = -mynb(mynb < 0);
%             mynb = mynb(mynb > 0);