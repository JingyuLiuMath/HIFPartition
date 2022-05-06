classdef HIFGraph < handle
    % HIFGraph HIF algorithm on general graphs.
    
    properties
        
        % Root properties.
        
        % The following information will be stored only in the root node.
        
        inputAxy; % Input Axy.
        active; % Whether a vertex is eliminated.
        inputVec; % Input vector.
        solution; % Solution.
        demoHIF = 0; % Whether to demo our HIF process.
        
        % Graph properties.
        
        Axy; % Adjacency matrix (and coordinates of vertices).
        vtx; % Vertices on the graph.
        sep; % Separator vertices.
        nb; % Neighbor vertices.
        int; % Interior vertices.
        nbA; % Adjacency matrix of sep (row) and nb (col).
        re; % Redundant sep. We also use check (c) to reprsent it.
        sk; % Skeleton sep. We also use hat (h) to represent it.
        nbre; % Redundant neighbor.
        singlesep = {}; % Sep which only interact with one node.
        complexsep; % Sep - singlesep;
        
        % Partition properties.
        
        numLevels; % Total number of levels.
        level; % Current level, start from 0.
        seqNum; % A node's order in its level.
        endFlag = 0; % Whether the partition ends.
        
        % Tree properties.
        
        parent; % Parent node.
        children = cell(1,2); % Children nodes.
        nbNode = {}; % Neighbor nodes.
        nbNodeSeqNum = []; % Neighbor nodes's seqNum.
        root; % Root node.
        nbInfo = struct([]); % Neighbor nodes' information when skeletonization.
        indexInfo = struct([]); % Index information when merge and split.
        
        % Matrices properties.
        
        % For the following matrices, the fist index is row, and the second
        % index is col.
        
        AII; % Interaction between int and int.
        ASI; % Interaction between sep and int.
        ASS; % Interaction between sep and sep.
        ANS; % Interaction between nb and sep.
        DI; % The D part of LDL factorization about AII.
        LI; % The L part of LDL factorization about AII.
        AIIinvAIS; % AIIinvAIS = AII^{-1} * ASI^{T}.
        
        % Vectors properties.
        
        xI; % The int part of a vector x.
        xS; % The sep part of a vector x.
        
    end
    
    methods
        
        function obj = HIFGraph(Axy,level,seqNum,vtx,sep,nb,nbA)
        % HIFGraph Create a HIF class.
        
        if nargin == 1
            level = 0;
            seqNum = 0;
            n = size(Axy.A,1);
            vtx = 1:1:n;
            sep = [];
            nb = [];
            nbA = [];
            
            active = ones(1,n);
            obj.root = obj;
            obj.inputAxy = Axy;
            obj.active = active;
        end
        
        obj.Axy = Axy;
        obj.level = level;
        obj.seqNum = seqNum;
        obj.vtx = vtx;
        obj.sep = sep;
        obj.nb = nb;
        obj.nbA = nbA;
        
        end
        
        function obj = BuildTree(obj,method)
        % BuildTree Build tree structure according to a graph partition algorithm.
        
        if nargin == 1
            method = "Specpart";
        end
        
        if obj.level == 0
            disp("  ");
            disp(" Start build tree ");
            disp("  ");
        end
        
        numtries = 30; % Only useful when using geopart.
        minvtx = 16; % Don't separate pieces smaller than this.
        
        n = size(obj.Axy.A,1);
        
        if n <= minvtx
            obj.numLevels = obj.level;
            obj.endFlag = 1;
            return
        end
        
        if method == "Geopart"            
            % Geopart.
            [p1,p2,sep1,sep2] = geopart(obj.Axy.A,obj.Axy.xy,numtries);
            sep1 = unique(sep1);
            sep2 = unique(sep2);
            p = {p1,p2};
            partsep = {sep1,sep2};
            childAxy(1).A = obj.Axy.A(p1,p1); childAxy(1).xy = obj.Axy.xy(p1,:);
            childAxy(2).A = obj.Axy.A(p2,p2); childAxy(2).xy = obj.Axy.xy(p2,:);
        else
            % Specpart.
            [p1,p2,sep1,sep2] = specpart(obj.Axy.A);
            sep1 = unique(sep1);
            sep2 = unique(sep2);
            p = {p1,p2};
            partsep = {sep1,sep2};
            childAxy(1).A = obj.Axy.A(p1,p1); childAxy(1).xy = [];
            childAxy(2).A = obj.Axy.A(p2,p2); childAxy(2).xy = [];
        end
        
        
        % Create children HIF.
        for iter = [1,2]
            obj.children{iter} = HIFGraph(childAxy(iter),obj.level+1,obj.seqNum*2+(iter-1),...
                obj.vtx(p{iter}),obj.vtx(partsep{iter}),obj.vtx(partsep{3-iter}),obj.Axy.A(partsep{iter},partsep{3-iter}));
            obj.children{iter}.root = obj.root;
            obj.children{iter}.parent = obj;
        end
        
        % Pass information to its children.
        obj = Pass(obj);
        
        % Clear information.
        obj.Axy = [];
        obj.nbA = [];
        
        % Recursively buildtree.
        for iter = [1,2]
            obj.children{iter} = BuildTree(obj.children{iter});
        end
        
        % Get numLevels from children when partition ends.
        obj.numLevels = max(obj.children{1}.numLevels,obj.children{2}.numLevels);
        
        if obj.level == 0
            disp("  ");
            disp(" End build tree ");
            disp("  ");
        end
        
        end
        
        function obj = Pass(obj)
        % PASS Send parent's sep, nb, nbA to children.
        
        for i = 1:length(obj.sep)
            sepi = obj.sep(i);
            for iter = [1, 2]
                obj_child = obj.children{iter};
                index_sepi_childnode = find(obj_child.vtx == sepi,1);
                if isempty(index_sepi_childnode)
                    continue;
                end
                % Now sepi is a vtx of child.
                index_sepi_childsep = find(obj_child.sep == sepi,1);
                num_childnb = length(obj_child.nb);
                num_childsep = length(obj_child.sep);
                if ~isempty(index_sepi_childsep)
                    % Now sepi is a sep of child, we need to pass nb and nbA.
                    index_addnb_nb = find(obj.nbA(i,:)~=0); % index_addnb_nb is always nonempty!
                    for j = 1: length(index_addnb_nb)
                        addnbj = obj.nb(index_addnb_nb(j)); % addnbj is a neigbor vtx.
                        index_addnbj_childnb = find(obj_child.nb == addnbj,1);
                        % index_addnb_childnb may be nonempty (addnbj has been added).
                        if isempty(index_addnbj_childnb)
                            % Now addnbj is NOT in child's nb,we need to add nb and nbA.
                            obj_child.nb = [obj_child.nb,addnbj];
                            num_childnb = num_childnb+1;
                            obj_child.nbA(index_sepi_childsep,num_childnb) = obj.nbA(i,index_addnb_nb(j));
                        else
                            % Now addnbj is in child's nb, we only need to add nbA.
                            obj_child.nbA(index_sepi_childsep,index_addnbj_childnb) = obj.nbA(i,index_addnb_nb(j));
                        end
                    end
                else
                    % Now sepi is NOT a sep of child, we need to pass sep,nb and nbA.
                    obj_child.sep = [obj_child.sep,sepi];
                    num_childsep = num_childsep+1;
                    index_addnb_nb = find(obj.nbA(i,:)~=0);% index_addnb_nb is always nonempty!
                    for j = 1: length(index_addnb_nb)
                        addnbj = obj.nb(index_addnb_nb(j)); % addnbj is a neighbor vtx.
                        index_addnbj_childnb = find(obj_child.nb == addnbj,1);
                        % index_addnb_childnb may be nonempty (addnbj has been added).
                        if isempty(index_addnbj_childnb)
                            % Now addnbj is NOT in child's nb, we need to add nb and nbA.
                            obj_child.nb = [obj_child.nb,addnbj];
                            num_childnb = num_childnb +1;
                            obj_child.nbA(num_childsep,num_childnb) = obj.nbA(i,index_addnb_nb(j));
                        else
                            % Now addnbj is in child's nb, we only need to add nbA.
                            obj_child.nbA(num_childsep,index_addnbj_childnb) = obj.nbA(i,index_addnb_nb(j));
                        end
                    end
                end
            end
        end
        
        end
        
        function obj = SetNbNode(obj)
        % SetNbNode Set nbNode.
        
        % We stand on the parent level to assign its children's nbNode.
        if obj.endFlag == 1
            return;
        end
        
        for iter = [1,2]
            obj_child = obj.children{iter};
            obj_child.nbNode{end+1} = obj.children{3-iter};
            obj_child.nbNodeSeqNum(end+1) = obj.children{3-iter}.seqNum;
            % We only need to find nbNode from the children node of
            % parent's nbNode.
            if ~isempty(obj.nbNode)
                for i = 1:length(obj.nbNode)
                    nbNodei = obj.nbNode{i};
                    % What we need is to check whether the vtx of nbNodei's
                    % chilren is in the nb of obj_child.
                    for it = [1,2]
                        nbNodei_child = nbNodei.children{it};
                        if isempty(nbNodei_child)
                            continue;
                        end
                        if ~isempty(intersect(obj_child.nb, nbNodei_child.vtx))
                            obj_child.nbNode{end+1} = nbNodei_child;
                            obj_child.nbNodeSeqNum(end+1) = nbNodei_child.seqNum;
                        end
                    end
                end
            end
        end
        
        % Recursively setNbNode.
        for iter = [1,2]
            obj.children{iter} = SetNbNode(obj.children{iter});
        end
        
        end
        
        function obj = FillTree(obj)
        % FillTree Fill tree structure with A.
        
        % We clear the following data: Axy, nbA and sort vtx, sep, nb.
        obj.vtx = sort(obj.vtx);
        obj.sep = sort(obj.sep);
        obj.nb = sort(obj.nb);
        
        if obj.endFlag == 0
            % We only fill the leaf nodes.
            for iter = [1,2]
                obj.children{iter} = FillTree(obj.children{iter});
            end
        else
            % First, we set int = vtx - sep (only holds on leaf nodes).
            obj.int = setdiff(obj.vtx,obj.sep,'sorted');
            % Then, we set the corresponding A**.
            rootHIF = obj.root;
            obj.AII = rootHIF.inputAxy.A(obj.int,obj.int);
            obj.ASI = rootHIF.inputAxy.A(obj.sep,obj.int);
            obj.ASS = rootHIF.inputAxy.A(obj.sep,obj.sep);
            obj.ANS = rootHIF.inputAxy.A(obj.nb,obj.sep);
            
            % Set sep type.
            obj = SetSepType(obj);
        end
        
        end
               
        function obj = SetSepType(obj)
        % SetSepType Set sep type.
       
        ordersep = zeros(length(obj.sep),1);
        
        for k = 1:length(obj.nbNode)
            nodek = obj.nbNode{k};
            [~,tmp,~] = intersect(obj.sep,nodek.nb);
            obj.singlesep{k} = tmp;
            ordersep(tmp) = ordersep(tmp) + 1;
        end
        
        for k = 1:length(obj.nbNode)
            obj.singlesep{k} = obj.sep(intersect(find(ordersep == 1),obj.singlesep{k}));
        end
        
        obj.complexsep = obj.sep(ordersep > 1);
        
        end
        
        function obj = Factorization(obj,tol,demoHIF)
        % Factorization HIF factorization.
        
        if nargin == 1
            disp(" ");
            disp(" Default tol :1e-2 ");
            disp(" ");
            tol = 1e-2;
        end
        
        if nargin > 2
            obj.demoHIF = demoHIF;
        end
        
        disp("  ");
        disp(" Start factorization ");
        disp("  ");
        
        for tmplevel = obj.numLevels:-1:1
            % Sparse elimination.
            obj = RecursiveSparseElim(obj,tmplevel);
            % Skeletonization.
            obj = RecursiveSkel(obj,tmplevel,tol);
            % Demo the process.
            if obj.demoHIF == 1
                DemoHIF(obj,tmplevel);
            end
            % Merge.
            obj = RecursiveMerge(obj,tmplevel-1);
        end
        
        % Root factorization.
        obj = RootFactorization(obj);
        % Demo the process.
        if obj.demoHIF == 1
            DemoHIF(obj,0);
        end
        
        disp("  ");
        disp(" End factorization ");
        disp("  ");
        
        end
        
        function obj = RecursiveSparseElim(obj,whatlevel)
        % RecursiveSparseElim Recursively sparse elimination.
        
        if obj.level == whatlevel
            obj = SparseElim(obj);
        else
            for iter = [1,2]
                if ~isempty(obj.children{iter})
                    obj.children{iter} = RecursiveSparseElim(obj.children{iter},whatlevel);
                end
            end
        end
        
        end
        
        function obj = SparseElim(obj)
        % SparseElim Sparse elimination.
        
        obj.root.active(obj.int) = 0;
        % AII = LI * DI * LI^{T}.
        [obj.LI,obj.DI] = ldl(obj.AII);
        % AIIinvAIS = AII^{-1} * ASI^{T}.
        obj.AIIinvAIS = obj.LI\(obj.ASI');
        obj.AIIinvAIS = obj.DI\obj.AIIinvAIS;
        obj.AIIinvAIS = obj.LI'\obj.AIIinvAIS;
        % ASS = ASS - ASI * AII^{-1} * ASI^{T}.
        obj.ASS = obj.ASS - obj.ASI*obj.AIIinvAIS;
        % ASI = 0;
        
        end
        
        function obj = RecursiveSkel(obj,whatlevel,tol)
        % RecursiveSkel Recursively skeletonization.
        
        if obj.level == whatlevel
            obj = Skel(obj,tol);
            % obj = NoSkel(obj); % No skeletonization.
        else
            for iter = [1,2]
                if ~isempty(obj.children{iter})
                    obj.children{iter} = RecursiveSkel(obj.children{iter},whatlevel,tol);
                end
            end
        end
        
        end
        
        function obj = Skel(obj,tol)
        % Skel Skeletonization.
        
        for k = 1:length(obj.nbNode)
            % We do skel according to nbNode.
            nodek = obj.nbNode{k};
            if nodek.seqNum < obj.seqNum
                obj.nbInfo(k).empty = 1;
                continue;
            end
            
            % The following data are vertices.            
            sep1 = obj.singlesep{k};
            % mysep1C = setdiff(obj.sep,sep1,'sorted');
            mysep1C = [];
            for nok = 1:length(obj.nbNode)
                if nok == k
                    continue;
                else
                    mysep1C = [mysep1C,obj.singlesep{nok}];
                end
            end
            mysep1C = [mysep1C,obj.complexsep];
            nodeksep1C = setdiff(nodek.nb,sep1,'sorted');
                        
            korder = find(nodek.nbNodeSeqNum == obj.seqNum,1);
            sep2 = nodek.singlesep{korder};
            % nodeksep2C = setdiff(nodek.sep,sep2,'sorted');
            nodeksep2C = [];
            for nok = 1:length(nodek.nbNode)
                if nok == korder
                    continue;
                else
                    nodeksep2C = [nodeksep2C,nodek.singlesep{nok}];
                end
            end
            nodeksep2C = [nodeksep2C,nodek.complexsep];
            mysep2C = setdiff(obj.nb,sep2,'sorted'); 
 
            % The following data are indices.
            [~,myindex_sep1] = ismember(sep1,obj.sep);
            [~,nodekindex_sep1] = ismember(sep1,nodek.nb);
            [~,myindex_sep1C] = ismember(mysep1C,obj.sep);
            [~,myindex_mysep2C] = ismember(mysep2C,obj.nb);
            
            [~,nodekindex_sep2] = ismember(sep2,nodek.sep);
            [~,myindex_sep2] = ismember(sep2,obj.nb);
            [~,nodekindex_sep2C] = ismember(nodeksep2C,nodek.sep);
            [~,nodekindex_nodeksep1C] = ismember(nodeksep1C,nodek.nb);
            
            % ID decomposition.
            
            % In the following process, the first "1" or "2" denotes my or
            % nodek, the second "1" or "2" denotes sk or re.
           
            skelmtx1 = [obj.ASS(myindex_sep1C,myindex_sep1);
                obj.ANS(myindex_mysep2C,myindex_sep1)];
            if isempty(skelmtx1)
                obj.nbInfo(k).empty = 1;
                continue
            end
            [T1,p11,p12] = ID(skelmtx1,tol); % skelmtx1(:,p12) = skelmtx1(:,p11) * T1.
            myindex_p11 = myindex_sep1(p11);
            myindex_p12 = myindex_sep1(p12);
            nodekindex_p11 = nodekindex_sep1(p11);
            nodekindex_p12 = nodekindex_sep1(p12);
            obj.nbInfo(k).myindex_p11 = myindex_p11;
            obj.nbInfo(k).myindex_p12 = myindex_p12;
            obj.nbInfo(k).nodekindex_p11 = nodekindex_p11;
            obj.nbInfo(k).nodekindex_p12 = nodekindex_p12;
            obj.re = [obj.re,sep1(p12)];
            obj.nbInfo(k).Th1c1 = T1;
            nodek.nbre = [nodek.nbre,sep1(p12)];
            
            skelmtx2 = [nodek.ASS(nodekindex_sep2C,nodekindex_sep2);
                nodek.ANS(nodekindex_nodeksep1C,nodekindex_sep2)];
            if isempty(skelmtx2)
                obj.nbInfo(k).empty = 1;
                continue
            end
            [T2,p21,p22] = ID(skelmtx2,tol); % skelmtx2(:,p22) = skelmtx1(:,p21) * T2.
            myindex_p21 = myindex_sep2(p21);
            myindex_p22 = myindex_sep2(p22);
            nodekindex_p21 = nodekindex_sep2(p21);
            nodekindex_p22 = nodekindex_sep2(p22);
            obj.nbInfo(k).myindex_p21 = myindex_p21;
            obj.nbInfo(k).myindex_p22 = myindex_p22;
            obj.nbInfo(k).nodekindex_p21 = nodekindex_p21;
            obj.nbInfo(k).nodekindex_p22 = nodekindex_p22;
            nodek.re = [nodek.re,sep2(p22)];
            obj.nbInfo(k).Th2c2 = T2;
            obj.nbre = [obj.nbre,sep2(p22)];
            
            % Step 1
            Ac1h1T1 = obj.ASS(myindex_p11,myindex_p12)'*T1;
            Ah1h1T1 = obj.ASS(myindex_p11,myindex_p11)*T1;
            Ac2h2T2 = nodek.ASS(nodekindex_p21,nodekindex_p22)'*T2;
            Ah2h2T2 = nodek.ASS(nodekindex_p21,nodekindex_p21)*T2;
            % Ac1c1 = Ac1c1 - Ah1c1^{T} * Th1c1 - Th1c1^{T} * Ah1c1 + Th1c1^{T} * Ah1h1 * Th1c1.
            obj.ASS(myindex_p12,myindex_p12) = obj.ASS(myindex_p12,myindex_p12) - Ac1h1T1 - Ac1h1T1' + T1'*Ah1h1T1;
            % Ah1c1 = Ah1c1 - Ah1h1 * Th1c1.
            obj.ASS(myindex_p11,myindex_p12) = obj.ASS(myindex_p11,myindex_p12) - Ah1h1T1;
            obj.ASS(myindex_p12,myindex_p11) = obj.ASS(myindex_p11,myindex_p12)';
            % Ac2c1 = Ac2c1 - Ac2h1 * Th1c1 - Th2c2^{T} * Ah2c1 + Th2c2^{T} * Ah2h1 * Th1c1.
            obj.ANS(myindex_p22,myindex_p12) = obj.ANS(myindex_p22,myindex_p12) - obj.ANS(myindex_p22,myindex_p11)*T1 - ...
                T2'*obj.ANS(myindex_p21,myindex_p12) + T2'*obj.ANS(myindex_p21,myindex_p11)*T1;
            nodek.ANS(nodekindex_p12,nodekindex_p22) = obj.ANS(myindex_p22,myindex_p12)';
            % Ah2c1 = Ah2c1 - Ah2h1 * Th1c1.
            obj.ANS(myindex_p21,myindex_p12) = obj.ANS(myindex_p21,myindex_p12) - obj.ANS(myindex_p21,myindex_p11)*T1;
            nodek.ANS(nodekindex_p12,nodekindex_p21) = obj.ANS(myindex_p21,myindex_p12)';
            % Ac2h1 = Ac2h1 - Th2c2^{T} * Ah2h1.
            obj.ANS(myindex_p22,myindex_p11) = obj.ANS(myindex_p22,myindex_p11) - T2'*obj.ANS(myindex_p21,myindex_p11);
            nodek.ANS(nodekindex_p11,nodekindex_p22) = obj.ANS(myindex_p22,myindex_p11)';
            % Ac2c2 = Ac2c2 - Ah2c2^{T} * Th2c2 - Th2c2^{T} * Ah2c2 + sTh2c2^{T} * Ah2h2 * Th2c2.
            nodek.ASS(nodekindex_p22,nodekindex_p22) = nodek.ASS(nodekindex_p22,nodekindex_p22) - Ac2h2T2 - Ac2h2T2' + T2'*Ah2h2T2;
            % Ah2c2 = Ah2c2 - Ah2h2 * Th2c2.
            nodek.ASS(nodekindex_p21,nodekindex_p22) = nodek.ASS(nodekindex_p21,nodekindex_p22) - Ah2h2T2;
            nodek.ASS(nodekindex_p22,nodekindex_p21) = nodek.ASS(nodekindex_p21,nodekindex_p22)';
            % As11c1 = As12c1 = As21c1 = As22c2 = As12c2 = As21c2 = 0.
            
            % Step 2
            % Ac1c1 = Lc1 * Dc1 * Lc1^{T}.
            [L1,D1] = ldl(obj.ASS(myindex_p12,myindex_p12));
            obj.root.active(sep1(p12)) = 0;
            obj.nbInfo(k).Lc1 = L1;
            obj.nbInfo(k).Dc1 = D1;
            % Ac1c1invAc1h1 = Ac1c1^{-1} * Ah1c1^{T}.
            obj.nbInfo(k).Ac1c1invAc1h1 = L1\(obj.ASS(myindex_p11,myindex_p12)');
            obj.nbInfo(k).Ac1c1invAc1h1 = D1\obj.nbInfo(k).Ac1c1invAc1h1;
            obj.nbInfo(k).Ac1c1invAc1h1 = L1'\obj.nbInfo(k).Ac1c1invAc1h1;
            % Ac1c1invAc1c2 = Ac1c1^{-1} * Ac2c1^{T}.
            obj.nbInfo(k).Ac1c1invAc1c2 = L1\(obj.ANS(myindex_p22,myindex_p12)');
            obj.nbInfo(k).Ac1c1invAc1c2 = D1\obj.nbInfo(k).Ac1c1invAc1c2;
            obj.nbInfo(k).Ac1c1invAc1c2 = L1'\obj.nbInfo(k).Ac1c1invAc1c2;
            % Ac1c1invAc1h2 = Ac1c1^{-1} * Ah2c1^{T}.
            obj.nbInfo(k).Ac1c1invAc1h2 = L1\(obj.ANS(myindex_p21,myindex_p12)');
            obj.nbInfo(k).Ac1c1invAc1h2 = D1\obj.nbInfo(k).Ac1c1invAc1h2;
            obj.nbInfo(k).Ac1c1invAc1h2 = L1'\obj.nbInfo(k).Ac1c1invAc1h2;
            % Ah1h1 = Ah1h1 - Ah1c1 * Ac1c1^{-1} * Ah1c1^{T}.
            obj.ASS(myindex_p11,myindex_p11) = obj.ASS(myindex_p11,myindex_p11) - obj.ASS(myindex_p11,myindex_p12)*obj.nbInfo(k).Ac1c1invAc1h1;
            % Ac2h1 = Ac2h1 - Ac2c1 * Ac1c1^{-1} * Ah1c1^{T}.
            obj.ANS(myindex_p22,myindex_p11) = obj.ANS(myindex_p22,myindex_p11) - obj.ANS(myindex_p22,myindex_p12)*obj.nbInfo(k).Ac1c1invAc1h1;
            nodek.ANS(nodekindex_p11,nodekindex_p22) = obj.ANS(myindex_p22,myindex_p11)';
            % Ah2h1 = Ah2h1 - Ah2c1 * Ac1c1^{-1} * Ah1c1^{T}.
            obj.ANS(myindex_p21,myindex_p11) = obj.ANS(myindex_p21,myindex_p11) - obj.ANS(myindex_p21,myindex_p12)*obj.nbInfo(k).Ac1c1invAc1h1;
            nodek.ANS(nodekindex_p11,nodekindex_p21) = obj.ANS(myindex_p21,myindex_p11)';
            % Ac2c2 = Ac2c2 - Ac2c1 * Ac1c1^{-1} * Ac2c1^{T}.
            nodek.ASS(nodekindex_p22,nodekindex_p22) = nodek.ASS(nodekindex_p22,nodekindex_p22) - obj.ANS(myindex_p22,myindex_p12)*obj.nbInfo(k).Ac1c1invAc1c2;
            % Ah2c2 = Ah2c2 - Ah2c1 * Ac1c1^{-1} * Ac2c1^{T}.
            nodek.ASS(nodekindex_p21,nodekindex_p22) = nodek.ASS(nodekindex_p21,nodekindex_p22) - obj.ANS(myindex_p21,myindex_p12)*obj.nbInfo(k).Ac1c1invAc1c2;
            nodek.ASS(nodekindex_p22,nodekindex_p21) = nodek.ASS(nodekindex_p21,nodekindex_p22)';
            % Ah2h2 = Ah2h2 - Ah2c1 * Ac1c1^{-1} * Ah2c1^{T}.
            nodek.ASS(nodekindex_p21,nodekindex_p21) = nodek.ASS(nodekindex_p21,nodekindex_p21) - obj.ANS(myindex_p21,myindex_p12)*obj.nbInfo(k).Ac1c1invAc1h2;
            % Ah1c1 = Ac2c1 = Ah2c1 = 0.
            
            % Step 3
            % Ac2c2 = Lc2 * Dc2 * Lc2^{T};
            [L2,D2] = ldl(nodek.ASS(nodekindex_p22,nodekindex_p22));
            obj.root.active(sep2(p22)) = 0;
            obj.nbInfo(k).Lc2 = L2;
            obj.nbInfo(k).Dc2 = D2;
            % Ac2c2invAc2h1 = Ac2c2^{-1} * Ac2h1.
            obj.nbInfo(k).Ac2c2invAc2h1 = L2\obj.ANS(myindex_p22,myindex_p11);
            obj.nbInfo(k).Ac2c2invAc2h1 = D2\obj.nbInfo(k).Ac2c2invAc2h1;
            obj.nbInfo(k).Ac2c2invAc2h1 = L2'\obj.nbInfo(k).Ac2c2invAc2h1;
            % Ac2c2invAc2h2 = Ac2c2^{-1} * Ah2c2^{T}.
            obj.nbInfo(k).Ac2c2invAc2h2 = L2\(nodek.ASS(nodekindex_p21,nodekindex_p22)');
            obj.nbInfo(k).Ac2c2invAc2h2 = D2\obj.nbInfo(k).Ac2c2invAc2h2;
            obj.nbInfo(k).Ac2c2invAc2h2 = L2'\obj.nbInfo(k).Ac2c2invAc2h2;
            % Ah1h1 = Ah1h1 - Ac2h1^{T} * Ac2c2^{-1} * Ac2h1.
            obj.ASS(myindex_p11,myindex_p11) = obj.ASS(myindex_p11,myindex_p11) - obj.ANS(myindex_p22,myindex_p11)'*obj.nbInfo(k).Ac2c2invAc2h1;
            % Ah2h1 = Ah2h1 - Ah2c2 * Ac2c2^{-1} * Ac2h1.
            obj.ANS(myindex_p21,myindex_p11) = obj.ANS(myindex_p21,myindex_p11) - nodek.ASS(nodekindex_p21,nodekindex_p22)*obj.nbInfo(k).Ac2c2invAc2h1;
            nodek.ANS(nodekindex_p11,nodekindex_p21) = obj.ANS(myindex_p21,myindex_p11)';
            % Ah2h2 = Ah2h2 - Ah2c2 * Ac2c2^{-1} * Ah2c2^{T}.
            nodek.ASS(nodekindex_p21,nodekindex_p21) = nodek.ASS(nodekindex_p21,nodekindex_p21)- nodek.ASS(nodekindex_p21,nodekindex_p22)*obj.nbInfo(k).Ac2c2invAc2h2;
            % Ah2c2 = Ac2h1 = 0.
            
            obj.nbInfo(k).empty = 0;
        end
        
        obj.re = sort(obj.re);
        obj.sk = setdiff(obj.sep,obj.re,'sorted');
        obj.nbre = sort(obj.nbre);
        
        end
        
        function obj = NoSkel(obj)
        % NoSkel No skeletonization.
        
        obj.sk = obj.sep;
        
        end
        
        function obj = RecursiveMerge(obj,whatlevel)
        % RecursiveMerge Recusively send information from children to parent.
        
        if obj.level == whatlevel
            obj = Merge(obj);
        else
            for iter = [1,2]
                if ~isempty(obj.children{iter})
                    obj.children{iter} = RecursiveMerge(obj.children{iter},whatlevel);
                end
            end
        end
        
        end
        
        function obj = Merge(obj)
        % Merge Send matrices' information from children to parent.
        
        % We stand on the parent level.
        
        if obj.endFlag == 1
            return;
        end
        
        % First we tell the parent what its int, sep, nb is after we
        % eliminate the children's vtx.
        % int: children's sk - sep.
        % sep: sep \cup children's sk.
        % nb: nb - children's nbre.
        
        for iter = [1,2]
            obj.int = [obj.int,obj.children{iter}.sk];
            obj.re = [obj.re,obj.children{iter}.re];
            obj.nbre = [obj.nbre,obj.children{iter}.nbre];
        end
        
        obj.sep = intersect(obj.sep,obj.int,'sorted');
        obj.nbre = setdiff(obj.nbre,obj.vtx);
        obj.int = setdiff(obj.int,obj.sep,'sorted');
        obj.nb = setdiff(obj.nb,obj.nbre,'sorted');
        
        
        % Next we assign the corresponding matrices.
        % From child to parent.
        % int: children's sep and children's nb
        % sep: children's sep and children's nb
        % nb: children's nb
        
        % We assign values blockly. The details can be found in the for-for
        % loop later.
        
        obj.AII = zeros(length(obj.int));
        [int1,myindex_int1,~] = intersect(obj.int,obj.children{1}.vtx);
        [~,cindex_int1] = ismember(int1,obj.children{1}.sep);
        [int2,myindex_int2,~] = intersect(obj.int,obj.children{2}.vtx);
        [~,cindex_int2] = ismember(int2,obj.children{2}.sep);
        [~,myindex_int21,cindex_int21] = intersect(obj.int,obj.children{1}.nb);
        obj.indexInfo(1).myindex_int = myindex_int1;
        obj.indexInfo(1).cindex_int = cindex_int1;
        obj.indexInfo(2).myindex_int = myindex_int2;
        obj.indexInfo(2).cindex_int = cindex_int2;
        obj.AII(myindex_int1,myindex_int1) = obj.children{1}.ASS(cindex_int1,cindex_int1);
        obj.AII(myindex_int2,myindex_int2) = obj.children{2}.ASS(cindex_int2,cindex_int2);
        obj.AII(myindex_int21,myindex_int1) = obj.children{1}.ANS(cindex_int21,cindex_int1);
        obj.AII(myindex_int1,myindex_int2) = obj.AII(myindex_int2,myindex_int1)';
                
        obj.ASI = zeros(length(obj.sep),length(obj.int));
        [~,myindex_sep1x,cindex_sep1x] = intersect(obj.sep,obj.children{1}.sep);
        [~,myindex_sep1y,cindex_sep1y] = intersect(obj.sep,obj.children{1}.nb);
        [~,myindex_sep2x,cindex_sep2x] = intersect(obj.sep,obj.children{2}.sep);
        [~,myindex_sep2y,cindex_sep2y] = intersect(obj.sep,obj.children{2}.nb);
        obj.ASI(myindex_sep1x,myindex_int1) = obj.children{1}.ASS(cindex_sep1x,cindex_int1);
        obj.ASI(myindex_sep2x,myindex_int2) = obj.children{2}.ASS(cindex_sep2x,cindex_int2);
        obj.ASI(myindex_sep1y,myindex_int1) = obj.children{1}.ANS(cindex_sep1y,cindex_int1);
        obj.ASI(myindex_sep2y,myindex_int2) = obj.children{2}.ANS(cindex_sep2y,cindex_int2);
        
        obj.ASS = zeros(length(obj.sep));
        [sep1,myindex_sep1,~] = intersect(obj.sep,obj.children{1}.vtx);
        [~,cindex_sep1] = ismember(sep1,obj.children{1}.sep);
        [sep2,myindex_sep2,~] = intersect(obj.sep,obj.children{2}.vtx);
        [~,cindex_sep2] = ismember(sep2,obj.children{2}.sep);
        [~,myindex_sep21,cindex_sep21] = intersect(obj.sep,obj.children{1}.nb);
        obj.indexInfo(1).myindex_sep = myindex_sep1;
        obj.indexInfo(1).cindex_sep = cindex_sep1;
        obj.indexInfo(2).myindex_sep = myindex_sep2;
        obj.indexInfo(2).cindex_sep = cindex_sep2;
        obj.ASS(myindex_sep1,myindex_sep1) = obj.children{1}.ASS(cindex_sep1,cindex_sep1);
        obj.ASS(myindex_sep2,myindex_sep2) = obj.children{2}.ASS(cindex_sep2,cindex_sep2);
        obj.ASS(myindex_sep21,myindex_sep1) = obj.children{1}.ANS(cindex_sep21,cindex_sep1);
        obj.ASS(myindex_sep1,myindex_sep2) = obj.ASS(myindex_sep2,myindex_sep1)';
               
        obj.ANS = zeros(length(obj.nb),length(obj.sep));
        [~,myindex_nb1x,cindex_nb1x] = intersect(obj.nb,obj.children{1}.nb);
        [~,myindex_nb2x,cindex_nb2x] = intersect(obj.nb,obj.children{2}.nb);
        obj.ANS(myindex_nb1x,myindex_sep1) = obj.children{1}.ANS(cindex_nb1x,cindex_sep1);
        obj.ANS(myindex_nb2x,myindex_sep2) = obj.children{2}.ANS(cindex_nb2x,cindex_sep2);
        
        %%% ------------------------Details---------------------------- %%%
%         obj.AII = zeros(length(obj.int));
%         % An int of the parent only belongs to the sep of one of its
%         % children. If two ints belong to the same child, we assign AII
%         % from the child's ASS. Otherwise, we assign AII from one child's
%         % ANS or 0.
%         for j = 1:length(obj.int)
%             intj = obj.int(j);
%             if find(obj.children{1}.vtx == intj)
%                 where_intj = 1;
%             else
%                 where_intj = 2;
%             end
%             index_intj = find(obj.children{where_intj}.sep == intj);
%             for i = 1:length(obj.int)
%                 inti = obj.int(i);
%                 index_inti = find(obj.children{where_intj}.sep == inti,1);
%                 if isempty(index_inti)
%                     index_inti = find(obj.children{where_intj}.nb == inti,1);
%                     if isempty(index_inti)
%                         obj.AII(i,j) = 0;
%                     else
%                         obj.AII(i,j) = obj.children{where_intj}.ANS(index_inti,index_intj);
%                     end
%                 else
%                     obj.AII(i,j) = obj.children{where_intj}.ASS(index_inti,index_intj);
%                 end
%             end
%         end
%         
%         obj.ASI = zeros(length(obj.sep),length(obj.int));
%         % A sep of the parent only belongs to the sep of one of its
%         % children. If an int and a sep belongs to the same child, we
%         % assign ASI from the child's ASS. Otherwise, we assign ASI from
%         % one child's ANS or 0.
%         for j = 1:length(obj.int)
%             intj = obj.int(j);
%             if find(obj.children{1}.vtx == intj)
%                 where_intj = 1;
%             else
%                 where_intj = 2;
%             end
%             index_intj = find(obj.children{where_intj}.sep == intj);
%             for i = 1:length(obj.sep)
%                 sepi = obj.sep(i);
%                 index_sepi = find(obj.children{where_intj}.sep == sepi,1);
%                 if isempty(index_sepi)
%                     index_sepi = find(obj.children{where_intj}.nb == sepi,1);
%                     if isempty(index_sepi)
%                         obj.ASI(i,j) = 0;
%                     else
%                         obj.ASI(i,j) = obj.children{where_intj}.ANS(index_sepi,index_intj);
%                     end
%                 else
%                     obj.ASI(i,j) = obj.children{where_intj}.ASS(index_sepi,index_intj);
%                 end
%             end
%         end
%         
%         obj.ASS = zeros(length(obj.sep));
%         % If two seps belongs to the same child, we assign ASS from the
%         % child's ASS. Otherwise, we assign ASS from the child's ANS or 0
%         % whose seqNum is smaller.
%         for j = 1:length(obj.sep)
%             sepj = obj.sep(j);
%             if find(obj.children{1}.vtx == sepj)
%                 where_sepj = 1;
%             else
%                 where_sepj = 2;
%             end
%             index_sepj = find(obj.children{where_sepj}.sep == sepj);
%             for i = 1:length(obj.sep)
%                 sepi = obj.sep(i);
%                 index_sepi = find(obj.children{where_sepj}.sep == sepi,1);
%                 if isempty(index_sepi)
%                     index_sepi = find(obj.children{where_sepj}.nb == sepi,1);
%                     if isempty(index_sepi)
%                         obj.ASS(i,j) = 0;
%                     else
%                         obj.ASS(i,j) = obj.children{where_sepj}.ANS(index_sepi,index_sepj);
%                     end
%                 else
%                     obj.ASS(i,j) = obj.children{where_sepj}.ASS(index_sepi,index_sepj);
%                 end
%             end
%         end
%         
%         obj.ANS = zeros(length(obj.nb),length(obj.sep));
%         % If a nb and a sep in the same child, we assign ANS from the
%         % child's ANS. Otherwise, ANS= 0.
%         for j = 1:length(obj.sep)
%             sepj = obj.sep(j);
%             if find(obj.children{1}.vtx == sepj)
%                 where_sepj = 1;
%             else
%                 where_sepj = 2;
%             end
%             index_sepj = find(obj.children{where_sepj}.sep == sepj);
%             for i = 1:length(obj.nb)
%                 nbi = obj.nb(i);
%                 index_nbi = find(obj.children{where_sepj}.nb == nbi,1);
%                 if isempty(index_nbi)
%                     obj.ANS(i,j) = 0;
%                 else
%                     obj.ANS(i,j) = obj.children{where_sepj}.ANS(index_nbi,index_sepj);
%                 end
%             end
%         end
        %%% ----------------------------------------------------------- %%%
        
        % Clear children information.
        for iter = [1,2]
            obj.children{iter} = HIFClear(obj.children{iter});
        end
        
        % Set sep type.
        obj = SetSepType(obj);
        
        end
        
        function obj = HIFClear(obj)
        % HIFClear Clear unnecessary information.
        
        % TODO: HIFClear.
        
        end
        
        function obj = RootFactorization(obj)
        % RootFactorization Factorization on the root.
        
        obj.root.active(obj.int) = 0;
        % AII = LI * DI * LI^{T}.
        [obj.LI,obj.DI] = ldl(obj.AII);
        
        end
        
        function obj = HIFSolve(obj,b)
        % HIFSolve Solve Ax = b through HIF.
        
        disp("  ");
        disp(" Start solve ");
        disp("  ");
        
        obj.inputVec = b;
        
        obj = BuildVecTree(obj,b);
        
        for tmplevel = obj.numLevels:-1:1
            obj = RecursiveApplySparseElimUp(obj,tmplevel);
            obj = RecursiveApplySkelUp(obj,tmplevel);
            obj = RecursiveApplyMerge(obj,tmplevel-1);
        end
        
        obj = RootApply(obj);
        
        for tmplevel = 1:1:obj.numLevels
            obj = RecursiveApplySplit(obj,tmplevel-1);
            obj = RecursiveApplySkelDown(obj,tmplevel);
            obj = RecursiveApplySparseElimDown(obj,tmplevel);
        end
        
        obj = GetSolution(obj);
        
        disp("  ");
        disp(" End solve ");
        disp(" The solution is stored in obj.solution.");
        disp("  ");
        
        end
        
        function obj = BuildVecTree(obj,b)
        % BuildVecTree Fill b in our HIF tree.
        
        if obj.level == 0
            disp("  ");
            disp(" Start build vector tree ");
            disp("  ");
        end
        
        if obj.endFlag == 0
            for iter = [1,2]
                obj.children{iter} = BuildVecTree(obj.children{iter},b);
            end
        else
            obj.xI = b(obj.int);
            obj.xS = b(obj.sep);
        end
        
        if obj.level == 0
            disp("  ");
            disp(" End build vector tree ");
            disp("  ");
        end
        
        end
        
        function obj = RecursiveApplySparseElimUp(obj,whatlevel)
        % RecursiveApplySparseElimUp Phase 1 for applying sparse elimination recusively.
        
        if obj.level == whatlevel
            obj = ApplySparseElimUp(obj);
        else
            for iter = [1,2]
                if ~isempty(obj.children{iter})
                    obj.children{iter} = RecursiveApplySparseElimUp(obj.children{iter},whatlevel);
                end
            end
        end
        
        end
        
        function obj = ApplySparseElimUp(obj)
        %ApplySparseElimUp Phase 1 for applying sparse elimination.
        
        % xS = xS - (AII^{-1} * ASI^{T})^{T} * xI.
        obj.xS = obj.xS - obj.AIIinvAIS'*obj.xI;
        % xI = LI^{-1} * xI.
        obj.xI = obj.LI\obj.xI;
        
        % xI = DI^{-1} * xI. We only apply D once.
        obj.xI = obj.DI\obj.xI;
        
        end
        
        function obj = RecursiveApplySkelUp(obj,whatlevel)
        % RecursiveApplySkelUp Phase 1 for applying skeletonization recusively.
        
        if obj.level == whatlevel
            obj = ApplySkelUp(obj);
        else
            for iter = [1,2]
                if ~isempty(obj.children{iter})
                    obj.children{iter} = RecursiveApplySkelUp(obj.children{iter},whatlevel);
                end
            end
        end
        
        end
        
        function obj = ApplySkelUp(obj)
        % ApplySkelUp Phase 1 for applying skeletonization.
        
        for k = 1:length(obj.nbNode)
            nbnodek = obj.nbNode{k};
            if isempty(obj.nbInfo)
                continue;
            end
            nbinfok = obj.nbInfo(k);
            if nbinfok.empty
                continue;
            end
            
            % Step 1
            % xc1 = xc1 - Th1c1^{T} * xh1.
            obj.xS(nbinfok.myindex_p12) = obj.xS(nbinfok.myindex_p12) - nbinfok.Th1c1'*obj.xS(nbinfok.myindex_p11);
            % xc2 = xc2 - Th2c2^{T} * xh2.
            nbnodek.xS(nbinfok.nodekindex_p22) = nbnodek.xS(nbinfok.nodekindex_p22) - nbinfok.Th2c2'*nbnodek.xS(nbinfok.nodekindex_p21);
            
            % Step 2
            % xh1 = xh1 - (Ac1c1^{-1} * Ah1c1^{T})^{T} * xc1.
            obj.xS(nbinfok.myindex_p11) = obj.xS(nbinfok.myindex_p11) - nbinfok.Ac1c1invAc1h1'*obj.xS(nbinfok.myindex_p12);
            % xc2 = xc2 - (Ac1c1^{-1} * Ac2c1^{T})^{T} * xc1.
            nbnodek.xS(nbinfok.nodekindex_p22) = nbnodek.xS(nbinfok.nodekindex_p22) - nbinfok.Ac1c1invAc1c2'*obj.xS(nbinfok.myindex_p12);
            % xh2 = xh2 - (Ac1c1^{-1} * Ah2c1^{T})^{T} * xc1.
            nbnodek.xS(nbinfok.nodekindex_p21) = nbnodek.xS(nbinfok.nodekindex_p21) - nbinfok.Ac1c1invAc1h2'*obj.xS(nbinfok.myindex_p12);
            % xc1 = Lc1^{-1} * xc1.
            obj.xS(nbinfok.myindex_p12) = nbinfok.Lc1\obj.xS(nbinfok.myindex_p12);
            
            % Step 3
            % xh1 = xh1 - (Ac2c2^{-1} * Ac2h1)^{T} * xc2.
            obj.xS(nbinfok.myindex_p11) = obj.xS(nbinfok.myindex_p11) - nbinfok.Ac2c2invAc2h1'*nbnodek.xS(nbinfok.nodekindex_p22);
            % xh2 = xh2 - (Ac2c2^{-1} * Ah2c2^{T})^{T} * xc2.
            nbnodek.xS(nbinfok.nodekindex_p21) = nbnodek.xS(nbinfok.nodekindex_p21) - nbinfok.Ac2c2invAc2h2'*nbnodek.xS(nbinfok.nodekindex_p22);
            % xc2 = Lc2^{-1} * xc2.
            nbnodek.xS(nbinfok.nodekindex_p22) = nbinfok.Lc2\nbnodek.xS(nbinfok.nodekindex_p22);
            
            % xc1 = Dc1^{-1} * xc1. xc2 = Dc2^{-1} * xc2. We only apply D once.
            obj.xS(nbinfok.myindex_p12) = nbinfok.Dc1\obj.xS(nbinfok.myindex_p12);
            nbnodek.xS(nbinfok.nodekindex_p22) = nbinfok.Dc2\nbnodek.xS(nbinfok.nodekindex_p22);
        end
        
        end
        
        function obj = RecursiveApplyMerge(obj,whatlevel)
        % RecursiveApplyMerge Recusively send vectors' information from children to parent.
        
        if obj.level == whatlevel
            obj = ApplyMerge(obj);
        else
            for iter = [1,2]
                if ~isempty(obj.children{iter})
                    obj.children{iter} = RecursiveApplyMerge(obj.children{iter},whatlevel);
                end
            end
        end
        
        end
        
        function obj = ApplyMerge(obj)
        % ApplyMerge Send vectors' information from children to parent.
        
        % We stand on the parent level.
        
        if obj.endFlag == 1
            return;
        end
        
        % We only need to assign the corresponding vectors.
        
        obj.xI = zeros(length(obj.int),1);
        for iter = [1,2]
            obj.xI(obj.indexInfo(iter).myindex_int) = obj.children{iter}.xS(obj.indexInfo(iter).cindex_int);
        end
        
        obj.xS = zeros(length(obj.sep),1);
        for iter = [1,2]
            obj.xS(obj.indexInfo(iter).myindex_sep) = obj.children{iter}.xS(obj.indexInfo(iter).cindex_sep);
        end
        
        %%% ------------------------Details---------------------------- %%%
%         obj.xI = zeros(length(obj.int),1);
%         % An int of the parent only belongs to the sep of one of its
%         % children. We get xI from the children's xS.
%         for j =1:length(obj.int)
%             intj = obj.int(j);
%             if find(obj.children{1}.vtx == intj)
%                 where_intj = 1;
%             else
%                 where_intj = 2;
%             end
%             index_intj = find(obj.children{where_intj}.sep == intj);
%             obj.xI(j) = obj.children{where_intj}.xS(index_intj);
%         end
%         
%         obj.xS = zeros(length(obj.sep),1);
%         % A sep of the parent only belongs to the sep of one of its
%         % children. We get xS from the children's xS.
%         for j = 1:length(obj.sep)
%             sepj = obj.sep(j);
%             if find(obj.children{1}.vtx == sepj)
%                 where_sepj = 1;
%             else
%                 where_sepj = 2;
%             end
%             index_sepj = find(obj.children{where_sepj}.sep == sepj);
%             obj.xS(j) = obj.children{where_sepj}.xS(index_sepj);
%         end
        %%% ----------------------------------------------------------- %%%
        
        end
        
        function obj = RootApply(obj)
        % RootApply Apply on the root.
        
        % xI = LI^{-1} * xI.
        obj.xI = obj.LI\obj.xI;
        
        % xI = DI^{-1} * xI.
        obj.xI = obj.DI\obj.xI;
        
        % xI = LI^{-T} * xI.
        obj.xI = obj.LI'\obj.xI;
        
        end
        
        function obj = RecursiveApplySplit(obj,whatlevel)
        % RecursiveApplySplit Recusively send vectors' information from parent to children.
        
        if obj.level == whatlevel
            obj = ApplySplit(obj);
        else
            for iter = [1,2]
                if ~isempty(obj.children{iter})
                    obj.children{iter} = RecursiveApplySplit(obj.children{iter},whatlevel);
                end
            end
        end
        
        end
        
        function obj = ApplySplit(obj)
        % ApplySplit Send vectors' information from parent to children.
        
        % We stand on the parent level.
        
        if obj.endFlag == 1
            return;
        end
        
        % We only need to assign the corresponding vectors of the children.
        
        for iter = [1,2]
            obj.children{iter}.xS(obj.indexInfo(iter).cindex_int) = obj.xI(obj.indexInfo(iter).myindex_int);
            obj.children{iter}.xS(obj.indexInfo(iter).cindex_sep) = obj.xS(obj.indexInfo(iter).myindex_sep);
        end
        
        %%% ------------------------Details---------------------------- %%%
%         % xI
%         for j = 1:length(obj.int)
%             intj = obj.int(j);
%             if find(obj.children{1}.vtx == intj)
%                 where_intj = 1;
%             else
%                 where_intj = 2;
%             end
%             index_intj = find(obj.children{where_intj}.sep == intj);
%             obj.children{where_intj}.xS(index_intj) = obj.xI(j);
%         end
%        
%         % xS
%         for j = 1:length(obj.sep)
%             sepj = obj.sep(j);
%             if find(obj.children{1}.vtx == sepj)
%                 where_sepj = 1;
%             else
%                 where_sepj = 2;
%             end
%             index_sepj = find(obj.children{where_sepj}.sep == sepj);
%             obj.children{where_sepj}.xS(index_sepj) = obj.xS(j);
%         end
        %%% ----------------------------------------------------------- %%%
        
        end
        
        function obj = RecursiveApplySkelDown(obj,whatlevel)
        % RecursiveApplySkelDown Phase 2 for applying skeletonization recursively.
        
        if obj.level == whatlevel
            obj = ApplySkelDown(obj);
        else
            for iter = [1,2]
                if ~isempty(obj.children{iter})
                    obj.children{iter} = RecursiveApplySkelDown(obj.children{iter},whatlevel);
                end
            end
        end
        
        end
        
        function obj = ApplySkelDown(obj)
        % ApplySkelDown Phase 2 for applying skeletonization.
        
        for k = 1:length(obj.nbNode)
            nbnodek = obj.nbNode{k};
            if isempty(obj.nbInfo)
                continue;
            end
            nbinfok = obj.nbInfo(k);
            if nbinfok.empty
                continue;
            end
            
            % Step 3
            % xc2 = Lc2^{-T} * xc2 - (Ac2c2^{-1} * Ac2h1) * xh1 - (Ac2c2^{-1} * Ah2c2^{T}) * xh2.
            nbnodek.xS(nbinfok.nodekindex_p22) = nbinfok.Lc2'\nbnodek.xS(nbinfok.nodekindex_p22) - ...
                nbinfok.Ac2c2invAc2h1*obj.xS(nbinfok.myindex_p11) - ...
                nbinfok.Ac2c2invAc2h2*nbnodek.xS(nbinfok.nodekindex_p21);
            
            % Step 2
            % xc1 = Lc1^{-T} * xc1 - (Ac1c1^{-1} * Ah1c1^{T}) * xh1 - (Ac1c1^{-1} * Ac2c1^{T}) * xc2 - (Ac1c1^{-1} * Ah2c1^{T}) * xh2.
            obj.xS(nbinfok.myindex_p12) = nbinfok.Lc1'\obj.xS(nbinfok.myindex_p12) - ...
                nbinfok.Ac1c1invAc1h1*obj.xS(nbinfok.myindex_p11) - ...
                nbinfok.Ac1c1invAc1c2*nbnodek.xS(nbinfok.nodekindex_p22) - ...
                nbinfok.Ac1c1invAc1h2*nbnodek.xS(nbinfok.nodekindex_p21);
            
            % Step 1
            % xh1 = xh1 - Th1c1 * xc1.
            obj.xS(nbinfok.myindex_p11) = obj.xS(nbinfok.myindex_p11) - nbinfok.Th1c1*obj.xS(nbinfok.myindex_p12);
            % xh2 = xh2 - Th2c2 * xc2.
            nbnodek.xS(nbinfok.nodekindex_p21) = nbnodek.xS(nbinfok.nodekindex_p21) - nbinfok.Th2c2*nbnodek.xS(nbinfok.nodekindex_p22);
        end
        
        
        end
        
        function obj = RecursiveApplySparseElimDown(obj,whatlevel)
        % RecursiveApplySparseElimDown Phase 2 for applying sparse elimination recursively.
        
        if obj.level == whatlevel
            obj = ApplySparseElimDown(obj);
        else
            for iter = [1,2]
                if ~isempty(obj.children{iter})
                    obj.children{iter} = RecursiveApplySparseElimDown(obj.children{iter},whatlevel);
                end
            end
        end
        
        end
        
        function obj = ApplySparseElimDown(obj)
        % ApplySparseElimDown Phase 2 for applying sparse elimination.
        
        % xI = LI^{-T} * xI - (AII^{-1} * ASI^{T})^{T} * xS.
        obj.xI = obj.LI'\obj.xI - obj.AIIinvAIS*obj.xS;
        
        end
        
        function obj = GetSolution(obj)
        % GetSolution Get the final soluion through the tree structure.
        
        if obj.level == 0
            obj.solution = zeros(length(obj.int),1);
        end
        
        if obj.endFlag == 0
            for iter = [1,2]
                obj.children{iter} = GetSolution(obj.children{iter});
            end
        else
            obj.root.solution(obj.int) = obj.xI;
            obj.root.solution(obj.sep) = obj.xS;
        end
        
        end
        
        function DemoPart(obj)
        % DemoPart Demo the partition process.
        
        assert(~isempty(obj.inputAxy.xy),"Need coordinates!")
        
        disp("  ");
        disp(" Start demo the partition ");
        disp("  ");
        
        figure(1);
        clf reset;
        colordef(1,'black');
        gplotg(obj.inputAxy.A,obj.inputAxy.xy);
        
        disp(" Hit space to continue ... ");
        disp("  ");
        
        for tmplevel = 0:obj.numLevels
            disp(" Current level: " + tmplevel);
            disp("  ");
            
            map = GetPartMap(obj,tmplevel);
            gplotmap(obj.inputAxy.A,obj.inputAxy.xy,map);
            if tmplevel ~= obj.numLevels
                disp(" Hit space to continue ... ");
            else
                disp(" Hit space to end ... ");
            end
            disp("  ");
            pause;
        end
        
        end
        
        function map = GetPartMap(obj,whatlevel,map)
        % GetPartMap Get the map of the partition.
        
        if nargin == 2
            map = [];
        end
        
        if obj.level == 0
            n = size(obj.inputAxy.xy,1);
            map = ones(1,n);
        end
        
        if obj.level == whatlevel
            map(obj.vtx) = obj.seqNum;
            return;
        else
            for iter = [1,2]
                if ~isempty(obj.children{iter})
                    map = GetPartMap(obj.children{iter},whatlevel,map);
                end
            end
        end
        
        end
        
        function DemoFinalPart(obj)
        % DemoFinalPart Demo of the final partition.
        
        assert(~isempty(obj.inputAxy.xy),"Need coordinates!")
        
        disp("  ");
        disp(" Start demo the final partition ");
        disp("  ");
        
        figure(1);
        clf reset;
        colordef(1,'black');
        gplotg(obj.inputAxy.A,obj.inputAxy.xy);
        
        disp(" Hit space to continue ... ");
        disp("  ");
        pause;
        map = GetPartMap(obj,obj.numLevels);
        gplotmap(obj.inputAxy.A,obj.inputAxy.xy,map);
        disp(" Hit space to end ... ");
        disp("  ");
        pause;
        
        end
        
        function DemoHIF(obj,whatlevel)
        % DemoHIF Demo the HIF process.
        
        assert(~isempty(obj.inputAxy.xy), "Need coordinates!")
        
        if whatlevel == obj.numLevels
            DemoFinalPart(obj);
        end
        
        disp("  ");
        disp(" Current level: " + whatlevel);
        disp("  ");
                
        map = GetHIFMap(obj,whatlevel);
        inactive = find(obj.active == 0);
        start = max(map) + 1;
        n = length(inactive);
        map(inactive) = start:start+n-1;
        gplotmap(obj.inputAxy.A,obj.inputAxy.xy,map);
        if whatlevel ~= 0
            disp(" Hit space to continue ...");
        else
            disp(" Hit space to end ...");
        end
        disp("  ");
        pause;
        
        end
        
        function map = GetHIFMap(obj,whatlevel,map)
        % GetHIFMap Recursively get HIF map.
        
        if nargin == 2
            map = [];
        end
        
        if obj.level == 0
            n = size(obj.inputAxy.xy,1);
            map = ones(1,n);
        end
        
        if obj.level == whatlevel
            map(obj.vtx) = obj.seqNum;
            return;
        else
            for iter = [1,2]
                map = GetHIFMap(obj.children{iter},whatlevel,map);
            end
        end
        
        end
        
    end
    
end
