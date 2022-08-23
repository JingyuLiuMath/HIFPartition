classdef HIFGraph < handle
    % HIFGraph HIF algorithm based on graph partition.
    
    properties
        
        % Root properties.
        % The following information will be stored only in the root node.
        inputAxy; % Input Axy.
        active; % Whether a vertex is eliminated.
        demoHIF = 0; % Whether to demo the HIF process.
        
        % Graph properties.
        vtx; % Vertices.
        sep; % Separator vertices.
        nb; % Neighbor vertices.
        int; % Interior vertices.
        sk; % Skeleton sep. We also use hat (h) to represent it.
        re; % Redundant sep. We also use check (c) to reprsent it.
        nbsk; % Skeleton nb.
        nbre; % Redundant nb.
        singlesep = {}; % Sep which only interact with one node.
        complexsep; % Sep which interact with more than one node.
        
        % Tree properties.
        numlevels; % The highest level in current subtree.
        level; % Current level, start from 0.
        seqnum; % A node's order in its level.
        endflag = 0; % Whether the partition ends.
        children = cell(1,2); % Children nodes.
        nbnode = {}; % Neighbor nodes.
        nbnodeseqnum = []; % Neighbor nodes' seqnum.
        nbnodelevel = []; % Neighbor nodes' level.
        root; % Root node.
        nbinfo = struct([]); % Information between a node and its nbnode when skeletonization.
        copy = 0; % Whether it is a copy HIFGraph of its parerent.
        havecopy = 0; % Whether it has a copy.
        indexinfo = struct([]); % Index information of a node and its children.
        
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
        
        function obj = HIFGraph(Axy,minvtx,method,level,seqnum,vtx,sep,nb)
        % HIFGraph Create a HIF class.
        
        if nargin == 1
            minvtx = 64;
            method = "metis";
        end
        
        if nargin == 2
            method = "metis";
        end
        
        if nargin <= 3
            level = 0;
            seqnum = 0;
            n = size(Axy.A,1);
            vtx = 1:1:n;
            sep = [];
            nb = [];
            
            obj.root = obj;
            obj.inputAxy = Axy;
            obj.active = ones(1,n);
        end
        
        obj.level = level;
        obj.seqnum = seqnum;
        obj.vtx = vtx;
        obj.sep = sep;
        obj.nb = nb;
        
        if obj.level == 0
            disp("  ");
            disp(" Start initialization ");
            disp("  ");
            obj = BuildTree(obj,Axy,minvtx,method);
            obj = SetNeighborNode(obj);
            obj = FillTree(obj,Axy.A);
            disp("  ");
            disp(" End initialization ");
            disp("  ");
        end
        
        end
        
        function obj = BuildTree(obj,Axy,minvtx,method)
        % BuildTree Build tree structure according to a graph partition algorithm.
        
        n = length(obj.vtx);
        
        if n <= minvtx
            obj.numlevels = obj.level;
            obj.endflag = 1;
            return
        end
        
        % Partition
        tmpAxy.A = Axy.A(obj.vtx,obj.vtx);
        if ~isempty(Axy.xy)
            tmpAxy.xy = Axy.xy(obj.vtx,:);
        end
        [p1,p2,sep1,sep2] = GraphPart(tmpAxy,method);
        % sep1 = unique(sep1);
        % sep2 = unique(sep2);
        p = {p1,p2};
        partsep = {sep1,sep2};
        
        % Create children HIF.
        for iter = [1,2]
            obj.children{iter} = HIFGraph(Axy,minvtx,method,obj.level+1,obj.seqnum*2+(iter-1),...
                obj.vtx(p{iter}),obj.vtx(partsep{iter}),obj.vtx(partsep{3-iter}));
            obj.children{iter}.root = obj.root;
        end
        
        % Pass information to its children.
        obj = PassSeparatorNeighbor(obj,Axy.A);
        
        % Recursively buildtree.
        for iter = [1,2]
            obj.children{iter} = BuildTree(obj.children{iter},Axy,minvtx,method);
        end
        
        % Get numlevels from children when partition ends.
        obj.numlevels = max(obj.children{1}.numlevels,obj.children{2}.numlevels);
        
        end
        
        function obj = PassSeparatorNeighbor(obj,A)
        % PassSeparatorNeighbor Send parent's sep, nb to children.
        
        nbA = A(obj.sep,obj.nb);
        addsep1 = [];
        addnb1 = [];
        addsep2 = [];
        addnb2 = [];
        for i = 1:length(obj.sep)
            sepi = obj.sep(i);
            for iter = [1, 2]
                obj_child = obj.children{iter};
                index_sepi_childnode = find(obj_child.vtx == sepi,1);
                if isempty(index_sepi_childnode)
                    continue;
                end
                index_addnb = find(nbA(i,:)~=0);
                addnb = obj.nb(index_addnb);
                if iter == 1
                    addsep1 = [addsep1,sepi];
                    addnb1 = [addnb1, addnb];
                else
                    addsep2 = [addsep2,sepi];
                    addnb2 = [addnb2, addnb];
                end
            end
        end
        obj.children{1}.sep = sort(unique([obj.children{1}.sep,addsep1]));
        obj.children{2}.sep = sort(unique([obj.children{2}.sep,addsep2]));
        obj.children{1}.nb = sort(unique([obj.children{1}.nb,addnb1]));
        obj.children{2}.nb = sort(unique([obj.children{2}.nb,addnb2]));
        
        end
        
        function obj = SetNeighborNode(obj)
        % SetNeighborNode Set nbnode.
        
        % We stand on the parent level to assign its children's nbnode.
        if obj.endflag == 1
            return;
        end
        
        for iter = [1,2]
            obj_child = obj.children{iter};
            obj_child.nbnode{end+1} = obj.children{3-iter};
            obj_child.nbnodeseqnum(end+1) = obj.children{3-iter}.seqnum;
            obj_child.nbnodelevel(end+1) = obj.children{3-iter}.level;
            % We only need to find nbnode from the children node of
            % parent's nbnode or parent's nbnode if it doesn't have a
            % child.
            for i = 1:length(obj.nbnode)
                nbnodei = obj.nbnode{i};
                % What we need is to check whether the vtx of nbnodei's
                % chilren is in the nb of obj_child.
                if nbnodei.endflag == 1
                    if nbnodei.havecopy == 0
                        nbnodei.children{1} = HIFGraph([],[],[],nbnodei.level+1,nbnodei.seqnum*2,...
                            nbnodei.vtx,nbnodei.sep,nbnodei.nb);
                        nbnodei.children{1}.endflag = 1;
                        nbnodei.children{1}.copy = 1;  
                        nbnodei.havecopy = 1;
                    end
                    nbnodei_child = nbnodei.children{1};
                    if ~isempty(intersect(obj_child.nb,nbnodei_child.sep))
                        obj_child.nbnode{end+1} = nbnodei_child;
                        obj_child.nbnodeseqnum(end+1) = nbnodei_child.seqnum;
                        obj_child.nbnodelevel(end+1) = nbnodei_child.level;
                    end
                else
                    for it = [1,2]
                        nbnodei_child = nbnodei.children{it};
                        if ~isempty(intersect(obj_child.nb,nbnodei_child.sep))
                            obj_child.nbnode{end+1} = nbnodei_child;
                            obj_child.nbnodeseqnum(end+1) = nbnodei_child.seqnum;
                            obj_child.nbnodelevel(end+1) = nbnodei_child.level;
                        end
                    end
                end
            end
        end
        
        % Recursively setNbNode.
        for iter = [1,2]
            obj.children{iter} = SetNeighborNode(obj.children{iter});
        end
        
        end
        
        function obj = FillTree(obj,A)
        % FillTree Fill tree structure with A.
        
        % Sort vtx, sep, nb.
        obj.vtx = sort(obj.vtx);
        obj.sep = sort(obj.sep);
        obj.nb = sort(obj.nb);
        
        % We only fill the leaf nodes.
        if obj.endflag == 1
            % First, we set int = vtx - sep (only holds on leaf nodes).
            obj.int = setdiff(obj.vtx,obj.sep,'sorted');
            
            % Then, we set the corresponding A**.
            obj.AII = full(A(obj.int,obj.int));
            obj.ASI = full(A(obj.sep,obj.int));
            obj.ASS = full(A(obj.sep,obj.sep));
            obj.ANS = full(A(obj.nb,obj.sep));
            
            % Set sep type.
            obj = SetSeparatorType(obj);
        else
            for iter = [1,2]
                obj.children{iter} = FillTree(obj.children{iter},A);
            end
        end
        
        end
        
        function obj = RecursiveSetSeparatorType(obj,whatlevel)
        % RecursiveSetSeparatorType Recusively set separator type.
        
        if obj.level == whatlevel
            obj = SetSeparatorType(obj);
        else
            if obj.endflag == 0
                for iter = [1,2]
                    obj.children{iter} = RecursiveSetSeparatorType(obj.children{iter},whatlevel);
                end
            end
        end
        
        end
        
        function obj = SetSeparatorType(obj)
        % SetSeparatorType Set separator type.

        ordersep = zeros(length(obj.sep),1);
        
        for k = 1:length(obj.nbnode)
            nodek = obj.nbnode{k};
            [~,tmp,~] = intersect(obj.sep,nodek.nb);
            obj.singlesep{k} = tmp;
            ordersep(tmp) = ordersep(tmp) + 1;
        end
        
        obj.complexsep = obj.sep;
        for k = 1:length(obj.nbnode)
            obj.singlesep{k} = obj.sep(intersect(find(ordersep == 1),obj.singlesep{k}));
            obj.complexsep = setdiff(obj.complexsep,obj.singlesep{k});
        end
        
        % obj.complexsep = obj.sep(ordersep > 1);
        % There exists the case where ordersep == 0 is not empty, so this
        % is wrong.
        
        end
        
        function obj = Factorization(obj,tol,demoHIF)
        % Factorization HIF factorization.
        
        if nargin == 1
            disp(" ");
            disp(" Default tol :1e-3 ");
            disp(" ");
            tol = 1e-3;
            demoHIF = 0;
        end
        
        if nargin > 2
            obj.demoHIF = demoHIF;
        end
        
        disp("  ");
        disp(" Start factorization ");
        disp("  ");
        
        for tmplevel = obj.numlevels:-1:1
            % Sparse elimination.
            obj = RecursiveSparseElim(obj,tmplevel);
            % Demo the process.
            if obj.demoHIF == 1
                DemoHIF(obj,tmplevel);
            end
            % Skeletonization.
            obj = RecursiveSkel(obj,tmplevel,tol);
            % Demo the process.
            if ((obj.demoHIF == 1) && (tmplevel ~= obj.numlevels))
                DemoHIF(obj,tmplevel);
            end
            % Merge.
            obj = RecursiveMerge(obj,tmplevel-1);
            % SetSeparatorType.
            obj = RecursiveSetSeparatorType(obj,tmplevel-1);
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
            if obj.endflag == 0
                for iter = [1,2]
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
        
        % DEBUG: Check whether the int is decoupled with other vertices.
        obj.root.inputAxy.A(obj.int,obj.int) = 0;
        obj.root.inputAxy.A(obj.int,obj.sep) = 0;
        obj.root.inputAxy.A(obj.sep,obj.int) = 0;
        
        end
        
        function obj = RecursiveSkel(obj,whatlevel,tol)
        % RecursiveSkel Recursively skeletonization.
        
        if obj.level == whatlevel
            obj = Skel(obj,tol);
            % obj = NoSkel(obj); % No skeletonization.
        else
            if obj.endflag == 0
                for iter = [1,2]
                    obj.children{iter} = RecursiveSkel(obj.children{iter},whatlevel,tol);
                end
            end
        end
        
        end
        
        function obj = Skel(obj,tol)
        % Skel Skeletonization.
        
        for k = 1:length(obj.nbnode)
            % We do skel according to nbnode.
            nodek = obj.nbnode{k};
            if nodek.level ~= obj.level || nodek.copy == 1
                obj.nbinfo(k).empty = 1;
                continue;
            end
            if nodek.seqnum < obj.seqnum
                obj.nbinfo(k).empty = 1;
                continue;
            end
            
            % The following data are vertices.
            sep1 = obj.singlesep{k};
            mysep1C = []; % sep - sep1.
            for nok = 1:length(obj.nbnode)
                if nok == k
                    continue;
                else
                    mysep1C = [mysep1C,obj.singlesep{nok}];
                end
            end
            mysep1C = [mysep1C,obj.complexsep];
            mysep1C = sort(mysep1C);
            nodeksep1C = setdiff(nodek.nb,sep1,'sorted');
                    
            % debug:
            nodeksep1C = intersect(nodeksep1C,obj.sep,'sorted');
            
            korder = find(nodek.nbnodeseqnum == obj.seqnum);
            if length(korder) > 1
                klevel = find(nodek.nbnodelevel == obj.level);
                korder = intersect(klevel,korder);
            end
            sep2 = nodek.singlesep{korder};
            nodeksep2C = [];
            for nok = 1:length(nodek.nbnode)
                if nok == korder
                    continue;
                else
                    nodeksep2C = [nodeksep2C,nodek.singlesep{nok}];
                end
            end
            nodeksep2C = [nodeksep2C,nodek.complexsep];
            nodeksep2C = sort(nodeksep2C);
            mysep2C = setdiff(obj.nb,sep2,'sorted');
            
            % debug:
            mysep2C = intersect(mysep2C,nodek.sep,'sorted');
            
            % The following data are indices.
            [~,myindex_sep1] = ismember(sep1,obj.sep);
            [~,nodekindex_sep1] = ismember(sep1,nodek.nb);
            [~,myindex_sep1C] = ismember(mysep1C,obj.sep);
            [~,myindex_mysep2C] = ismember(mysep2C,obj.nb);
            
            [~,nodekindex_mysep2C] = ismember(mysep2C,nodek.sep);
            
            [~,nodekindex_sep2] = ismember(sep2,nodek.sep);
            [~,myindex_sep2] = ismember(sep2,obj.nb);
            [~,nodekindex_sep2C] = ismember(nodeksep2C,nodek.sep);
            [~,nodekindex_nodeksep1C] = ismember(nodeksep1C,nodek.nb);
            
            [~,myindex_nodeksep1C] = ismember(nodeksep1C,obj.sep);
            
            % ID decomposition.
            
            % In the following process, the first "1" or "2" denotes my or
            % nodek, the second "1" or "2" denotes sk or re.
            
            skelmtx1 = [obj.ASS(myindex_sep1C,myindex_sep1);
                obj.ANS(myindex_mysep2C,myindex_sep1)];
            skelmtx2 = [nodek.ASS(nodekindex_sep2C,nodekindex_sep2);
                nodek.ANS(nodekindex_nodeksep1C,nodekindex_sep2)];
            if isempty(skelmtx1) || isempty(skelmtx2)
                obj.nbinfo(k).empty = 1;
                continue;
            end
            
            [T1,p11,p12] = ID(skelmtx1,tol); % skelmtx1(:,p12) = skelmtx1(:,p11) * T1.
            
            myindex_p11 = myindex_sep1(p11);
            myindex_p12 = myindex_sep1(p12);
            nodekindex_p11 = nodekindex_sep1(p11);
            nodekindex_p12 = nodekindex_sep1(p12);
            
            obj.re = [obj.re,sep1(p12)];
            nodek.nbre = [nodek.nbre,sep1(p12)];
            obj.nbinfo(k).Th1c1 = T1;
            obj.nbinfo(k).myindex_p11 = myindex_p11;
            obj.nbinfo(k).myindex_p12 = myindex_p12;
            obj.nbinfo(k).nodekindex_p11 = nodekindex_p11;
            obj.nbinfo(k).nodekindex_p12 = nodekindex_p12;
            
            [T2,p21,p22] = ID(skelmtx2,tol); % skelmtx2(:,p22) = skelmtx1(:,p21) * T2.
            
            myindex_p21 = myindex_sep2(p21);
            myindex_p22 = myindex_sep2(p22);
            nodekindex_p21 = nodekindex_sep2(p21);
            nodekindex_p22 = nodekindex_sep2(p22);
            
            nodek.re = [nodek.re,sep2(p22)];
            obj.nbre = [obj.nbre,sep2(p22)];
            obj.nbinfo(k).Th2c2 = T2;
            obj.nbinfo(k).myindex_p21 = myindex_p21;
            obj.nbinfo(k).myindex_p22 = myindex_p22;
            obj.nbinfo(k).nodekindex_p21 = nodekindex_p21;
            obj.nbinfo(k).nodekindex_p22 = nodekindex_p22;
            
            obj.nbinfo(k).empty = 0;
            
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
%             % Ad1c1 = Ac1d1 = 0, At1c1 = Ac1t1 = 0.
%             obj.ASS(myindex_sep1C,myindex_p12) = 0;
%             obj.ASS(myindex_p12,myindex_sep1C) = 0;
%             % Ad2c1 = Ac1d2 = 0.
%             obj.ANS(myindex_mysep2C,myindex_p12) = 0;
%             nodek.ANS(nodekindex_p12,nodekindex_mysep2C) = 0;
%             % Ad1c2 = Ac2d1 = 0, At2c2 = Ac2t2 = 0;
%             nodek.ASS(nodekindex_sep2C,nodekindex_p22) = 0;
%             nodek.ASS(nodekindex_p22,nodekindex_sep2C) = 0;
%             % Ad1c2 = Ac2d1 = 0.
%             nodek.ANS(nodekindex_nodeksep1C,nodekindex_p22) = 0;
%             obj.ANS(myindex_p12,myindex_nodeksep1C) = 0;
            
            % DEBUG: Check whether the int is decoupled with other vertices.
            sep1_p11 = sep1(p11);
            sep1_p12 = sep1(p12);
            sep2_p21 = sep2(p21);
            sep2_p22 = sep2(p22);
            obj.root.inputAxy.A(mysep1C,sep1_p12) = 0;
            obj.root.inputAxy.A(sep1_p12,mysep1C) = 0;
            obj.root.inputAxy.A(mysep2C,sep1_p12) = 0;
            obj.root.inputAxy.A(sep1_p12,mysep2C) = 0;
            obj.root.inputAxy.A(nodeksep2C,sep2_p22) = 0;
            obj.root.inputAxy.A(sep2_p22,nodeksep2C) = 0;
            obj.root.inputAxy.A(nodeksep1C,sep2_p22) = 0;
            obj.root.inputAxy.A(sep2_p22,nodeksep1C) = 0;
        
            % Step 2
            % Ac1c1 = Lc1 * Dc1 * Lc1^{T}.
            [L1,D1] = ldl(obj.ASS(myindex_p12,myindex_p12));      
            obj.root.active(sep1(p12)) = 0;
            obj.nbinfo(k).Lc1 = L1;
            obj.nbinfo(k).Dc1 = D1;
            % Ac1c1invAc1h1 = Ac1c1^{-1} * Ah1c1^{T}.
            obj.nbinfo(k).Ac1c1invAc1h1 = L1\(obj.ASS(myindex_p11,myindex_p12)');
            obj.nbinfo(k).Ac1c1invAc1h1 = D1\obj.nbinfo(k).Ac1c1invAc1h1;
            obj.nbinfo(k).Ac1c1invAc1h1 = L1'\obj.nbinfo(k).Ac1c1invAc1h1;
            % Ac1c1invAc1c2 = Ac1c1^{-1} * Ac2c1^{T}.
            obj.nbinfo(k).Ac1c1invAc1c2 = L1\(obj.ANS(myindex_p22,myindex_p12)');
            obj.nbinfo(k).Ac1c1invAc1c2 = D1\obj.nbinfo(k).Ac1c1invAc1c2;
            obj.nbinfo(k).Ac1c1invAc1c2 = L1'\obj.nbinfo(k).Ac1c1invAc1c2;
            % Ac1c1invAc1h2 = Ac1c1^{-1} * Ah2c1^{T}.
            obj.nbinfo(k).Ac1c1invAc1h2 = L1\(obj.ANS(myindex_p21,myindex_p12)');
            obj.nbinfo(k).Ac1c1invAc1h2 = D1\obj.nbinfo(k).Ac1c1invAc1h2;
            obj.nbinfo(k).Ac1c1invAc1h2 = L1'\obj.nbinfo(k).Ac1c1invAc1h2;
            % Ah1h1 = Ah1h1 - Ah1c1 * Ac1c1^{-1} * Ah1c1^{T}.
            obj.ASS(myindex_p11,myindex_p11) = obj.ASS(myindex_p11,myindex_p11) - obj.ASS(myindex_p11,myindex_p12)*obj.nbinfo(k).Ac1c1invAc1h1;
            % Ac2h1 = Ac2h1 - Ac2c1 * Ac1c1^{-1} * Ah1c1^{T}.
            obj.ANS(myindex_p22,myindex_p11) = obj.ANS(myindex_p22,myindex_p11) - obj.ANS(myindex_p22,myindex_p12)*obj.nbinfo(k).Ac1c1invAc1h1;
            nodek.ANS(nodekindex_p11,nodekindex_p22) = obj.ANS(myindex_p22,myindex_p11)';
            % Ah2h1 = Ah2h1 - Ah2c1 * Ac1c1^{-1} * Ah1c1^{T}.
            obj.ANS(myindex_p21,myindex_p11) = obj.ANS(myindex_p21,myindex_p11) - obj.ANS(myindex_p21,myindex_p12)*obj.nbinfo(k).Ac1c1invAc1h1;
            nodek.ANS(nodekindex_p11,nodekindex_p21) = obj.ANS(myindex_p21,myindex_p11)';
            % Ac2c2 = Ac2c2 - Ac2c1 * Ac1c1^{-1} * Ac2c1^{T}.
            nodek.ASS(nodekindex_p22,nodekindex_p22) = nodek.ASS(nodekindex_p22,nodekindex_p22) - obj.ANS(myindex_p22,myindex_p12)*obj.nbinfo(k).Ac1c1invAc1c2;
            % Ah2c2 = Ah2c2 - Ah2c1 * Ac1c1^{-1} * Ac2c1^{T}.
            nodek.ASS(nodekindex_p21,nodekindex_p22) = nodek.ASS(nodekindex_p21,nodekindex_p22) - obj.ANS(myindex_p21,myindex_p12)*obj.nbinfo(k).Ac1c1invAc1c2;
            nodek.ASS(nodekindex_p22,nodekindex_p21) = nodek.ASS(nodekindex_p21,nodekindex_p22)';
            % Ah2h2 = Ah2h2 - Ah2c1 * Ac1c1^{-1} * Ah2c1^{T}.
            nodek.ASS(nodekindex_p21,nodekindex_p21) = nodek.ASS(nodekindex_p21,nodekindex_p21) - obj.ANS(myindex_p21,myindex_p12)*obj.nbinfo(k).Ac1c1invAc1h2;
            % Ah1c1 = Ac2c1 = Ah2c1 = 0.
%             obj.ASS(myindex_p11,myindex_p12) = 0;
%             obj.ASS(myindex_p12,myindex_p11) = 0;
%             obj.ANS(myindex_p22,myindex_p12) = 0;
%             nodek.ANS(nodekindex_p12,nodekindex_p22) = 0;
%             obj.ANS(myindex_p21,myindex_p12) = 0;
%             nodek.ANS(nodekindex_p12,nodekindex_p21) = 0;
            
            % DEBUG: Check whether the int is decoupled with other vertices.
            obj.root.inputAxy.A(sep1_p12,sep1_p12) = 0;
            obj.root.inputAxy.A(sep1_p11,sep1_p12) = 0;
            obj.root.inputAxy.A(sep1_p12,sep1_p11) = 0;
            obj.root.inputAxy.A(sep2,sep1_p12) = 0;
            obj.root.inputAxy.A(sep1_p12,sep2) = 0;
            
            % Step 3
            % Ac2c2 = Lc2 * Dc2 * Lc2^{T};
            [L2,D2] = ldl(nodek.ASS(nodekindex_p22,nodekindex_p22));
            obj.root.active(sep2(p22)) = 0;
            obj.nbinfo(k).Lc2 = L2;
            obj.nbinfo(k).Dc2 = D2;
            % Ac2c2invAc2h1 = Ac2c2^{-1} * Ac2h1.
            obj.nbinfo(k).Ac2c2invAc2h1 = L2\obj.ANS(myindex_p22,myindex_p11);
            obj.nbinfo(k).Ac2c2invAc2h1 = D2\obj.nbinfo(k).Ac2c2invAc2h1;
            obj.nbinfo(k).Ac2c2invAc2h1 = L2'\obj.nbinfo(k).Ac2c2invAc2h1;
            % Ac2c2invAc2h2 = Ac2c2^{-1} * Ah2c2^{T}.
            obj.nbinfo(k).Ac2c2invAc2h2 = L2\(nodek.ASS(nodekindex_p21,nodekindex_p22)');
            obj.nbinfo(k).Ac2c2invAc2h2 = D2\obj.nbinfo(k).Ac2c2invAc2h2;
            obj.nbinfo(k).Ac2c2invAc2h2 = L2'\obj.nbinfo(k).Ac2c2invAc2h2;
            % Ah1h1 = Ah1h1 - Ac2h1^{T} * Ac2c2^{-1} * Ac2h1.
            obj.ASS(myindex_p11,myindex_p11) = obj.ASS(myindex_p11,myindex_p11) - obj.ANS(myindex_p22,myindex_p11)'*obj.nbinfo(k).Ac2c2invAc2h1;
            % Ah2h1 = Ah2h1 - Ah2c2 * Ac2c2^{-1} * Ac2h1.
            obj.ANS(myindex_p21,myindex_p11) = obj.ANS(myindex_p21,myindex_p11) - nodek.ASS(nodekindex_p21,nodekindex_p22)*obj.nbinfo(k).Ac2c2invAc2h1;
            nodek.ANS(nodekindex_p11,nodekindex_p21) = obj.ANS(myindex_p21,myindex_p11)';
            % Ah2h2 = Ah2h2 - Ah2c2 * Ac2c2^{-1} * Ah2c2^{T}.
            nodek.ASS(nodekindex_p21,nodekindex_p21) = nodek.ASS(nodekindex_p21,nodekindex_p21)- nodek.ASS(nodekindex_p21,nodekindex_p22)*obj.nbinfo(k).Ac2c2invAc2h2;
            % Ah2c2 = Ac2h1 = 0.
%             nodek.ASS(nodekindex_p21,nodekindex_p22) = 0;
%             nodek.ASS(nodekindex_p22,nodekindex_p21) = 0;
%             nodek.ANS(nodekindex_p11,nodekindex_p22) = 0;
%             obj.ANS(myindex_p22,myindex_p11) = 0;
            
%             obj.ASS(myindex_p12,:) = 0;
%             obj.ASS(:,myindex_p12) = 0;
%             obj.ANS(:,myindex_p12) = 0;
%             nodek.ANS(nodekindex_p12,:) = 0;
%             
%             nodek.ASS(nodekindex_p22,:) = 0;
%             nodek.ASS(:,nodekindex_p22) = 0;
%             nodek.ANS(:,nodekindex_p22) = 0;
%             obj.ANS(myindex_p22,:) = 0;
            
            % DEBUG: Check whether the int is decoupled with other vertices.
            obj.root.inputAxy.A(sep2_p22,sep2_p22) = 0;
            obj.root.inputAxy.A(sep2_p21,sep2_p22) = 0;
            obj.root.inputAxy.A(sep2_p22,sep2_p21) = 0;
            obj.root.inputAxy.A(sep1_p11,sep2_p22) = 0;
            obj.root.inputAxy.A(sep2_p22,sep1_p11) = 0;
        end
        
        obj.re = sort(obj.re);
        obj.sk = setdiff(obj.sep,obj.re,'sorted');
        obj.nbre = sort(obj.nbre);
        obj.nbsk = setdiff(obj.nb,obj.nbre,'sorted');
        
        end
        
        function obj = NoSkel(obj)
        % NoSkel No skeletonization.
        
        obj.sk = obj.sep;
        obj.nbsk = obj.nb;
        
        end
        
        function obj = RecursiveMerge(obj,whatlevel)
        % RecursiveMerge Recusively send matrices' information from children to parent.
        
        if obj.level == whatlevel
            obj = Merge(obj);
        else
            if obj.endflag == 0
                for iter = [1,2]
                    obj.children{iter} = RecursiveMerge(obj.children{iter},whatlevel);
                end
            end
        end
        
        end
        
        function obj = Merge(obj)
        % Merge Send matrices' information from children to parent.
        
        % We stand on the parent level.
        
        if obj.endflag == 1
            return;
        end
        
        % First we tell the parent what its int, sep, nb is after we
        % eliminate the children's vtx.
        % int: children's sk - sep.
        % sep: sep \cup children's sk.
        % nb: nb \cup children's nbsk.
        tmp = [];
        for iter = [1,2]
            obj.int = [obj.int,obj.children{iter}.sk];
            tmp = [tmp,obj.children{iter}.nbsk];
        end
        obj.sep = intersect(obj.sep,obj.int,'sorted');
        obj.int = setdiff(obj.int,obj.sep,'sorted');
        obj.nb = intersect(obj.nb,tmp,'sorted');
        
        % Next we assign the corresponding matrices.
        % From child to parent.
        % int: children's sep and children's nb
        % sep: children's sep and children's nb
        % nb: children's nb
        
        % We assign values blockly.
        
        % AII
        % An int of the parent only belongs to the sep of one of its
        % children. If two ints belong to the same child, we assign AII
        % from the child's ASS. Otherwise, we assign AII from one child's
        % ANS or 0.
        obj.AII = zeros(length(obj.int));
        [~,myindex_int1,cindex_int1] = intersect(obj.int,obj.children{1}.sep);
        [~,myindex_int2,cindex_int2] = intersect(obj.int,obj.children{2}.sep);
        [~,myindex_int21,cindex_int21] = intersect(obj.int,obj.children{1}.nb);
        obj.indexinfo(1).myindex_int = myindex_int1;
        obj.indexinfo(1).cindex_int = cindex_int1;
        obj.indexinfo(2).myindex_int = myindex_int2;
        obj.indexinfo(2).cindex_int = cindex_int2;
        obj.AII(myindex_int1,myindex_int1) = obj.children{1}.ASS(cindex_int1,cindex_int1);
        obj.AII(myindex_int2,myindex_int2) = obj.children{2}.ASS(cindex_int2,cindex_int2);
        obj.AII(myindex_int21,myindex_int1) = obj.children{1}.ANS(cindex_int21,cindex_int1);
        obj.AII(myindex_int1,myindex_int2) = obj.AII(myindex_int2,myindex_int1)';
        
        % ASI
        % A sep of the parent only belongs to the sep of one of its
        % children. If an int and a sep belongs to the same child, we
        % assign ASI from the child's ASS. Otherwise, we assign ASI from
        % one child's ANS or 0.
        obj.ASI = zeros(length(obj.sep),length(obj.int));
        [~,myindex_sep1x,cindex_sep1x] = intersect(obj.sep,obj.children{1}.sep);
        [~,myindex_sep1y,cindex_sep1y] = intersect(obj.sep,obj.children{1}.nb);
        [~,myindex_sep2x,cindex_sep2x] = intersect(obj.sep,obj.children{2}.sep);
        [~,myindex_sep2y,cindex_sep2y] = intersect(obj.sep,obj.children{2}.nb);
        obj.ASI(myindex_sep1x,myindex_int1) = obj.children{1}.ASS(cindex_sep1x,cindex_int1);
        obj.ASI(myindex_sep2x,myindex_int2) = obj.children{2}.ASS(cindex_sep2x,cindex_int2);
        obj.ASI(myindex_sep1y,myindex_int1) = obj.children{1}.ANS(cindex_sep1y,cindex_int1);
        obj.ASI(myindex_sep2y,myindex_int2) = obj.children{2}.ANS(cindex_sep2y,cindex_int2);
        
        % ASS
        % If two seps belongs to the same child, we assign ASS from the
        % child's ASS. Otherwise, we assign ASS from one child's ANS or 0.
        obj.ASS = zeros(length(obj.sep));
        [~,myindex_sep1,cindex_sep1] = intersect(obj.sep,obj.children{1}.sep);
        [~,myindex_sep2,cindex_sep2] = intersect(obj.sep,obj.children{2}.sep);
        [~,myindex_sep21,cindex_sep21] = intersect(obj.sep,obj.children{1}.nb);
        obj.indexinfo(1).myindex_sep = myindex_sep1;
        obj.indexinfo(1).cindex_sep = cindex_sep1;
        obj.indexinfo(2).myindex_sep = myindex_sep2;
        obj.indexinfo(2).cindex_sep = cindex_sep2;
        obj.ASS(myindex_sep1,myindex_sep1) = obj.children{1}.ASS(cindex_sep1,cindex_sep1);
        obj.ASS(myindex_sep2,myindex_sep2) = obj.children{2}.ASS(cindex_sep2,cindex_sep2);
        obj.ASS(myindex_sep21,myindex_sep1) = obj.children{1}.ANS(cindex_sep21,cindex_sep1);
        obj.ASS(myindex_sep1,myindex_sep2) = obj.ASS(myindex_sep2,myindex_sep1)';
        
        % ANS.
        % If a nb and a sep in the same child, we assign ANS from the
        % child's ANS. Otherwise, ANS= 0.
        obj.ANS = zeros(length(obj.nb),length(obj.sep));
        [~,myindex_nb1x,cindex_nb1x] = intersect(obj.nb,obj.children{1}.nb);
        [~,myindex_nb2x,cindex_nb2x] = intersect(obj.nb,obj.children{2}.nb);
        obj.ANS(myindex_nb1x,myindex_sep1) = obj.children{1}.ANS(cindex_nb1x,cindex_sep1);
        obj.ANS(myindex_nb2x,myindex_sep2) = obj.children{2}.ANS(cindex_nb2x,cindex_sep2);
        
        % Clear children information.
        for iter = [1,2]
            obj.children{iter} = HIFClear(obj.children{iter});
        end
        
        end
        
        function obj = HIFClear(obj)
        % HIFClear Clear unnecessary information.
        
        obj.AII = [];
        obj.ASI = [];
        obj.ASS = [];
        obj.ANS = [];
        
        end
        
        function obj = RootFactorization(obj)
        % RootFactorization Factorization on the root.
        
        obj.root.active(obj.int) = 0;
        % AII = LI * DI * LI^{T}.
        [obj.LI,obj.DI] = ldl(obj.AII);
        
        % DEBUG: Check whether the int is decoupled with other vertices.
        obj.root.inputAxy.A(obj.int,obj.int) = 0;
        
        end
        
        function x = HIFSolve(obj,b)
        % HIFSolve Solve Ax = b through HIF.
        
        disp("  ");
        disp(" Start solve ");
        disp("  ");
        
        obj = BuildVecTree(obj,b);
        
        for tmplevel = obj.numlevels:-1:1
            obj = RecursiveApplySparseElimUp(obj,tmplevel);
            obj = RecursiveApplySkelUp(obj,tmplevel);
            obj = RecursiveApplyMerge(obj,tmplevel-1);
        end
        
        obj = RootApply(obj);
        
        for tmplevel = 1:1:obj.numlevels
            obj = RecursiveApplySplit(obj,tmplevel-1);
            obj = RecursiveApplySkelDown(obj,tmplevel);
            obj = RecursiveApplySparseElimDown(obj,tmplevel);
        end
        
        x = zeros(size(b));
        x = GetSolution(obj,x);
        
        disp("  ");
        disp(" End solve ");
        disp("  ");
        
        end
        
        function obj = BuildVecTree(obj,b)
        % BuildVecTree Fill b in our HIF tree.
        
        if obj.level == 0
            disp("  ");
            disp(" Start build vector tree ");
            disp("  ");
        end
        
        if obj.endflag == 0
            for iter = [1,2]
                obj.children{iter} = BuildVecTree(obj.children{iter},b);
            end
        else
            obj.xI = b(obj.int,:);
            obj.xS = b(obj.sep,:);
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
            if obj.endflag == 0
                for iter = [1,2]
                    obj.children{iter} = RecursiveApplySparseElimUp(obj.children{iter},whatlevel);
                end
            end
        end
        
        end
        
        function obj = ApplySparseElimUp(obj)
        % ApplySparseElimUp Phase 1 for applying sparse elimination.
        
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
            if obj.endflag == 0
                for iter = [1,2]
                    obj.children{iter} = RecursiveApplySkelUp(obj.children{iter},whatlevel);
                end
            end
        end
        
        end
        
        function obj = ApplySkelUp(obj)
        % ApplySkelUp Phase 1 for applying skeletonization.
        
        for k = 1:length(obj.nbnode)
            nbnodek = obj.nbnode{k};
            if isempty(obj.nbinfo)
                continue;
            end
            nbinfok = obj.nbinfo(k);
            if nbinfok.empty
                continue;
            end
            
            
            % Step 1
            % xc1 = xc1 - Th1c1^{T} * xh1.
            obj.xS(nbinfok.myindex_p12,:) = obj.xS(nbinfok.myindex_p12,:) - nbinfok.Th1c1'*obj.xS(nbinfok.myindex_p11,:);
            % xc2 = xc2 - Th2c2^{T} * xh2.
            nbnodek.xS(nbinfok.nodekindex_p22,:) = nbnodek.xS(nbinfok.nodekindex_p22,:) - nbinfok.Th2c2'*nbnodek.xS(nbinfok.nodekindex_p21,:);
            
            % Step 2
            % xh1 = xh1 - (Ac1c1^{-1} * Ah1c1^{T})^{T} * xc1.
            obj.xS(nbinfok.myindex_p11,:) = obj.xS(nbinfok.myindex_p11,:) - nbinfok.Ac1c1invAc1h1'*obj.xS(nbinfok.myindex_p12,:);
            % xc2 = xc2 - (Ac1c1^{-1} * Ac2c1^{T})^{T} * xc1.
            nbnodek.xS(nbinfok.nodekindex_p22,:) = nbnodek.xS(nbinfok.nodekindex_p22,:) - nbinfok.Ac1c1invAc1c2'*obj.xS(nbinfok.myindex_p12,:);
            % xh2 = xh2 - (Ac1c1^{-1} * Ah2c1^{T})^{T} * xc1.
            nbnodek.xS(nbinfok.nodekindex_p21,:) = nbnodek.xS(nbinfok.nodekindex_p21,:) - nbinfok.Ac1c1invAc1h2'*obj.xS(nbinfok.myindex_p12,:);
            % xc1 = Lc1^{-1} * xc1.
            obj.xS(nbinfok.myindex_p12,:) = nbinfok.Lc1\obj.xS(nbinfok.myindex_p12,:);
            
            % Step 3
            % xh1 = xh1 - (Ac2c2^{-1} * Ac2h1)^{T} * xc2.
            obj.xS(nbinfok.myindex_p11,:) = obj.xS(nbinfok.myindex_p11,:) - nbinfok.Ac2c2invAc2h1'*nbnodek.xS(nbinfok.nodekindex_p22,:);
            % xh2 = xh2 - (Ac2c2^{-1} * Ah2c2^{T})^{T} * xc2.
            nbnodek.xS(nbinfok.nodekindex_p21,:) = nbnodek.xS(nbinfok.nodekindex_p21,:) - nbinfok.Ac2c2invAc2h2'*nbnodek.xS(nbinfok.nodekindex_p22,:);
            % xc2 = Lc2^{-1} * xc2.
            nbnodek.xS(nbinfok.nodekindex_p22,:) = nbinfok.Lc2\nbnodek.xS(nbinfok.nodekindex_p22,:);
            
            % xc1 = Dc1^{-1} * xc1. xc2 = Dc2^{-1} * xc2. We only apply D once.
            obj.xS(nbinfok.myindex_p12,:) = nbinfok.Dc1\obj.xS(nbinfok.myindex_p12,:);
            nbnodek.xS(nbinfok.nodekindex_p22,:) = nbinfok.Dc2\nbnodek.xS(nbinfok.nodekindex_p22,:);
        end
        
        end
        
        function obj = RecursiveApplyMerge(obj,whatlevel)
        % RecursiveApplyMerge Recusively send vectors' information from children to parent.
        
        if obj.level == whatlevel
            obj = ApplyMerge(obj);
        else
            if obj.endflag == 0
                for iter = [1,2]
                    obj.children{iter} = RecursiveApplyMerge(obj.children{iter},whatlevel);
                end
            end
        end
        
        end
        
        function obj = ApplyMerge(obj)
        % ApplyMerge Send vectors' information from children to parent.
        
        % We stand on the parent level.
        
        if obj.endflag == 1
            return;
        end
        
        % We only need to assign the corresponding vectors.
        width = size(obj.children{1}.xS,2);
        
        % xI.
        % An int of the parent only belongs to the sep of one of its
        % children. We assign xI from the child's xS.
        obj.xI = zeros(length(obj.int),width);
        for iter = [1,2]
            obj.xI(obj.indexinfo(iter).myindex_int,:) = obj.children{iter}.xS(obj.indexinfo(iter).cindex_int,:);
        end
        
        % xS.
        % A sep of the parent only belongs to the sep of one of its
        % children. We assign xS from the child's xS.
        obj.xS = zeros(length(obj.sep),width);
        for iter = [1,2]
            obj.xS(obj.indexinfo(iter).myindex_sep,:) = obj.children{iter}.xS(obj.indexinfo(iter).cindex_sep,:);
        end
        
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
            if obj.endflag == 0
                for iter = [1,2]
                    obj.children{iter} = RecursiveApplySplit(obj.children{iter},whatlevel);
                end
            end
        end
        
        end
        
        function obj = ApplySplit(obj)
        % ApplySplit Send vectors' information from parent to children.
        
        % We stand on the parent level.
        
        if obj.endflag == 1
            return;
        end
        
        % We only need to assign the corresponding vectors of the children.
        
        for iter = [1,2]
            obj.children{iter}.xS(obj.indexinfo(iter).cindex_int,:) = obj.xI(obj.indexinfo(iter).myindex_int,:);
            obj.children{iter}.xS(obj.indexinfo(iter).cindex_sep,:) = obj.xS(obj.indexinfo(iter).myindex_sep,:);
        end
        
        end
        
        function obj = RecursiveApplySkelDown(obj,whatlevel)
        % RecursiveApplySkelDown Phase 2 for applying skeletonization recursively.
        
        if obj.level == whatlevel
            obj = ApplySkelDown(obj);
        else
            if obj.endflag == 0
                for iter = [1,2]
                    obj.children{iter} = RecursiveApplySkelDown(obj.children{iter},whatlevel);
                end
            end
        end
        
        end
        
        function obj = ApplySkelDown(obj)
        % ApplySkelDown Phase 2 for applying skeletonization.
        
        for k = 1:length(obj.nbnode)
            nbnodek = obj.nbnode{k};
            if isempty(obj.nbinfo)
                continue;
            end
            nbinfok = obj.nbinfo(k);
            if nbinfok.empty
                continue;
            end
            
            % Step 3
            % xc2 = Lc2^{-T} * xc2 - (Ac2c2^{-1} * Ac2h1) * xh1 - (Ac2c2^{-1} * Ah2c2^{T}) * xh2.
            nbnodek.xS(nbinfok.nodekindex_p22,:) = nbinfok.Lc2'\nbnodek.xS(nbinfok.nodekindex_p22,:) - ...
                nbinfok.Ac2c2invAc2h1*obj.xS(nbinfok.myindex_p11,:) - ...
                nbinfok.Ac2c2invAc2h2*nbnodek.xS(nbinfok.nodekindex_p21,:);
            
            % Step 2
            % xc1 = Lc1^{-T} * xc1 - (Ac1c1^{-1} * Ah1c1^{T}) * xh1 - (Ac1c1^{-1} * Ac2c1^{T}) * xc2 - (Ac1c1^{-1} * Ah2c1^{T}) * xh2.
            obj.xS(nbinfok.myindex_p12,:) = nbinfok.Lc1'\obj.xS(nbinfok.myindex_p12,:) - ...
                nbinfok.Ac1c1invAc1h1*obj.xS(nbinfok.myindex_p11,:) - ...
                nbinfok.Ac1c1invAc1c2*nbnodek.xS(nbinfok.nodekindex_p22,:) - ...
                nbinfok.Ac1c1invAc1h2*nbnodek.xS(nbinfok.nodekindex_p21,:);
            
            % Step 1
            % xh1 = xh1 - Th1c1 * xc1.
            obj.xS(nbinfok.myindex_p11,:) = obj.xS(nbinfok.myindex_p11,:) - nbinfok.Th1c1*obj.xS(nbinfok.myindex_p12,:);
            % xh2 = xh2 - Th2c2 * xc2.
            nbnodek.xS(nbinfok.nodekindex_p21,:) = nbnodek.xS(nbinfok.nodekindex_p21,:) - nbinfok.Th2c2*nbnodek.xS(nbinfok.nodekindex_p22,:);
        end
        
        
        end
        
        function obj = RecursiveApplySparseElimDown(obj,whatlevel)
        % RecursiveApplySparseElimDown Phase 2 for applying sparse elimination recursively.
        
        if obj.level == whatlevel
            obj = ApplySparseElimDown(obj);
        else
            if obj.endflag == 0
                for iter = [1,2]
                    obj.children{iter} = RecursiveApplySparseElimDown(obj.children{iter},whatlevel);
                end
            end
        end
        
        end
        
        function obj = ApplySparseElimDown(obj)
        % ApplySparseElimDown Phase 2 for applying sparse elimination.
        
        % xI = LI^{-T} * xI - (AII^{-1} * ASI^{T}) * xS.
        obj.xI = obj.LI'\obj.xI - obj.AIIinvAIS*obj.xS;
        
        end
        
        function x = GetSolution(obj,x)
        % GetSolution Get the final soluion through the tree structure.
        
        if obj.endflag == 1
            x(obj.int,:) = obj.xI;
            x(obj.sep,:) = obj.xS;
            obj.xI = [];
            obj.xS = [];
        else
            obj.xI = [];
            obj.xS = [];
            for iter = [1,2]
                x = GetSolution(obj.children{iter},x);
            end
        end
        
        end
        
        function DemoPart(obj)
        % DemoPart Demo the partition process.
        
        assert(~isempty(obj.inputAxy.xy),"Need coordinates!")
        
        disp("  ");
        disp(" Start demo the partition ");
        disp("  ");
        
        fig = figure();
        clf reset;
        colordef(fig,'black');
        gplotg(obj.inputAxy.A,obj.inputAxy.xy);
        
        disp(" Hit space to continue ... ");
        disp("  ");
        
        for tmplevel = 0:obj.numlevels
            disp(" Current level: " + tmplevel);
            disp("  ");
            
            map = GetPartMap(obj,tmplevel);
            gplotmap(obj.inputAxy.A,obj.inputAxy.xy,map);
            if tmplevel ~= obj.numlevels
                disp(" Hit space to continue ... ");
            else
                disp(" Hit space to end ... ");
            end
            disp("  ");
            pause;
        end
        
        end
        
        function DemoLevelPart(obj,whatlevel)
        % DemoLevelPart Demo the specified level partition.
        
        assert(~isempty(obj.inputAxy.xy),"Need coordinates!")
        
        disp("  ");
        disp(" Current level: " + whatlevel);
        disp("  ");
        
        fig = figure();
        clf reset;
        colordef(fig,'black');
        gplotg(obj.inputAxy.A,obj.inputAxy.xy);
        
        disp(" Hit space to continue ... ");
        disp("  ");
        pause;
        map = GetPartMap(obj,whatlevel);
        gplotmap(obj.inputAxy.A,obj.inputAxy.xy,map);
        disp(" Hit space to end ... ");
        disp("  ");
        pause;
        
        end
        
        function map = GetPartMap(obj,whatlevel,map)
        % GetPartMap Get the map of the partition.
        
        if nargin == 2
            map = [];
        end
        
        if obj.level == 0
            n = size(obj.inputAxy.xy,1);
            map = zeros(1,n);
        end
        
        if obj.level == whatlevel || obj.endflag == 1
            map(obj.vtx) = max(map)+1;
            return;
        else
            for iter = [1,2]
                map = GetPartMap(obj.children{iter},whatlevel,map);
            end
        end
        
        end
        
        function DemoFinalPart(obj)
        % DemoFinalPart Demo of the final partition.
        
        assert(~isempty(obj.inputAxy.xy),"Need coordinates!")
        
        disp("  ");
        disp(" Start demo the final partition ");
        disp("  ");
        
        fig = figure();
        clf reset;
        colordef(fig,'black');
        gplotg(speye(size(obj.inputAxy.A)),obj.inputAxy.xy);
        
        map = GetPartMap(obj,obj.numlevels);
        parts = unique(map);
        nparts = length(parts);
        tmpA = obj.inputAxy.A;
        for i = 1:nparts
            for j = i + 1:nparts
                idi = find(map == i);
                idj = find(map == j);
                tmpA(idi,idj) = 0;
                tmpA(idj,idi) = 0;
            end
        end
        gplotmap(tmpA,obj.inputAxy.xy,map);
        disp(" Hit space to end ... ");
        disp("  ");
        pause;
        
        end
        
        function DemoHIF(obj,whatlevel)
        % DemoHIF Demo the HIF process.
        
        assert(~isempty(obj.inputAxy.xy), "Need coordinates!")
        
        if whatlevel == obj.numlevels
            DemoFinalPart(obj);
        end
        
        disp("  ");
        disp(" Current level: " + whatlevel);
        disp("  ");
        
        fig = figure();
        clf reset;
        colordef(fig,'black');
        % gplotg(speye(size(obj.inputAxy.A)),obj.inputAxy.xy);
        
        map = GetPartMap(obj,whatlevel);
        tmpactive = find(obj.active > 0);
        parts = unique(map);
        nparts = length(parts);
        tmpA = obj.inputAxy.A;
        for i = 1:nparts
            for j = i + 1:nparts
                idi = find(map == i);
                idj = find(map == j);
                tmpA(idi,idj) = 0;
                tmpA(idj,idi) = 0;
            end
        end
        tmpA = tmpA(tmpactive,tmpactive);
        tmpxy = obj.inputAxy.xy(tmpactive,:);
        tmpmap = map(1,tmpactive);
        gplotmap(tmpA,tmpxy,tmpmap);
        if size(obj.inputAxy.xy,2) == 3
            view(3);
        end
        if whatlevel ~= 0
            disp(" Hit space to continue ...");
        else
            disp(" Hit space to end ...");
        end
        disp("  ");
        pause;
        
        end
        
        function levelVec = ReadLevel(obj,levelVec)
        % ReadLevel Obtain level vectors.
        
        if nargin == 1
            levelVec = [];
        end
        
        if obj.level == 0
            levelVec = zeros(obj.numlevels,1);
        end
        
        if obj.endflag == 1
            levelVec(obj.level) = levelVec(obj.level) + 1;
            return;
        else
            for iter = [1,2]
                levelVec = ReadLevel(obj.children{iter},levelVec);
            end
        end
        
        end
        
    end
    
end
