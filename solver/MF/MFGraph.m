classdef MFGraph < handle
    % MFGraph MF algorithm based on graph partition.
    
    properties
        
        % Root properties.
        
        % The following information will be stored only in the root node.
        
        inputAxy; % Input Axy.
        active; % Whether a vertex is eliminated.
        demoMF = 0; % Whether to demo the MF process.
        
        % Graph properties.
        
        vtx; % Vertices on the graph.
        sep; % Separator vertices.
        nb; % Neighbor vertices.
        int; % Interior vertices.
        
        % Tree properties.
        
        numlevels; % Total number of levels.
        level; % Current level, start from 0.
        seqnum; % A node's order in its level.
        endflag = 0; % Whether the partition ends.
        children = cell(1,2); % Children nodes.
        nbnode = {}; % Neighbor nodes. In fact, we don't need this in MF.
        nbnodeseqnum = []; % Neighbor nodes's seqnum.
        nbnodelevel = []; % Neighbor nodes's level.
        root; % Root node.
        indexinfo = struct([]); % Index information when merge and split.
        
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
        
        function obj = MFGraph(Axy,minvtx,method,level,seqnum,vtx,sep,nb)
        % MFGraph Create a MF class.
        
        if nargin == 1
            minvtx = 16;
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
        sep1 = unique(sep1);
        sep2 = unique(sep2);
        p = {p1,p2};
        partsep = {sep1,sep2};
        
        % Create children MF.
        for iter = [1,2]
            obj.children{iter} = MFGraph(Axy,minvtx,method,obj.level+1,obj.seqnum*2+(iter-1),...
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
                if isempty(index_sepi_childsep)
                    obj_child.sep = [obj_child.sep,sepi];
                end
                index_addnb_nb = find(nbA(i,:)~=0);% index_addnb_nb is always nonempty!
                for j = 1: length(index_addnb_nb)
                    addnbj = obj.nb(index_addnb_nb(j)); % addnbj is a neighbour vtx.
                    index_addnbj_childnb = find(obj_child.nb == addnbj,1);
                    % index_addnb_childnb may be nonempty (addnbj has been added).
                    if isempty(index_addnbj_childnb)
                        % Now addnbj is NOT in child's nb, we need to add nb.
                        obj_child.nb = [obj_child.nb,addnbj];
                    end
                end
            end
        end
        
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
            if ~isempty(obj.nbnode)
                for i = 1:length(obj.nbnode)
                    nbnodei = obj.nbnode{i};
                    % What we need is to check whether the vtx of nbnodei's
                    % chilren is in the nb of obj_child.
                    for it = [1,2]
                        nbnodei_child = nbnodei.children{it};
                        if isempty(nbnodei_child)
                            % The nbnodei doesn't have a child. We should
                            % look it as a nbnode.
                            if ~isempty(intersect(obj_child.nb,nbnodei.vtx))
                                % NOTE: We have to avoid add one's ancestor as its nbnode.
                                dlevel = obj_child.level - nbnodei.level;
                                myseqnum = obj_child.seqnum;
                                for k = 1:dlevel
                                    myseqnum = floor(myseqnum/2);
                                end
                                if myseqnum == nbnodei.seqnum
                                    break;
                                end
                                
                                obj_child.nbnode{end+1} = nbnodei;
                                obj_child.nbnodeseqnum(end+1) = nbnodei.seqnum;
                                obj_child.nbnodelevel(end+1) = nbnodei.level;
                                nbnodei.nbnode{end+1} = obj_child;
                                nbnodei.nbnodeseqnum(end+1) = obj_child.seqnum;
                                nbnodei.nbnodelevel(end+1) = obj_child.level;
                            end
                            break;
                        else
                            if ~isempty(intersect(obj_child.nb, nbnodei_child.vtx))
                                obj_child.nbnode{end+1} = nbnodei_child;
                                obj_child.nbnodeseqnum(end+1) = nbnodei_child.seqnum;
                                obj_child.nbnodelevel(end+1) = nbnodei_child.level;
                            end
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
        
        if obj.endflag == 0
            % We only fill the leaf nodes.
            for iter = [1,2]
                obj.children{iter} = FillTree(obj.children{iter},A);
            end
        else
            % First, we set int = vtx - sep (only holds on leaf nodes).
            obj.int = setdiff(obj.vtx,obj.sep,'sorted');
            % Then, we set the corresponding A**.
            obj.AII = A(obj.int,obj.int);
            obj.ASI = A(obj.sep,obj.int);
            obj.ASS = A(obj.sep,obj.sep);
            obj.ANS = A(obj.nb,obj.sep);
        end
        
        end
        
        function obj = Factorization(obj,demoMF)
        % Factorization MF factorization.
        
        if nargin > 1
            obj.demoMF = demoMF;
        end
        
        disp("  ");
        disp(" Start factorization ");
        disp("  ");
        
        for tmplevel = obj.numlevels:-1:1
            % Sparse elimination.
            obj = RecursiveSparseElim(obj,tmplevel);
            % Demo the process.
            if obj.demoMF == 1
                DemoMF(obj,tmplevel);
            end
            % Merge.
            obj = RecursiveMerge(obj,tmplevel-1);
        end
        
        % Root factorization.
        obj = RootFactorization(obj);
        % Demo the process.
        if obj.demoMF == 1
            DemoMF(obj,0);
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
        % ASI = 0.
        
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
        
        % First we tell the parent what its int is after we eliminate
        % the children's vtx.
        % int: children's sep - sep.
        for iter = [1,2]
            obj.int = [obj.int,obj.children{iter}.sep];
        end
        obj.int = setdiff(obj.int,obj.sep,'sorted');
        
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
        [int1,myindex_int1,~] = intersect(obj.int,obj.children{1}.vtx);
        [~,cindex_int1] = ismember(int1,obj.children{1}.sep);
        [int2,myindex_int2,~] = intersect(obj.int,obj.children{2}.vtx);
        [~,cindex_int2] = ismember(int2,obj.children{2}.sep);
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
        [sep1,myindex_sep1,~] = intersect(obj.sep,obj.children{1}.vtx);
        [~,cindex_sep1] = ismember(sep1,obj.children{1}.sep);
        [sep2,myindex_sep2,~] = intersect(obj.sep,obj.children{2}.vtx);
        [~,cindex_sep2] = ismember(sep2,obj.children{2}.sep);
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
        % child's ANS. Otherwise, ANS = 0.
        obj.ANS = zeros(length(obj.nb),length(obj.sep));
        [~,myindex_nb1x,cindex_nb1x] = intersect(obj.nb,obj.children{1}.nb);
        [~,myindex_nb2x,cindex_nb2x] = intersect(obj.nb,obj.children{2}.nb);
        obj.ANS(myindex_nb1x,myindex_sep1) = obj.children{1}.ANS(cindex_nb1x,cindex_sep1);
        obj.ANS(myindex_nb2x,myindex_sep2) = obj.children{2}.ANS(cindex_nb2x,cindex_sep2);
        
        % Clear children information.
        for iter = [1,2]
            obj.children{iter} = MFClear(obj.children{iter});
        end
        
        end
        
        function obj = MFClear(obj)
        % HIFClear Clear unnecessary information.
        
        % TODO: MFClear.
        
        end
        
        function obj = RootFactorization(obj)
        % RootFactorization Factorization on the root.
        
        obj.root.active(obj.int) = 0;
        % AII = LI * DI * LI^{T}.
        [obj.LI,obj.DI] = ldl(obj.AII);
        
        end
        
        function x = MFSolve(obj,b)
        % MFSolve Solve Ax = b through MF.
        
        disp("  ");
        disp(" Start solve ");
        disp("  ");
        
        obj = BuildVecTree(obj,b);
        
        for tmplevel = obj.numlevels:-1:1
            obj = RecursiveApplySparseElimUp(obj,tmplevel);
            obj = RecursiveApplyMerge(obj,tmplevel-1);
        end
        
        obj = RootApply(obj);
        
        for tmplevel = 1:1:obj.numlevels
            obj = RecursiveApplySplit(obj,tmplevel-1);
            obj = RecursiveApplySparseElimDown(obj,tmplevel);
        end
        
        x = zeros(size(b));
        x = GetSolution(obj,x);
        
        disp("  ");
        disp(" End solve ");
        disp("  ");
        
        end
        
        function obj = BuildVecTree(obj,b)
        % BuildVecTree Fill b in our MF tree.
        
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
        
        if obj.endflag == 0
            for iter = [1,2]
                x = GetSolution(obj.children{iter},x);
            end
        else
            x(obj.int,:) = obj.xI;
            x(obj.sep,:) = obj.xS;
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
        gplotg(obj.inputAxy.A,obj.inputAxy.xy);
        
        disp(" Hit space to continue ... ");
        disp("  ");
        pause;
        map = GetPartMap(obj,obj.numlevels);
        gplotmap(obj.inputAxy.A,obj.inputAxy.xy,map);
        disp(" Hit space to end ... ");
        disp("  ");
        pause;
        
        end
        
        function DemoMF(obj,whatlevel)
        % DemoMF Demo the MF process.
        
        assert(~isempty(obj.inputAxy.xy), "Need coordinates!")
        
        if whatlevel == obj.numlevels
            DemoFinalPart(obj);
        end
        
        disp("  ");
        disp(" Current level: " + whatlevel);
        disp("  ");
        
        
        map = GetPartMap(obj,whatlevel);
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
