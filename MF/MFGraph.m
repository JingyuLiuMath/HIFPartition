classdef MFGraph < handle
    % MFGraph Multifrontal algorithm on general graphs.
    
    properties
        
        % Root properties.
        
        % The following information will be stored only in the root node.
        
        inputAxy; % Input Axy.
        active; % Whether a vtx is eliminated.
        inputVec; % Input vec.
        solution; % Solution.
        demoMF = 0; % Whether to demo our MF process.
        
        % Graph properties.
        
        Axy; % Adjacency matrix (and coordinates of vertices).
        vtx; % Vertices on the graph.
        sep; % Separator vertices.
        nb; % Neighbor vertices.
        int; % Interior vertices.
        nbA; % Adjacency matrix of sep (row) and nb (col).
        
        % Partition properties.
        
        numLevels; % Total number of levels.
        level; % Current level, start from 0.
        seqNum; % A node's order in its level.
        endFlag = 0; % Whether the partition ends.
        
        % Tree properties.
        
        parent; % Parent node.
        children = cell(1,2); % Children nodes.
        nbNode = {}; % Neighbor nodes. In fact, we don't need this in MF.
        root; % Root node.
        
        % Matrices properties.
        
        % For the following matrices, the fist index is row, and the second
        % index is col.
        
        AII; % Interaction between int and int.
        ASI; % Interaction between sep and int.
        ASS; % Interaction between sep and sep.
        ANS; % Interaction between nb and sep.
        AIIinv; % AIIinv = L^{-T} where AII = L * L^{T}.
        AIIinvAIS; % AIIinvAIS = AII^{-1} * ASI^{T}
        
        % Vectors properties.
        
        xI; % The int part of a vector x.
        xS; % The sep part of a vector x.
        
    end
    
    methods
        
        function obj = MFGraph(Axy,level,seqNum,vtx,sep,nb,nbA)
        % MFGraph Create a MF class.
        
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
        
        function obj = BuildTree(obj)
        % BuildTree Build tree structure according to a graph partition algorithm.
        
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
        
        if isempty(obj.Axy.xy)
            % Specpart.
            [p1,p2,sep1,sep2] = specpart(obj.Axy.A);
            sep1 = unique(sep1);
            sep2 = unique(sep2);
            child1Axy.A = obj.Axy.A(p1,p1); child1Axy.xy = [];
            child2Axy.A = obj.Axy.A(p2,p2); child2Axy.xy = [];
        else
            % Geopart.
            [p1,p2,sep1,sep2] = geopart(obj.Axy.A,obj.Axy.xy,numtries);
            sep1 = unique(sep1);
            sep2 = unique(sep2);
            child1Axy.A = obj.Axy.A(p1,p1); child1Axy.xy = obj.Axy.xy(p1,:);
            child2Axy.A = obj.Axy.A(p2,p2); child2Axy.xy = obj.Axy.xy(p2,:);
        end
        
        % Create children MF.
        obj.children{1} = MFGraph(child1Axy,obj.level+1,obj.seqNum*2,...
            obj.vtx(p1),obj.vtx(sep1),obj.vtx(sep2),obj.Axy.A(sep1,sep2));
        obj.children{1}.root = obj.root;
        obj.children{2} = MFGraph(child2Axy,obj.level+1,obj.seqNum*2+1,...
            obj.vtx(p2),obj.vtx(sep2),obj.vtx(sep1),obj.Axy.A(sep2,sep1));
        obj.children{2}.root = obj.root;
        
        % Pass information to its children.
        obj = Pass(obj);
        
        % Recursively buildtree.
        for iter = [1,2]
            obj.children{iter} = BuildTree(obj.children{iter});
            obj.children{iter}.parent = obj;
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
            % We only need to find nbNode from the children node of
            % parent's nbNode.
            if ~isempty(obj.nbNode)
                for i = 1:length(obj.nbNode)
                    nbNodei = obj.nbNode{i};
                    % What we need is to check whether the vtx of nbNodei's
                    % chilren is in the nb of obj_child.
                    for it = [1,2]
                        nbNodei_child = nbNodei.children{it};
                        if ~isempty(intersect(obj_child.nb, nbNodei_child.vtx))
                            obj_child.nbNode{end+1} = nbNodei_child;
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
        % In fact, they can be clear after being passed to children.
        obj.Axy = []; obj.nbA = [];
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
            rootMF = obj.root;
            obj.AII = rootMF.inputAxy.A(obj.int,obj.int);
            obj.ASI = rootMF.inputAxy.A(obj.sep,obj.int);
            obj.ASS = rootMF.inputAxy.A(obj.sep,obj.sep);
            obj.ANS = rootMF.inputAxy.A(obj.nb,obj.sep);
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
        
        for tmplevel = obj.numLevels:-1:1
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
            for iter = [1,2]
                obj.children{iter} = RecursiveSparseElim(obj.children{iter},whatlevel);
            end
        end
        
        end
        
        function obj = SparseElim(obj)
        % SparseElim Sparse elimination.
        
        obj.root.active(obj.int) = 0;
        L = chol(obj.AII,'lower'); % AII = L * L^T.
        obj.AIIinv = L'\eye(size(L,1)); % AIIinv = L^{-T}.
        obj.AIIinvAIS = L\(obj.ASI');
        obj.AIIinvAIS = L'\obj.AIIinvAIS; % AIIinvAIS = AII^{-1} * ASI^{T}.
        obj.ASS = obj.ASS - obj.ASI*obj.AIIinvAIS; % ASS = ASS - ASI * AII^{-1} * ASI^{T}.
        
        end
        
        function obj = RecursiveMerge(obj,whatlevel)
        % RecursiveMerge Recusively send information from children to parent.
        
        if obj.level == whatlevel
            obj = Merge(obj);
        else
            for iter = [1,2]
                obj.children{iter} = RecursiveMerge(obj.children{iter},whatlevel);
            end
        end
        
        end
        
        function obj = Merge(obj)
        % Merge Send matrices' information from children to parent.
        
        % We stand on the parent level.
        
        % First we tell the parent what its int is after we eliminate
        % the children's vtx.
        % int: children's sep - sep
        
        for iter = [1,2]
            obj.int = [obj.int,obj.children{iter}.sep];
        end
        obj.int = setdiff(obj.int,obj.sep,'sorted');
        
        % Next we assign the corresponding matrices.
        % From child to parent.
        % int: children's sep and children's nb
        % sep: children's sep and children's nb
        % nb: children's nb
        
        obj.AII = zeros(length(obj.int));
        % An int of the parent only belongs to the sep of one of its
        % children. If two ints belong to the same child, we assign AII
        % from the child's ASS. Otherwise, we assign AII from one child's
        % ANS or 0.
        for j = 1:length(obj.int)
            intj = obj.int(j);
            if find(obj.children{1}.vtx == intj)
                where_intj = 1;
            else
                where_intj = 2;
            end
            index_intj = find(obj.children{where_intj}.sep == intj);
            for i = 1:length(obj.int)
                inti = obj.int(i);
                index_inti = find(obj.children{where_intj}.sep == inti,1);
                if isempty(index_inti)
                    index_inti = find(obj.children{where_intj}.nb == inti,1);
                    if isempty(index_inti)
                        obj.AII(i,j) = 0;
                    else
                        obj.AII(i,j) = obj.children{where_intj}.ANS(index_inti,index_intj);
                    end
                else
                    obj.AII(i,j) = obj.children{where_intj}.ASS(index_inti,index_intj);
                end
            end
        end
        
        obj.ASI = zeros(length(obj.sep),length(obj.int));
        % A sep of the parent only belongs to the sep of one of its
        % children. If an int and a sep belongs to the same child, we
        % assign ASI from the child's ASS. Otherwise, we assign ASI from
        % the int child's ANS or 0.
        for j = 1:length(obj.int)
            intj = obj.int(j);
            if find(obj.children{1}.vtx == intj)
                where_intj = 1;
            else
                where_intj = 2;
            end
            index_intj = find(obj.children{where_intj}.sep == intj);
            for i = 1:length(obj.sep)
                sepi = obj.sep(i);
                index_sepi = find(obj.children{where_intj}.sep == sepi,1);
                if isempty(index_sepi)
                    index_sepi = find(obj.children{where_intj}.nb == sepi,1);
                    if isempty(index_sepi)
                        obj.ASI(i,j) = 0;
                    else
                        obj.ASI(i,j) = obj.children{where_intj}.ANS(index_sepi,index_intj);
                    end
                else
                    obj.ASI(i,j) = obj.children{where_intj}.ASS(index_sepi,index_intj);
                end
            end
        end
        
        obj.ASS = zeros(length(obj.sep));
        % If two seps belongs to the same child, we assign ASS from the
        % child's ASS. Otherwise, we assign ASS from one child's ANN or 0.
        for j = 1:length(obj.sep)
            sepj = obj.sep(j);
            if find(obj.children{1}.vtx == sepj)
                where_sepj = 1;
            else
                where_sepj = 2;
            end
            index_sepj = find(obj.children{where_sepj}.sep == sepj);
            for i = 1:length(obj.sep)
                sepi = obj.sep(i);
                index_sepi = find(obj.children{where_sepj}.sep == sepi,1);
                if isempty(index_sepi)
                    index_sepi = find(obj.children{where_sepj}.nb == sepi,1);
                    if isempty(index_sepi)
                        obj.ASS(i,j) = 0;
                    else
                        obj.ASS(i,j) = obj.children{where_sepj}.ANS(index_sepi,index_sepj);
                    end
                else
                    obj.ASS(i,j) = obj.children{where_sepj}.ASS(index_sepi,index_sepj);
                end
            end
        end
        
        obj.ANS = zeros(length(obj.nb),length(obj.sep));
        % If a nb and a sep in the same child, we assign ANS from the
        % child's ANS. Otherwise, ANS= 0.
        for j = 1:length(obj.sep)
            sepj = obj.sep(j);
            if find(obj.children{1}.vtx == sepj)
                where_sepj = 1;
            else
                where_sepj = 2;
            end
            index_sepj = find(obj.children{where_sepj}.sep == sepj);
            for i = 1:length(obj.nb)
                nbi = obj.nb(i);
                index_nbi = find(obj.children{where_sepj}.nb == nbi,1);
                if isempty(index_nbi)
                    obj.ANS(i,j) = 0;
                else
                    obj.ANS(i,j) = obj.children{where_sepj}.ANS(index_nbi,index_sepj);
                end
            end
        end
        
        end
        
        function obj = RootFactorization(obj)
        % RootFactorization Factorization on the root.
        
        obj.root.active(obj.int) = 0;
        L = chol(obj.AII,'lower'); % AII = L * L^T.
        obj.AIIinv = L'\eye(size(L,1)); % AIIinv = L^{-T}
        
        end
        
        function obj = MFSolve(obj,b)
        % MFSolve Solve Ax = b through MF.
        
        disp("  ");
        disp(" Start solve ");
        disp("  ");
        
        obj.inputVec = b;
        
        obj = BuildVecTree(obj,b);
        
        for tmplevel = obj.numLevels:-1:1
            obj = RecursiveApplyUp(obj,tmplevel);
            obj = RecursiveApplyMerge(obj,tmplevel-1);
        end
        
        obj = RootApply(obj);
        
        for tmplevel = 1:1:obj.numLevels
            obj = RecursiveApplySplit(obj,tmplevel-1);
            obj = RecursiveApplyDown(obj,tmplevel);
        end
        
        obj = GetSolution(obj);
        
        disp("  ");
        disp(" End solve ");
        disp(" The solution is stored in obj.solution.");
        disp("  ");
        
        end
        
        function obj = BuildVecTree(obj,b)
        % BuildVecTree Fill b in our MF tree.
        
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
        
        function obj = RecursiveApplyUp(obj,whatlevel)
        % RecursiveApplyUp Phase 1 for applying MF recusively.
        
        if obj.level == whatlevel
            obj = ApplyUp(obj);
        else
            for iter = [1,2]
                obj.children{iter} = RecursiveApplyUp(obj.children{iter}, whatlevel);
            end
        end
        
        end
        
        function obj = ApplyUp(obj)
        %ApplyUp Phase 1 for applying MF
        
        obj.xS = obj.xS - obj.AIIinvAIS'*obj.xI;
        obj.xI = obj.AIIinv'*obj.xI;
        
        end
        
        function obj = RecursiveApplyMerge(obj,whatlevel)
        % RecursiveApplyMerge Recusively send vectors' information from children to parent.
        
        if obj.level == whatlevel
            obj = ApplyMerge(obj);
        else
            for iter = [1,2]
                obj.children{iter} = RecursiveApplyMerge(obj.children{iter},whatlevel);
            end
        end
        
        end
        
        function obj = ApplyMerge(obj)
        % ApplyMerge Send vectors' information from children to parent.
        
        % We stand on the parent level.
        
        % We have specified the parent's int. So we only to assign the
        % corresponding vectors.
        
        obj.xI = zeros(length(obj.int),1);
        % An int of the parent only belongs to the sep of one of its
        % children. We get xI from the children's xS.
        for j =1:length(obj.int)
            intj = obj.int(j);
            if find(obj.children{1}.vtx == intj)
                where_intj = 1;
            else
                where_intj = 2;
            end
            index_intj = find(obj.children{where_intj}.sep == intj);
            obj.xI(j) = obj.children{where_intj}.xS(index_intj);
        end
        
        obj.xS = zeros(length(obj.sep),1);
        % A sep of the parent only belongs to the sep of one of its
        % children. We get xS from the children's xS.
        for j = 1:length(obj.sep)
            sepj = obj.sep(j);
            if find(obj.children{1}.vtx == sepj)
                where_sepj = 1;
            else
                where_sepj = 2;
            end
            index_sepj = find(obj.children{where_sepj}.sep == sepj);
            obj.xS(j) = obj.children{where_sepj}.xS(index_sepj);
        end
        
        end
        
        function obj = RootApply(obj)
        % RootApply Apply on the root.
        
        obj.xI = obj.AIIinv * (obj.AIIinv' * obj.xI);
        
        end
        
        function obj = RecursiveApplySplit(obj,whatlevel)
        % RecursiveApplySplit Recusively send vectors' information from parent to children.
        
        if obj.level == whatlevel
            obj = ApplySplit(obj);
        else
            for iter = [1,2]
                obj.children{iter} = RecursiveApplySplit(obj.children{iter},whatlevel);
            end
        end
        
        end
        
        function obj = ApplySplit(obj)
        % ApplySplit Send vectors' information from parent to children.
        
        % We stand on the parent level.
        
        % We only need to assign the corresponding vectors of the children.
        
        % xI
        for j = 1:length(obj.int)
            intj = obj.int(j);
            if find(obj.children{1}.vtx == intj)
                where_intj = 1;
            else
                where_intj = 2;
            end
            index_intj = find(obj.children{where_intj}.sep == intj);
            obj.children{where_intj}.xS(index_intj) = obj.xI(j);
        end
        
        % xS
        for j = 1:length(obj.sep)
            sepj = obj.sep(j);
            if find(obj.children{1}.vtx == sepj)
                where_sepj = 1;
            else
                where_sepj = 2;
            end
            index_sepj = find(obj.children{where_sepj}.sep == sepj);
            obj.children{where_sepj}.xS(index_sepj) = obj.xS(j);
        end
        
        end
        
        function obj = RecursiveApplyDown(obj,whatlevel)
        % RecursiveApplyDown Phase 2 for applying MF recursively.
        
        if obj.level == whatlevel
            obj = ApplyDown(obj);
        else
            for iter = [1,2]
                obj.children{iter} = RecursiveApplyDown(obj.children{iter}, whatlevel);
            end
        end
        
        end
        
        function obj = ApplyDown(obj)
        % ApplyDown Phase 2 for applying MF.
        
        obj.xI = obj.AIIinv*obj.xI - obj.AIIinvAIS*obj.xS;
        
        end
        
        function obj = GetSolution(obj)
        % GetSolution Get the final soluion through the tree structure.
        
        if obj.level == 0
            obj.solution = zeros(length(obj.int),1);
        end
        
        if obj.endFlag == 0
            obj.root.solution(obj.int) = obj.xI;
            for iter = [1,2]
                obj.children{iter} = GetSolution(obj.children{iter});
            end
        else
            obj.root.solution(obj.int) = obj.xI;
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

        function DemoMF(obj,whatlevel)
        % DemoMF Demo the MF process.
        
        assert(~isempty(obj.inputAxy.xy), "Need coordinates!")
        
        if whatlevel == obj.numLevels
            DemoFinalPart(obj);
        end
        
        disp("  ");
        disp(" Current level: " + whatlevel);
        disp("  ");
        
        
        map = GetMFMap(obj,whatlevel);
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
        
        function map = GetMFMap(obj,whatlevel,map)
        % GetMFMap Recursively get MF map.
        
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
                map = GetMFMap(obj.children{iter},whatlevel,map);
            end
        end
        
        end
        
    end
    
end
