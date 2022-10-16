classdef HIFGraph < handle
    % HIFGraph HIF algorithm based on graph partition.
    
    properties
        % Root properties.
        % The following information will be stored only in the root node.
        inputAxy_; % Input Axy.
        active_; % Whether a vertex is eliminated.
        demoHIF_ = 0; % Demonstration button.
        HIFalg_ = 1; % HIF or MF.
        
        % Graph properties.
        vtx_; % Vertices.
        sep_; % Separator vertices.
        nb_; % Neighbor vertices.
        int_; % Interior vertices.
        sk_; % Skeleton sep. We also use hat (h) to represent it.
        re_; % Redundant sep. We also use check (c) to reprsent it.
        nbsk_; % Skeleton nb.
        nbre_; % Redundant nb.
        singlesep_ = {}; % Sep which only interact with one node.
        complexsep_; % Sep which interact with more than one node.
        
        % Tree properties.
        numlevels_; % The highest level in current subtree.
        level_; % Current level, start from 0.
        seqnum_; % A node's order in its level.
        endflag_ = 0; % Whether the partition ends.
        children_ = cell(1,2); % Children nodes.
        nbnode_ = {}; % Neighbor nodes.
        nbnodeseqnum_ = []; % Neighbor nodes' seqnum.
        nbnodelevel_ = []; % Neighbor nodes' level.
        root_; % Root node.
        nbinfo_ = struct([]); % Information between a node and its nbnode after skeletonization.
        indexinfo_ = struct([]); % Index information of a node and its children.
        
        % Matrices properties.
        % For the following matrices, the fist index is row, and the second
        % index is col.
        AII_; % Interaction between int and int.
        ASI_; % Interaction between sep and int.
        ASS_; % Interaction between sep and sep.
        ANS_; % Interaction between nb and sep.
        DI_; % The D part of LDL factorization about AII.
        LI_; % The L part of LDL factorization about AII.
        AIIinvAIS_; % AIIinvAIS = AII^{-1} * ASI^{T}.
        
        % Vectors properties.
        xI_; % The int part of a vector x.
        xS_; % The sep part of a vector x.
        
    end
    
    methods
        function obj = HIFGraph(Axy, minvtx, method, level, seqnum, vtx, sep, nb)
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
            n = size(Axy.A, 1);
            vtx = 1 : 1 : n;
            sep = [];
            nb = [];
            
            obj.root_ = obj;
            obj.inputAxy_ = Axy;
            obj.active_ = ones(1, n);
        end
        
        obj.level_ = level;
        obj.seqnum_ = seqnum;
        obj.vtx_ = vtx;
        obj.sep_ = sep;
        obj.nb_ = nb;
        
        if obj.level_ == 0
            disp("Start initialization");
            obj = BuildTree(obj, Axy, minvtx, method);
            obj = SetNeighborNode(obj);
            obj = FillTree(obj, Axy.A);
            disp("End initialization");
        end
        
        end
        
        function obj = BuildTree(obj, Axy, minvtx, method)
        % BuildTree Build tree structure according to a graph partition algorithm.
        
        n = length(obj.vtx_);
        
        if n <= minvtx
            obj.numlevels_ = obj.level_;
            obj.endflag_ = 1;
            return
        end
        
        % Partition
        tmpAxy.A = Axy.A(obj.vtx_, obj.vtx_);
        if ~isempty(Axy.xy)
            tmpAxy.xy = Axy.xy(obj.vtx_,:);
        end
        [p1, p2, sep1, sep2] = GraphPart(tmpAxy, method);
        % Maybe need the following two lines.
        % sep1 = unique(sep1);
        % sep2 = unique(sep2);
        p = {p1, p2};
        partsep = {sep1, sep2};
        
        % Create children HIF.
        for iter = [1, 2]
            obj.children_{iter} = HIFGraph(Axy, minvtx, method, obj.level_ + 1, obj.seqnum_ * 2 + (iter - 1),...
                obj.vtx_(p{iter}), obj.vtx_(partsep{iter}), obj.vtx_(partsep{3 - iter}));
            obj.children_{iter}.root_ = obj.root_;
        end
        
        % Pass information to its children.
        obj = PassSeparatorNeighbor(obj, Axy.A);
        
        % Recursively buildtree.
        for iter = [1, 2]
            obj.children_{iter} = BuildTree(obj.children_{iter}, Axy, minvtx, method);
        end
        
        % Get numlevels from children when partition ends.
        obj.numlevels_ = max(obj.children_{1}.numlevels_, obj.children_{2}.numlevels_);
        
        end
        
        function obj = PassSeparatorNeighbor(obj, A)
        % PassSeparatorNeighbor Send parent's sep, nb to children.
        
        nbA = A(obj.sep_, obj.nb_);
        addsep1 = [];
        addnb1 = [];
        addsep2 = [];
        addnb2 = [];
        for i = 1 : length(obj.sep_)
            sepi = obj.sep_(i);
            for iter = [1, 2]
                if isempty(find(obj.children_{iter}.vtx_ == sepi, 1))
                    continue;
                end
                index_addnb = find(nbA(i, :) ~= 0);
                addnb = obj.nb_(index_addnb);
                if iter == 1
                    addsep1 = [addsep1, sepi];
                    addnb1 = [addnb1, addnb];
                else
                    addsep2 = [addsep2, sepi];
                    addnb2 = [addnb2, addnb];
                end
            end
        end
        obj.children_{1}.sep_ = sort(unique([obj.children_{1}.sep_, addsep1]));
        obj.children_{2}.sep_ = sort(unique([obj.children_{2}.sep_, addsep2]));
        obj.children_{1}.nb_ = sort(unique([obj.children_{1}.nb_, addnb1]));
        obj.children_{2}.nb_ = sort(unique([obj.children_{2}.nb_, addnb2]));
        
        end
        
        function obj = SetNeighborNode(obj)
        % SetNeighborNode Set nbnode.
        
        % We stand on the parent level to assign its children's nbnode.
        if obj.endflag_ == 1
            return;
        end
        
        for iter = [1, 2]
            obj_child = obj.children_{iter};
            obj_child.nbnode_{end + 1} = obj.children_{3 - iter};
            obj_child.nbnodeseqnum_(end + 1) = obj.children_{3 - iter}.seqnum_;
            obj_child.nbnodelevel_(end + 1) = obj.children_{3 - iter}.level_;
            % We only need to find nbnode from the children node of
            % parent's nbnode or parent's nbnode if it doesn't have a
            % child.
            for i = 1 : length(obj.nbnode_)
                nbnodei = obj.nbnode_{i};
                % What we need is to check whether the vtx of nbnodei's
                % chilren is in the nb of obj_child.
                if nbnodei.endflag_ == 1
                    % The nbnodei doesn't have a child. We should treat
                    % it as a nbnode.
                    if ~isempty(intersect(obj_child.nb_, nbnodei.sep_))
                        % We have to avoid add one's ancestor
                        % as its nbnode.
                        dlevel = obj_child.level_ - nbnodei.level_;
                        myseqnum = obj_child.seqnum_ / 2^dlevel;
                        if myseqnum == nbnodei.seqnum_
                            continue;
                        end
                        if isempty(intersect(find(obj_child.nbnodeseqnum_ == nbnodei.seqnum_), ...
                                find(obj_child.level_ == nbnodei.level_)))
                            obj_child.nbnode_{end+1} = nbnodei;
                            obj_child.nbnodeseqnum_(end+1) = nbnodei.seqnum_;
                            obj_child.nbnodelevel_(end+1) = nbnodei.level_;
%                             nbnodei.nbnode{end+1} = obj_child;
%                             nbnodei.nbnodeseqnum(end+1) = obj_child.seqnum;
%                             nbnodei.nbnodelevel(end+1) = obj_child.level;
                        end
                    end
                else
                    for it = [1, 2]
                        nbnodei_child = nbnodei.children_{it};
                        if ~isempty(intersect(obj_child.nb_, nbnodei_child.sep_))
                            if isempty(intersect(find(obj_child.nbnodeseqnum_ == nbnodei_child.seqnum_), ...
                                    find(obj_child.level_ == nbnodei_child.level_)))
                                obj_child.nbnode_{end + 1} = nbnodei_child;
                                obj_child.nbnodeseqnum_(end + 1) = nbnodei_child.seqnum_;
                                obj_child.nbnodelevel_(end + 1) = nbnodei_child.level_;
                            end
                        end
                    end
                end
            end
        end
        
        % Recursively setNbNode.
        for iter = [1, 2]
            obj.children_{iter} = SetNeighborNode(obj.children_{iter});
        end
        
        end
        
        function obj = FillTree(obj, A)
        % FillTree Fill tree structure with A.
        
        % Sort vtx, sep, nb.
        obj.vtx_ = sort(obj.vtx_);
        obj.sep_ = sort(obj.sep_);
        obj.nb_ = sort(obj.nb_);
        
        % We only fill the leaf nodes.
        if obj.endflag_ == 1
            % First, we set int = vtx - sep (only holds on leaf nodes).
            obj.int_ = setdiff(obj.vtx_, obj.sep_, 'sorted');
            
            % Then, we set the corresponding A**.
            obj.AII_ = full(A(obj.int_, obj.int_));
            obj.ASI_ = full(A(obj.sep_, obj.int_));
            obj.ASS_ = full(A(obj.sep_, obj.sep_));
            obj.ANS_ = full(A(obj.nb_, obj.sep_));
            
            % Set sep type.
            obj = SetSeparatorType(obj);
        else
            for iter = [1, 2]
                obj.children_{iter} = FillTree(obj.children_{iter}, A);
            end
        end
        
        end
        
        function obj = SetSeparatorType(obj)
        % SetSeparatorType Set separator type.
        
        ordersep = zeros(length(obj.sep_), 1);
        
        for k = 1 : length(obj.nbnode_)
            nodek = obj.nbnode_{k};
            [~, tmp, ~] = intersect(obj.sep_, nodek.nb_);
            obj.singlesep_{k} = tmp;
            ordersep(tmp) = ordersep(tmp) + 1;
        end
        
        obj.complexsep_ = obj.sep_;
        for k = 1 : length(obj.nbnode_)
            obj.singlesep_{k} = obj.sep_(intersect(find(ordersep == 1), obj.singlesep_{k}));
            obj.complexsep_ = setdiff(obj.complexsep_, obj.singlesep_{k});
        end
        
        end
        
        function obj = Factorization(obj, tol, HIFalg, demoHIF)
        % Factorization HIF factorization.
        
        if nargin == 1
            disp("Default tol :1e-3");
            tol = 1e-3;
            demoHIF = 0;
            HIF = 1;
        end
        
        if nargin == 2
            HIF = 1;
            demoHIF = 0;
        end
        
        if nargin == 3
            demoHIF = 0;
        end
        
        obj.HIFalg_ = HIFalg;
        obj.demoHIF_ = demoHIF;
        
        disp("Start factorization");
        
        for tmplevel = obj.numlevels_: -1 : 1
            % Sparse elimination.
            obj = RecursiveSparseElim(obj, tmplevel);
            % Demo the process.
            if obj.demoHIF_ == 1
                DemoHIF(obj, tmplevel);
            end
            % Skeletonization.
            obj = RecursiveSkel(obj, tmplevel, tol);
            % Demo the process.
            if obj.HIFalg_ ==1 && obj.demoHIF_ == 1 && tmplevel ~= obj.numlevels_
                DemoHIF(obj, tmplevel);
            end
            % Merge.
            obj = RecursiveMerge(obj, tmplevel - 1);
            % SetSeparatorType.
            obj = RecursiveSetSeparatorType(obj, tmplevel - 1);
        end
        
        % Root factorization.
        obj = RootFactorization(obj);
        % Demo the process.
        if obj.demoHIF_ == 1
            DemoHIF(obj, 0);
        end

        disp("End factorization");
        
        end
        
        function obj = RecursiveSparseElim(obj, whatlevel)
        % RecursiveSparseElim Recursively sparse elimination.
        
        if obj.level_ == whatlevel
            obj = SparseElim(obj);
        else
            if obj.endflag_ == 0
                for iter = [1, 2]
                    obj.children_{iter} = RecursiveSparseElim(obj.children_{iter}, whatlevel);
                end
            end
        end
        
        end
        
        function obj = SparseElim(obj)
        % SparseElim Sparse elimination.
        
        obj.root_.active_(obj.int_) = 0;
        % AII = LI * DI * LI^{T}.
        [obj.LI_, obj.DI_] = ldl(obj.AII_);
        % AIIinvAIS = AII^{-1} * ASI^{T}.
        obj.AIIinvAIS_ = obj.LI_ \ (obj.ASI_');
        obj.AIIinvAIS_ = obj.DI_ \ obj.AIIinvAIS_;
        obj.AIIinvAIS_ = obj.LI_' \ obj.AIIinvAIS_;
        % ASS = ASS - ASI * AII^{-1} * ASI^{T}.
        obj.ASS_ = obj.ASS_ - obj.ASI_ * obj.AIIinvAIS_;
        % ASI = 0;
        
        end
        
        function obj = RecursiveSkel(obj, whatlevel, tol)
        % RecursiveSkel Recursively skeletonization.
        
        if obj.level_ == whatlevel
            if obj.root_.HIFalg_ == 1
                obj = Skel(obj, tol);
            else
                obj = NoSkel(obj);
            end
        else
            if obj.endflag_ == 0
                for iter = [1, 2]
                    obj.children_{iter} = RecursiveSkel(obj.children_{iter} ,whatlevel, tol);
                end
            end
        end
        
        end
        
        function obj = Skel(obj, tol)
        % Skel Skeletonization.
        
        for k = 1 : length(obj.nbnode_)
            % We do skel according to nbnode.
            nodek = obj.nbnode_{k};
            if nodek.level_ ~= obj.level_ || nodek.seqnum_ < obj.seqnum_
                obj.nbinfo_(k).empty = 1;
                continue;
            end
            
            % The following data are vertices.
            sep1 = obj.singlesep_{k};
            mysep1C = setdiff(obj.sep_, sep1, 'sorted');
            nodeksep1C = setdiff(nodek.nb_, sep1, 'sorted');
            
            korder = find(nodek.nbnodeseqnum_ == obj.seqnum_);
            if length(korder) > 1
                klevel = find(nodek.nbnodelevel_ == obj.level_);
                korder = intersect(klevel, korder);
            end
            sep2 = nodek.singlesep_{korder};
            nodeksep2C = setdiff(nodek.sep_, sep2, 'sorted');
            mysep2C = setdiff(obj.nb_, sep2, 'sorted');
            
            % The following data are indices.
            [~, myindex_sep1] = ismember(sep1, obj.sep_);
            [~, nodekindex_sep1] = ismember(sep1, nodek.nb_);
            [~, myindex_sep1C] = ismember(mysep1C, obj.sep_);
            [~, myindex_mysep2C] = ismember(mysep2C, obj.nb_);
            
            [~, nodekindex_sep2] = ismember(sep2, nodek.sep_);
            [~, myindex_sep2] = ismember(sep2, obj.nb_);
            [~, nodekindex_sep2C] = ismember(nodeksep2C, nodek.sep_);
            [~, nodekindex_nodeksep1C] = ismember(nodeksep1C, nodek.nb_);
            
            % ID decomposition.
            % In the following process, the first "1" and "2" represent the
            % node or its neighbor node, the second "1" or "2" denotes sk 
            % or re.
            
            skelmtx1 = [obj.ASS_(myindex_sep1C, myindex_sep1);
                obj.ANS_(myindex_mysep2C, myindex_sep1)];
            skelmtx2 = [nodek.ASS_(nodekindex_sep2C, nodekindex_sep2);
                nodek.ANS_(nodekindex_nodeksep1C, nodekindex_sep2)];
            if isempty(skelmtx1) || isempty(skelmtx2)
                obj.nbinfo_(k).empty = 1;
                continue;
            end
            
            [T1, p11, p12] = ID(skelmtx1, tol); % skelmtx1(:, p12) = skelmtx1(:, p11) * T1.
            myindex_p11 = myindex_sep1(p11);
            myindex_p12 = myindex_sep1(p12);
            nodekindex_p11 = nodekindex_sep1(p11);
            nodekindex_p12 = nodekindex_sep1(p12);
            obj.re_ = [obj.re_, sep1(p12)];
            nodek.nbre_ = [nodek.nbre_, sep1(p12)];
            obj.nbinfo_(k).Th1c1 = T1;
            obj.nbinfo_(k).myindex_p11 = myindex_p11;
            obj.nbinfo_(k).myindex_p12 = myindex_p12;
            obj.nbinfo_(k).nodekindex_p11 = nodekindex_p11;
            obj.nbinfo_(k).nodekindex_p12 = nodekindex_p12;
            
            [T2, p21, p22] = ID(skelmtx2, tol); % skelmtx2(:, p22) = skelmtx1(:, p21) * T2.
            myindex_p21 = myindex_sep2(p21);
            myindex_p22 = myindex_sep2(p22);
            nodekindex_p21 = nodekindex_sep2(p21);
            nodekindex_p22 = nodekindex_sep2(p22);
            nodek.re_ = [nodek.re_, sep2(p22)];
            obj.nbre_ = [obj.nbre_, sep2(p22)];
            obj.nbinfo_(k).Th2c2 = T2;
            obj.nbinfo_(k).myindex_p21 = myindex_p21;
            obj.nbinfo_(k).myindex_p22 = myindex_p22;
            obj.nbinfo_(k).nodekindex_p21 = nodekindex_p21;
            obj.nbinfo_(k).nodekindex_p22 = nodekindex_p22;
            
            obj.nbinfo_(k).empty = 0;
            
            % Step 1
            Ac1h1T1 = obj.ASS_(myindex_p11, myindex_p12)' * T1;
            Ah1h1T1 = obj.ASS_(myindex_p11, myindex_p11) * T1;
            Ac2h2T2 = nodek.ASS_(nodekindex_p21, nodekindex_p22)' * T2;
            Ah2h2T2 = nodek.ASS_(nodekindex_p21, nodekindex_p21) * T2;
            % Ac1c1 = Ac1c1 - Ah1c1^{T} * Th1c1 - Th1c1^{T} * Ah1c1 + Th1c1^{T} * Ah1h1 * Th1c1.
            obj.ASS_(myindex_p12, myindex_p12) = ...
                obj.ASS_(myindex_p12, myindex_p12) - Ac1h1T1 - Ac1h1T1' + T1' * Ah1h1T1;
            % Ah1c1 = Ah1c1 - Ah1h1 * Th1c1.
            obj.ASS_(myindex_p11, myindex_p12) = ...
                obj.ASS_(myindex_p11, myindex_p12) - Ah1h1T1;
            obj.ASS_(myindex_p12, myindex_p11) = obj.ASS_(myindex_p11, myindex_p12)';
            % Ac2c1 = Ac2c1 - Ac2h1 * Th1c1 - Th2c2^{T} * Ah2c1 + Th2c2^{T} * Ah2h1 * Th1c1.
            obj.ANS_(myindex_p22, myindex_p12) = ...
                obj.ANS_(myindex_p22, myindex_p12) - obj.ANS_(myindex_p22, myindex_p11) * T1 - ...
                T2' * obj.ANS_(myindex_p21, myindex_p12) + T2' * obj.ANS_(myindex_p21,myindex_p11) * T1;
            nodek.ANS_(nodekindex_p12, nodekindex_p22) = obj.ANS_(myindex_p22, myindex_p12)';
            % Ah2c1 = Ah2c1 - Ah2h1 * Th1c1.
            obj.ANS_(myindex_p21, myindex_p12) = ...
                obj.ANS_(myindex_p21, myindex_p12) - obj.ANS_(myindex_p21, myindex_p11) * T1;
            nodek.ANS_(nodekindex_p12, nodekindex_p21) = obj.ANS_(myindex_p21, myindex_p12)';
            % Ac2h1 = Ac2h1 - Th2c2^{T} * Ah2h1.
            obj.ANS_(myindex_p22, myindex_p11) = ...
                obj.ANS_(myindex_p22, myindex_p11) - T2' * obj.ANS_(myindex_p21, myindex_p11);
            nodek.ANS_(nodekindex_p11, nodekindex_p22) = obj.ANS_(myindex_p22, myindex_p11)';
            % Ac2c2 = Ac2c2 - Ah2c2^{T} * Th2c2 - Th2c2^{T} * Ah2c2 + sTh2c2^{T} * Ah2h2 * Th2c2.
            nodek.ASS_(nodekindex_p22, nodekindex_p22) = ...
                nodek.ASS_(nodekindex_p22, nodekindex_p22) - Ac2h2T2 - Ac2h2T2' + T2' * Ah2h2T2;
            % Ah2c2 = Ah2c2 - Ah2h2 * Th2c2.
            nodek.ASS_(nodekindex_p21, nodekindex_p22) = ...
                nodek.ASS_(nodekindex_p21, nodekindex_p22) - Ah2h2T2;
            nodek.ASS_(nodekindex_p22, nodekindex_p21) = nodek.ASS_(nodekindex_p21, nodekindex_p22)';
            % Ad1c1 = Ac1d1 = 0, At1c1 = Ac1t1 = 0.
            % Ad2c1 = Ac1d2 = 0.
            % Ad1c2 = Ac2d1 = 0, At2c2 = Ac2t2 = 0;
            % Ad1c2 = Ac2d1 = 0.
            % TODO: make it zero.
            
            % Step 2
            % Ac1c1 = Lc1 * Dc1 * Lc1^{T}.
            [L1, D1] = ldl(obj.ASS_(myindex_p12, myindex_p12));
            obj.root_.active_(sep1(p12)) = 0;
            obj.nbinfo_(k).Lc1 = L1;
            obj.nbinfo_(k).Dc1 = D1;
            % Ac1c1invAc1h1 = Ac1c1^{-1} * Ah1c1^{T}.
            obj.nbinfo_(k).Ac1c1invAc1h1 = L1 \ obj.ASS_(myindex_p11, myindex_p12)';
            obj.nbinfo_(k).Ac1c1invAc1h1 = D1 \ obj.nbinfo_(k).Ac1c1invAc1h1;
            obj.nbinfo_(k).Ac1c1invAc1h1 = L1' \ obj.nbinfo_(k).Ac1c1invAc1h1;
            % Ac1c1invAc1c2 = Ac1c1^{-1} * Ac2c1^{T}.
            obj.nbinfo_(k).Ac1c1invAc1c2 = L1 \ obj.ANS_(myindex_p22,myindex_p12)';
            obj.nbinfo_(k).Ac1c1invAc1c2 = D1 \ obj.nbinfo_(k).Ac1c1invAc1c2;
            obj.nbinfo_(k).Ac1c1invAc1c2 = L1' \ obj.nbinfo_(k).Ac1c1invAc1c2;
            % Ac1c1invAc1h2 = Ac1c1^{-1} * Ah2c1^{T}.
            obj.nbinfo_(k).Ac1c1invAc1h2 = L1 \ obj.ANS_(myindex_p21,myindex_p12)';
            obj.nbinfo_(k).Ac1c1invAc1h2 = D1 \ obj.nbinfo_(k).Ac1c1invAc1h2;
            obj.nbinfo_(k).Ac1c1invAc1h2 = L1' \ obj.nbinfo_(k).Ac1c1invAc1h2;
            % Ah1h1 = Ah1h1 - Ah1c1 * Ac1c1^{-1} * Ah1c1^{T}.
            obj.ASS_(myindex_p11, myindex_p11) = ...
                obj.ASS_(myindex_p11, myindex_p11) - obj.ASS_(myindex_p11, myindex_p12) * obj.nbinfo_(k).Ac1c1invAc1h1;
            % Ac2h1 = Ac2h1 - Ac2c1 * Ac1c1^{-1} * Ah1c1^{T}.
            obj.ANS_(myindex_p22, myindex_p11) = ...
                obj.ANS_(myindex_p22, myindex_p11) - obj.ANS_(myindex_p22, myindex_p12) * obj.nbinfo_(k).Ac1c1invAc1h1;
            nodek.ANS_(nodekindex_p11, nodekindex_p22) = obj.ANS_(myindex_p22, myindex_p11)';
            % Ah2h1 = Ah2h1 - Ah2c1 * Ac1c1^{-1} * Ah1c1^{T}.
            obj.ANS_(myindex_p21, myindex_p11) = ...
                obj.ANS_(myindex_p21, myindex_p11) - obj.ANS_(myindex_p21, myindex_p12) * obj.nbinfo_(k).Ac1c1invAc1h1;
            nodek.ANS_(nodekindex_p11, nodekindex_p21) = obj.ANS_(myindex_p21, myindex_p11)';
            % Ac2c2 = Ac2c2 - Ac2c1 * Ac1c1^{-1} * Ac2c1^{T}.
            nodek.ASS_(nodekindex_p22, nodekindex_p22) = ...
                nodek.ASS_(nodekindex_p22, nodekindex_p22) - obj.ANS_(myindex_p22, myindex_p12) * obj.nbinfo_(k).Ac1c1invAc1c2;
            % Ah2c2 = Ah2c2 - Ah2c1 * Ac1c1^{-1} * Ac2c1^{T}.
            nodek.ASS_(nodekindex_p21, nodekindex_p22) = ...
                nodek.ASS_(nodekindex_p21, nodekindex_p22) - obj.ANS_(myindex_p21, myindex_p12) * obj.nbinfo_(k).Ac1c1invAc1c2;
            nodek.ASS_(nodekindex_p22, nodekindex_p21) = nodek.ASS_(nodekindex_p21, nodekindex_p22)';
            % Ah2h2 = Ah2h2 - Ah2c1 * Ac1c1^{-1} * Ah2c1^{T}.
            nodek.ASS_(nodekindex_p21, nodekindex_p21) = ...
                nodek.ASS_(nodekindex_p21, nodekindex_p21) - obj.ANS_(myindex_p21, myindex_p12) * obj.nbinfo_(k).Ac1c1invAc1h2;
            % Ah1c1 = Ac2c1 = Ah2c1 = 0.
            % TODO: make it zero.
            
            % Step 3
            % Ac2c2 = Lc2 * Dc2 * Lc2^{T};
            [L2, D2] = ldl(nodek.ASS_(nodekindex_p22, nodekindex_p22));
            obj.root_.active_(sep2(p22)) = 0;
            obj.nbinfo_(k).Lc2 = L2;
            obj.nbinfo_(k).Dc2 = D2;
            % Ac2c2invAc2h1 = Ac2c2^{-1} * Ac2h1.
            obj.nbinfo_(k).Ac2c2invAc2h1 = L2 \ obj.ANS_(myindex_p22, myindex_p11);
            obj.nbinfo_(k).Ac2c2invAc2h1 = D2 \ obj.nbinfo_(k).Ac2c2invAc2h1;
            obj.nbinfo_(k).Ac2c2invAc2h1 = L2' \ obj.nbinfo_(k).Ac2c2invAc2h1;
            % Ac2c2invAc2h2 = Ac2c2^{-1} * Ah2c2^{T}.
            obj.nbinfo_(k).Ac2c2invAc2h2 = L2  \ nodek.ASS_(nodekindex_p21, nodekindex_p22)';
            obj.nbinfo_(k).Ac2c2invAc2h2 = D2 \ obj.nbinfo_(k).Ac2c2invAc2h2;
            obj.nbinfo_(k).Ac2c2invAc2h2 = L2' \ obj.nbinfo_(k).Ac2c2invAc2h2;
            % Ah1h1 = Ah1h1 - Ac2h1^{T} * Ac2c2^{-1} * Ac2h1.
            obj.ASS_(myindex_p11, myindex_p11) = ...
                obj.ASS_(myindex_p11, myindex_p11) - obj.ANS_(myindex_p22, myindex_p11)' * obj.nbinfo_(k).Ac2c2invAc2h1;
            % Ah2h1 = Ah2h1 - Ah2c2 * Ac2c2^{-1} * Ac2h1.
            obj.ANS_(myindex_p21, myindex_p11) = ...
                obj.ANS_(myindex_p21, myindex_p11) - nodek.ASS_(nodekindex_p21, nodekindex_p22) * obj.nbinfo_(k).Ac2c2invAc2h1;
            nodek.ANS_(nodekindex_p11, nodekindex_p21) = obj.ANS_(myindex_p21,myindex_p11)';
            % Ah2h2 = Ah2h2 - Ah2c2 * Ac2c2^{-1} * Ah2c2^{T}.
            nodek.ASS_(nodekindex_p21, nodekindex_p21) = ...
                nodek.ASS_(nodekindex_p21, nodekindex_p21) - nodek.ASS_(nodekindex_p21, nodekindex_p22) * ...
                obj.nbinfo_(k).Ac2c2invAc2h2;
            % Ah2c2 = Ac2h1 = 0.
            % TODO: make it zero.
        end
        
        obj.re_ = sort(obj.re_);
        obj.sk_ = setdiff(obj.sep_, obj.re_, 'sorted');
        obj.nbre_ = sort(obj.nbre_);
        obj.nbsk_ = setdiff(obj.nb_, obj.nbre_, 'sorted');
        
        end
        
        function obj = NoSkel(obj)
        % NoSkel No skeletonization.
        
        obj.sk_ = obj.sep_;
        obj.nbsk_ = obj.nb_;
        
        end
        
        function obj = RecursiveMerge(obj, whatlevel)
        % RecursiveMerge Recusively send matrices' information from children to parent.
        
        if obj.level_ == whatlevel
            obj = Merge(obj);
        else
            if obj.endflag_ == 0
                for iter = [1, 2]
                    obj.children_{iter} = RecursiveMerge(obj.children_{iter}, whatlevel);
                end
            end
        end
        
        end
        
        function obj = Merge(obj)
        % Merge Send matrices' information from children to parent.
        
        % We stand on the parent level.
        if obj.endflag_ == 1
            return;
        end
        
        % First we tell the parent what its int, sep, nb is after we
        % eliminate the children's vtx.
        % int: children's sk - sep.
        % sep: sep \cup children's sk.
        % nb: nb \cup children's nbsk.
        tmp = [];
        for iter = [1, 2]
            obj.int_ = [obj.int_, obj.children_{iter}.sk_];
            tmp = [tmp, obj.children_{iter}.nbsk_];
        end
        obj.sep_ = intersect(obj.sep_, obj.int_, 'sorted');
        obj.int_ = setdiff(obj.int_, obj.sep_, 'sorted');
        obj.nb_ = intersect(obj.nb_, tmp, 'sorted');
        
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
        obj.AII_ = zeros(length(obj.int_));
        [~, myindex_int1, cindex_int1] = intersect(obj.int_, obj.children_{1}.sep_);
        [~, myindex_int2, cindex_int2] = intersect(obj.int_, obj.children_{2}.sep_);
        [~, myindex_int21, cindex_int21] = intersect(obj.int_, obj.children_{1}.nb_);
        obj.indexinfo_(1).myindex_int = myindex_int1;
        obj.indexinfo_(1).cindex_int = cindex_int1;
        obj.indexinfo_(2).myindex_int = myindex_int2;
        obj.indexinfo_(2).cindex_int = cindex_int2;
        obj.AII_(myindex_int1, myindex_int1) = obj.children_{1}.ASS_(cindex_int1, cindex_int1);
        obj.AII_(myindex_int2, myindex_int2) = obj.children_{2}.ASS_(cindex_int2, cindex_int2);
        obj.AII_(myindex_int21, myindex_int1) = obj.children_{1}.ANS_(cindex_int21, cindex_int1);
        obj.AII_(myindex_int1, myindex_int2) = obj.AII_(myindex_int2, myindex_int1)';
        
        % ASI
        % A sep of the parent only belongs to the sep of one of its
        % children. If an int and a sep belongs to the same child, we
        % assign ASI from the child's ASS. Otherwise, we assign ASI from
        % one child's ANS or 0.
        obj.ASI_ = zeros(length(obj.sep_), length(obj.int_));
        [~, myindex_sep1x, cindex_sep1x] = intersect(obj.sep_, obj.children_{1}.sep_);
        [~, myindex_sep1y, cindex_sep1y] = intersect(obj.sep_, obj.children_{1}.nb_);
        [~, myindex_sep2x, cindex_sep2x] = intersect(obj.sep_, obj.children_{2}.sep_);
        [~, myindex_sep2y, cindex_sep2y] = intersect(obj.sep_, obj.children_{2}.nb_);
        obj.ASI_(myindex_sep1x, myindex_int1) = obj.children_{1}.ASS_(cindex_sep1x, cindex_int1);
        obj.ASI_(myindex_sep2x, myindex_int2) = obj.children_{2}.ASS_(cindex_sep2x, cindex_int2);
        obj.ASI_(myindex_sep1y, myindex_int1) = obj.children_{1}.ANS_(cindex_sep1y, cindex_int1);
        obj.ASI_(myindex_sep2y, myindex_int2) = obj.children_{2}.ANS_(cindex_sep2y, cindex_int2);
        
        % ASS
        % If two seps belongs to the same child, we assign ASS from the
        % child's ASS. Otherwise, we assign ASS from one child's ANS or 0.
        obj.ASS_ = zeros(length(obj.sep_));
        [~, myindex_sep1, cindex_sep1] = intersect(obj.sep_, obj.children_{1}.sep_);
        [~, myindex_sep2, cindex_sep2] = intersect(obj.sep_, obj.children_{2}.sep_);
        [~, myindex_sep21, cindex_sep21] = intersect(obj.sep_, obj.children_{1}.nb_);
        obj.indexinfo_(1).myindex_sep = myindex_sep1;
        obj.indexinfo_(1).cindex_sep = cindex_sep1;
        obj.indexinfo_(2).myindex_sep = myindex_sep2;
        obj.indexinfo_(2).cindex_sep = cindex_sep2;
        obj.ASS_(myindex_sep1, myindex_sep1) = obj.children_{1}.ASS_(cindex_sep1, cindex_sep1);
        obj.ASS_(myindex_sep2, myindex_sep2) = obj.children_{2}.ASS_(cindex_sep2, cindex_sep2);
        obj.ASS_(myindex_sep21, myindex_sep1) = obj.children_{1}.ANS_(cindex_sep21, cindex_sep1);
        obj.ASS_(myindex_sep1, myindex_sep2) = obj.ASS_(myindex_sep2, myindex_sep1)';
        
        % ANS.
        % If a nb and a sep in the same child, we assign ANS from the
        % child's ANS. Otherwise, ANS= 0.
        obj.ANS_ = zeros(length(obj.nb_), length(obj.sep_));
        [~, myindex_nb1x, cindex_nb1x] = intersect(obj.nb_, obj.children_{1}.nb_);
        [~, myindex_nb2x, cindex_nb2x] = intersect(obj.nb_, obj.children_{2}.nb_);
        obj.ANS_(myindex_nb1x, myindex_sep1) = obj.children_{1}.ANS_(cindex_nb1x, cindex_sep1);
        obj.ANS_(myindex_nb2x, myindex_sep2) = obj.children_{2}.ANS_(cindex_nb2x, cindex_sep2);
        
        % Clear children information.
        for iter = [1, 2]
            obj.children_{iter} = HIFClear(obj.children_{iter});
        end
        
        end
        
        function obj = RecursiveSetSeparatorType(obj, whatlevel)
        % RecursiveSetSeparatorType Recusively set separator type.
        
        if obj.level_ == whatlevel
            obj = SetSeparatorType(obj);
        else
            if obj.endflag_ == 0
                for iter = [1, 2]
                    obj.children_{iter} = RecursiveSetSeparatorType(obj.children_{iter}, whatlevel);
                end
            end
        end
        
        end
        
        function obj = HIFClear(obj)
        % HIFClear Clear unnecessary information.
        
        obj.AII_ = [];
        obj.ASI_ = [];
        obj.ASS_ = [];
        obj.ANS_ = [];
        
        end
        
        function obj = RootFactorization(obj)
        % RootFactorization Factorization on the root.
        
        obj.root_.active_(obj.int_) = 0;
        % AII = LI * DI * LI^{T}.
        [obj.LI_, obj.DI_] = ldl(obj.AII_);
        
        end
        
        function x = HIFSolve(obj, b)
        % HIFSolve Solve Ax = b through HIF.
        
        disp("Start solve");
        
        obj = BuildVecTree(obj, b);
        
        for tmplevel = obj.numlevels_ : -1 : 1
            obj = RecursiveApplySparseElimUp(obj, tmplevel);
            obj = RecursiveApplySkelUp(obj, tmplevel);
            obj = RecursiveApplyMerge(obj, tmplevel - 1);
        end
        
        obj = RootApply(obj);
        
        for tmplevel = 1 : 1 : obj.numlevels_
            obj = RecursiveApplySplit(obj, tmplevel - 1);
            obj = RecursiveApplySkelDown(obj, tmplevel);
            obj = RecursiveApplySparseElimDown(obj, tmplevel);
        end
        
        x = zeros(size(b));
        x = GetSolution(obj, x);

        disp("End solve");
        
        end
        
        function obj = BuildVecTree(obj, b)
        % BuildVecTree Fill b in our HIF tree.
        
        if obj.endflag_ == 1
            obj.xI_ = b(obj.int_,:);
            obj.xS_ = b(obj.sep_,:);
        else
            for iter = [1, 2]
                obj.children_{iter} = BuildVecTree(obj.children_{iter}, b);
            end
        end
        
        end
        
        function obj = RecursiveApplySparseElimUp(obj, whatlevel)
        % RecursiveApplySparseElimUp Phase 1 for applying sparse elimination recusively.
        
        if obj.level_ == whatlevel
            obj = ApplySparseElimUp(obj);
        else
            if obj.endflag_ == 0
                for iter = [1, 2]
                    obj.children_{iter} = RecursiveApplySparseElimUp(obj.children_{iter}, whatlevel);
                end
            end
        end
        
        end
        
        function obj = ApplySparseElimUp(obj)
        % ApplySparseElimUp Phase 1 for applying sparse elimination.
        
        % xS = xS - (AII^{-1} * ASI^{T})^{T} * xI.
        obj.xS_ = obj.xS_ - obj.AIIinvAIS_' * obj.xI_;
        % xI = LI^{-1} * xI.
        obj.xI_ = obj.LI_ \ obj.xI_;
        
        % xI = DI^{-1} * xI. We only apply D once.
        obj.xI_ = obj.DI_ \ obj.xI_;
        
        end
        
        function obj = RecursiveApplySkelUp(obj, whatlevel)
        % RecursiveApplySkelUp Phase 1 for applying skeletonization recusively.
        
        if obj.level_ == whatlevel
            obj = ApplySkelUp(obj);
        else
            if obj.endflag_ == 0
                for iter = [1, 2]
                    obj.children_{iter} = RecursiveApplySkelUp(obj.children_{iter}, whatlevel);
                end
            end
        end
        
        end
        
        function obj = ApplySkelUp(obj)
        % ApplySkelUp Phase 1 for applying skeletonization.
        
        for k = 1 : length(obj.nbnode_)
            nbnodek = obj.nbnode_{k};
            if isempty(obj.nbinfo_)
                continue;
            end
            nbinfok = obj.nbinfo_(k);
            if nbinfok.empty
                continue;
            end
                
            % Step 1
            % xc1 = xc1 - Th1c1^{T} * xh1.
            obj.xS_(nbinfok.myindex_p12, :) = ...
                obj.xS_(nbinfok.myindex_p12, :) - nbinfok.Th1c1' * obj.xS_(nbinfok.myindex_p11, :);
            % xc2 = xc2 - Th2c2^{T} * xh2.
            nbnodek.xS_(nbinfok.nodekindex_p22, :) = ...
                nbnodek.xS_(nbinfok.nodekindex_p22, :) - nbinfok.Th2c2' * nbnodek.xS_(nbinfok.nodekindex_p21, :);
            
            % Step 2
            % xh1 = xh1 - (Ac1c1^{-1} * Ah1c1^{T})^{T} * xc1.
            obj.xS_(nbinfok.myindex_p11, :) = ...
                obj.xS_(nbinfok.myindex_p11, :) - nbinfok.Ac1c1invAc1h1' * obj.xS_(nbinfok.myindex_p12, :);
            % xc2 = xc2 - (Ac1c1^{-1} * Ac2c1^{T})^{T} * xc1.
            nbnodek.xS_(nbinfok.nodekindex_p22, :) = ...
                nbnodek.xS_(nbinfok.nodekindex_p22, :) - nbinfok.Ac1c1invAc1c2' * obj.xS_(nbinfok.myindex_p12, :);
            % xh2 = xh2 - (Ac1c1^{-1} * Ah2c1^{T})^{T} * xc1.
            nbnodek.xS_(nbinfok.nodekindex_p21, :) = ...
                nbnodek.xS_(nbinfok.nodekindex_p21, :) - nbinfok.Ac1c1invAc1h2' * obj.xS_(nbinfok.myindex_p12, :);
            % xc1 = Lc1^{-1} * xc1.
            obj.xS_(nbinfok.myindex_p12, :) = ...
                nbinfok.Lc1 \ obj.xS_(nbinfok.myindex_p12, :);
            
            % Step 3
            % xh1 = xh1 - (Ac2c2^{-1} * Ac2h1)^{T} * xc2.
            obj.xS_(nbinfok.myindex_p11, :) = ...
                obj.xS_(nbinfok.myindex_p11, :) - nbinfok.Ac2c2invAc2h1' * nbnodek.xS_(nbinfok.nodekindex_p22, :);
            % xh2 = xh2 - (Ac2c2^{-1} * Ah2c2^{T})^{T} * xc2.
            nbnodek.xS_(nbinfok.nodekindex_p21, :) = ...
                nbnodek.xS_(nbinfok.nodekindex_p21, :) - nbinfok.Ac2c2invAc2h2' * nbnodek.xS_(nbinfok.nodekindex_p22, :);
            % xc2 = Lc2^{-1} * xc2.
            nbnodek.xS_(nbinfok.nodekindex_p22, :) = ...
                nbinfok.Lc2 \ nbnodek.xS_(nbinfok.nodekindex_p22, :);
            
            % xc1 = Dc1^{-1} * xc1. xc2 = Dc2^{-1} * xc2. We only apply D once.
            obj.xS_(nbinfok.myindex_p12, :) = ...
                nbinfok.Dc1 \ obj.xS_(nbinfok.myindex_p12, :);
            nbnodek.xS_(nbinfok.nodekindex_p22, :) = ...
                nbinfok.Dc2 \ nbnodek.xS_(nbinfok.nodekindex_p22, :);
        end
        
        end
        
        function obj = RecursiveApplyMerge(obj, whatlevel)
        % RecursiveApplyMerge Recusively send vectors' information from children to parent.
        
        if obj.level_ == whatlevel
            obj = ApplyMerge(obj);
        else
            if obj.endflag_ == 0
                for iter = [1, 2]
                    obj.children_{iter} = RecursiveApplyMerge(obj.children_{iter}, whatlevel);
                end
            end
        end
        
        end
        
        function obj = ApplyMerge(obj)
        % ApplyMerge Send vectors' information from children to parent.
        
        % We stand on the parent level.
        if obj.endflag_ == 1
            return;
        end
        
        % We only need to assign the corresponding vectors.
        width = size(obj.children_{1}.xS_, 2);
        
        % xI.
        % An int of the parent only belongs to the sep of one of its
        % children. We assign xI from the child's xS.
        obj.xI_ = zeros(length(obj.int_), width);
        for iter = [1, 2]
            obj.xI_(obj.indexinfo_(iter).myindex_int, :) = obj.children_{iter}.xS_(obj.indexinfo_(iter).cindex_int, :);
        end
        
        % xS.
        % A sep of the parent only belongs to the sep of one of its
        % children. We assign xS from the child's xS.
        obj.xS_ = zeros(length(obj.sep_), width);
        for iter = [1, 2]
            obj.xS_(obj.indexinfo_(iter).myindex_sep, :) = obj.children_{iter}.xS_(obj.indexinfo_(iter).cindex_sep, :);
        end
        
        end
        
        function obj = RootApply(obj)
        % RootApply Apply on the root.
        
        % xI = LI^{-1} * xI.
        obj.xI_ = obj.LI_ \ obj.xI_;
        
        % xI = DI^{-1} * xI.
        obj.xI_ = obj.DI_ \ obj.xI_;
        
        % xI = LI^{-T} * xI.
        obj.xI_ = obj.LI_' \ obj.xI_;
        
        end
        
        function obj = RecursiveApplySplit(obj, whatlevel)
        % RecursiveApplySplit Recusively send vectors' information from parent to children.
        
        if obj.level_ == whatlevel
            obj = ApplySplit(obj);
        else
            if obj.endflag_ == 0
                for iter = [1, 2]
                    obj.children_{iter} = RecursiveApplySplit(obj.children_{iter}, whatlevel);
                end
            end
        end
        
        end
        
        function obj = ApplySplit(obj)
        % ApplySplit Send vectors' information from parent to children.
        
        % We stand on the parent level.
        if obj.endflag_ == 1
            return;
        end
        
        % We only need to assign the corresponding vectors of the children.
        
        for iter = [1, 2]
            obj.children_{iter}.xS_(obj.indexinfo_(iter).cindex_int, :) = obj.xI_(obj.indexinfo_(iter).myindex_int, :);
            obj.children_{iter}.xS_(obj.indexinfo_(iter).cindex_sep, :) = obj.xS_(obj.indexinfo_(iter).myindex_sep, :);
        end
        
        end
        
        function obj = RecursiveApplySkelDown(obj, whatlevel)
        % RecursiveApplySkelDown Phase 2 for applying skeletonization recursively.
        
        if obj.level_ == whatlevel
            obj = ApplySkelDown(obj);
        else
            if obj.endflag_ == 0
                for iter = [1, 2]
                    obj.children_{iter} = RecursiveApplySkelDown(obj.children_{iter}, whatlevel);
                end
            end
        end
        
        end
        
        function obj = ApplySkelDown(obj)
        % ApplySkelDown Phase 2 for applying skeletonization.
        
        for k = 1 : length(obj.nbnode_)
            nbnodek = obj.nbnode_{k};
            if isempty(obj.nbinfo_)
                continue;
            end
            nbinfok = obj.nbinfo_(k);
            if nbinfok.empty
                continue;
            end
            
            % Step 3
            % xc2 = Lc2^{-T} * xc2 - (Ac2c2^{-1} * Ac2h1) * xh1 - (Ac2c2^{-1} * Ah2c2^{T}) * xh2.
            nbnodek.xS_(nbinfok.nodekindex_p22, :) = ...
                nbinfok.Lc2' \ nbnodek.xS_(nbinfok.nodekindex_p22, :) - ...
                nbinfok.Ac2c2invAc2h1 * obj.xS_(nbinfok.myindex_p11, :) - ...
                nbinfok.Ac2c2invAc2h2 * nbnodek.xS_(nbinfok.nodekindex_p21, :);
            
            % Step 2
            % xc1 = Lc1^{-T} * xc1 - (Ac1c1^{-1} * Ah1c1^{T}) * xh1 - (Ac1c1^{-1} * Ac2c1^{T}) * xc2 - (Ac1c1^{-1} * Ah2c1^{T}) * xh2.
            obj.xS_(nbinfok.myindex_p12, :) = ...
                nbinfok.Lc1' \ obj.xS_(nbinfok.myindex_p12, :) - ...
                nbinfok.Ac1c1invAc1h1 * obj.xS_(nbinfok.myindex_p11, :) - ...
                nbinfok.Ac1c1invAc1c2 * nbnodek.xS_(nbinfok.nodekindex_p22, :) - ...
                nbinfok.Ac1c1invAc1h2 * nbnodek.xS_(nbinfok.nodekindex_p21, :);
            
            % Step 1
            % xh1 = xh1 - Th1c1 * xc1.
            obj.xS_(nbinfok.myindex_p11, :) = ...
                obj.xS_(nbinfok.myindex_p11, :) - nbinfok.Th1c1 * obj.xS_(nbinfok.myindex_p12, :);
            % xh2 = xh2 - Th2c2 * xc2.
            nbnodek.xS_(nbinfok.nodekindex_p21, :) = ...
                nbnodek.xS_(nbinfok.nodekindex_p21, :) - nbinfok.Th2c2 * nbnodek.xS_(nbinfok.nodekindex_p22, :);
        end
        
        
        end
        
        function obj = RecursiveApplySparseElimDown(obj, whatlevel)
        % RecursiveApplySparseElimDown Phase 2 for applying sparse elimination recursively.
        
        if obj.level_ == whatlevel
            obj = ApplySparseElimDown(obj);
        else
            if obj.endflag_ == 0
                for iter = [1, 2]
                    obj.children_{iter} = RecursiveApplySparseElimDown(obj.children_{iter}, whatlevel);
                end
            end
        end
        
        end
        
        function obj = ApplySparseElimDown(obj)
        % ApplySparseElimDown Phase 2 for applying sparse elimination.
        
        % xI = LI^{-T} * xI - (AII^{-1} * ASI^{T}) * xS.
        obj.xI_ = ...
            obj.LI_' \ obj.xI_ - obj.AIIinvAIS_ * obj.xS_;
        
        end
        
        function x = GetSolution(obj, x)
        % GetSolution Get the final soluion through the tree structure.
        
        if obj.endflag_ == 1
            x(obj.int_, :) = obj.xI_;
            x(obj.sep_, :) = obj.xS_;
            obj.xI_ = [];
            obj.xS_ = [];
        else
            obj.xI_ = [];
            obj.xS_ = [];
            for iter = [1, 2]
                x = GetSolution(obj.children_{iter}, x);
            end
        end
        
        end
        
        function DemoPart(obj)
        % DemoPart Demo the partition process.
        
        assert(~isempty(obj.inputAxy_.xy), "Need coordinates!")
        
        disp("Start demo the partition");
        
        fig = figure();
        clf reset;
        colordef(fig, 'black');
        gplotg(obj.inputAxy_.A, obj.inputAxy_.xy);
        
        disp("Hit space to continue...");
        for tmplevel = 0 : obj.numlevels_
            disp(" Current level: " + tmplevel);
            
            map = GetPartMap(obj, tmplevel);
            gplotmap(obj.inputAxy_.A, obj.inputAxy_.xy, map);
            if tmplevel ~= obj.numlevels_
                disp("Hit space to continue...");
            else
                disp("Hit space to end...");
            end
            pause;
        end
        
        end
        
        function DemoLevelPart(obj, whatlevel)
        % DemoLevelPart Demo the specified level partition.
        
        assert(~isempty(obj.inputAxy_.xy), "Need coordinates!")
        
        disp(" Current level: " + whatlevel);
        
        fig = figure();
        clf reset;
        colordef(fig, 'black');
        gplotg(obj.inputAxy_.A, obj.inputAxy_.xy);
        
        disp("Hit space to continue...");
        pause;
        map = GetPartMap(obj, whatlevel);
        gplotmap(obj.inputAxy_.A, obj.inputAxy_.xy, map);
        disp("Hit space to end ...");
        pause;
        
        end
        
        function map = GetPartMap(obj, whatlevel, map)
        % GetPartMap Get the map of the partition.
        
        if nargin == 2
            map = [];
        end
        
        if obj.level_ == 0
            n = size(obj.inputAxy_.xy, 1);
            map = zeros(1, n);
        end
        
        if obj.level_ == whatlevel || obj.endflag_ == 1
            map(obj.vtx_) = max(map)+1;
            return;
        else
            for iter = [1, 2]
                map = GetPartMap(obj.children_{iter}, whatlevel, map);
            end
        end
        
        end
        
        function DemoFinalPart(obj)
        % DemoFinalPart Demo of the final partition.
        
        assert(~isempty(obj.inputAxy_.xy), "Need coordinates!")
        
        disp(" Start demo the final partition ");
        
        fig = figure();
        clf reset;
        colordef(fig, 'black');
        gplotg(speye(size(obj.inputAxy_.A)), obj.inputAxy_.xy);
        
        map = GetPartMap(obj, obj.numlevels_);
        parts = unique(map);
        nparts = length(parts);
        tmpA = obj.inputAxy_.A;
        for i = 1 : nparts
            for j = i + 1 : nparts
                idi = find(map == i);
                idj = find(map == j);
                tmpA(idi, idj) = 0;
                tmpA(idj, idi) = 0;
            end
        end
        gplotmap(tmpA, obj.inputAxy_.xy, map);
        disp("Hit space to end ...");
        pause;
        
        end
        
        function DemoHIF(obj, whatlevel)
        % DemoHIF Demo the HIF process.
        
        assert(~isempty(obj.inputAxy_.xy), "Need coordinates!")
        
        if whatlevel == obj.numlevels_
            DemoFinalPart(obj);
        end
        
        disp("Current level: " + whatlevel);
        
        fig = figure();
        clf reset;
        colordef(fig, 'black');
        
        map = GetPartMap(obj, whatlevel);
        tmpactive = find(obj.active_ > 0);
        parts = unique(map);
        nparts = length(parts);
        tmpA = obj.inputAxy_.A;
        for i = 1 : nparts
            for j = i + 1 : nparts
                idi = find(map == i);
                idj = find(map == j);
                tmpA(idi, idj) = 0;
                tmpA(idj, idi) = 0;
            end
        end
        tmpA = tmpA(tmpactive, tmpactive);
        tmpxy = obj.inputAxy_.xy(tmpactive, :);
        tmpmap = map(1, tmpactive);
        gplotmap(tmpA, tmpxy, tmpmap);
        if size(obj.inputAxy_.xy, 2) == 3
            view(3);
        end
        if whatlevel ~= 0
            disp("Hit space to continue...");
        else
            disp("Hit space to end ...");
        end
        pause;
        
        end
        
        function levelVec = ReadLevel(obj, levelVec)
        % ReadLevel Obtain level vectors.
        
        if nargin == 1
            levelVec = [];
        end
        
        if obj.level_ == 0
            levelVec = zeros(obj.numlevels_, 1);
        end
        
        if obj.endflag_ == 1
            levelVec(obj.level_) = levelVec(obj.level_) + 1;
            return;
        else
            for iter = [1, 2]
                levelVec = ReadLevel(obj.children_{iter}, levelVec);
            end
        end
        
        end
        
    end
    
end
