classdef mesh < handle
    %GEOMESH Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Mesh node map
        nodMap;
        % Meh element map
        tri3Map;
        tri6Map;
        % Mesh bc map
        bcVMap;
        % DoF counting
        nextDoF = 1;
        % Flags
        isPresMain = false;
        isVeloMain = false;
        % VTK info
        vtkNodCor;
        vtkTri3Con;
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        function obj = mesh()
            % Node maps
            obj.nodMap  = containers.Map('KeyType', 'double', 'ValueType', 'any');
            % Element maps
            obj.tri3Map = containers.Map('KeyType', 'double', 'ValueType', 'any');
            obj.tri6Map = containers.Map('KeyType', 'double', 'ValueType', 'any');
            % BCs maps
            obj.bcVMap   = containers.Map('KeyType', 'double', 'ValueType', 'any');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Share tri3Map to velocity mesh
        function shareMeshP(obj,tri3Map)
            if obj.isVeloMain == true
                obj.tri3Map = tri3Map;
            else
                disp('Share not allowed! Mesh main physics is pressure!');
                return
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Share tri6Map to pressure mesh
        function shareMeshV(obj,tri6Map)
            if obj.isPresMain == true
                obj.tri6Map = tri6Map;
            else
                disp('Share not allowed! Mesh main physics is velocity');
                return
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Create map of nodes objects
        function setNodes(obj,pt)
            % Save all coordinates for vtk
            obj.vtkNodCor        = zeros(size(pt,1),3);
            obj.vtkNodCor(:,1:2) = pt(:,2:3);
            % Set nodes objects
            for i = 1:size(pt,1)
                obj.nodMap(pt(i,1)) = node(pt(i,1),pt(i,2:3));
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Create map of velocity Dirichlet BCs
        function setBC(obj)
            % Get number of BCs
            nBCs = size(obj.bcVMap,1);
            key  = keys(obj.bcVMap);
            for i = 1:nBCs
                % Get BC iter
                iterBC = obj.bcVMap(key{i});
                nodID  = iterBC.nodId;
                % Get node with ID
                iterNod = obj.nodMap(nodID);
                % Set BC val and flag;
                iterNod.bcVal    = iterBC.bcVal;
                iterNod.bcFlag   = iterBC.bcFlag;
                iterNod.hasDirBC = true;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Set Dirichlet BCs
        function readBC(obj,bcs)
            for i = 1:size(bcs,1)
                obj.bcVMap(bcs(i,1)) = bcinfo(bcs(i,1),bcs(i,2:4),bcs(i,5:7));
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Set Poiseuille profile at equilibrium as initial condition
        function setPoiseuilleProf(obj,umax,H)
            % Get number of BCs
            nBCs = size(obj.bcVMap,1);
            key  = keys(obj.bcVMap);
            for i = 1:nBCs
                % Get BC iter
                iterBC = obj.bcVMap(key{i});
                nodID  = iterBC.nodId;
                % Get node with ID
                iterNod = obj.nodMap(nodID);
                % Get nodeBc and bcFlag
                iterNodBcVal  = iterNod.bcVal;
                iterNodBcFlag = iterNod.bcFlag;
                if isequal(iterNodBcFlag,[1 1]) && isequal(iterNodBcVal,[1 0])
                    % Calculate Poiseuille profile
                    iterNodCoor = iterNod.coor;
                    % Get y-coordinate
                    y = iterNodCoor(2);
                    % Calculate profile
                    perturb = 0.0001;
                    ux      = -(umax/H^2)*(y-perturb)^2+umax;
                    % Set new BC val
                    iterNod.bcVal = [ux 0];
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Create tri3 elements
        function setTri3(obj,mat,conn,intRule)
            % Save connectivity for vtk
            obj.vtkTri3Con        = zeros(size(conn,1),3);
            obj.vtkTri3Con(:,1:3) = conn(:,2:4)-1;
            % Set flags
            %obj.isPresMain = true;
            % Create elements
            for i = 1:size(conn,1)
                obj.tri3Map(conn(i,1)) = tri3(...
                    conn(i,1),...
                    conn(i,2:4),...
                    mat,...
                    obj.nodMap(conn(i,2)),...
                    obj.nodMap(conn(i,3)),...
                    obj.nodMap(conn(i,4)),...
                    intRule...
                    );
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Create tri6 elements
        function setTri6(obj,mat,conn,intRule)
            % Set flags
            obj.isVeloMain = true;
            % Create elements
            for i = 1:size(conn,1)
                obj.tri6Map(conn(i,1)) = tri6(...
                    conn(i,1), ...
                    conn(i,2:7), ...
                    mat, ...
                    obj.nodMap(conn(i,2)),...
                    obj.nodMap(conn(i,3)),...
                    obj.nodMap(conn(i,4)),...
                    obj.nodMap(conn(i,5)),...
                    obj.nodMap(conn(i,6)),...
                    obj.nodMap(conn(i,7)),...
                    intRule...
                    );
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Set DoF for velocity problem
        function setDoFV(obj)
            nEle = size(obj.tri3Map,1);
            key  = keys(obj.tri3Map);
            for i = 1:nEle
                iterTri3 = obj.tri3Map(key{i});
                % Node 1
                obj.nextDoF = iterTri3.n1.setDoFTri3V(obj.nextDoF);
                % Node 2
                obj.nextDoF = iterTri3.n2.setDoFTri3V(obj.nextDoF);
                % Node 3
                obj.nextDoF = iterTri3.n3.setDoFTri3V(obj.nextDoF);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Set DoF for pressure problem
        function setDoFP(obj)
            nEle = size(obj.tri3Map,1);
            key  = keys(obj.tri3Map);
            for i = 1:nEle
                iterTri3 = obj.tri3Map(key{i});
                % Node 1
                obj.nextDoF = iterTri3.n1.setDoFTri3P(obj.nextDoF);
                % Node 2
                obj.nextDoF = iterTri3.n2.setDoFTri3P(obj.nextDoF);
                % Node 3
                obj.nextDoF = iterTri3.n3.setDoFTri3P(obj.nextDoF);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Set DoF for velocity AND pressure problem
        function setDoFVP(obj)
            nEle = size(obj.tri3Map,1);
            key  = keys(obj.tri3Map);
            for i = 1:nEle
                iterTri3 = obj.tri3Map(key{i});
                % Node 1
                obj.nextDoF = iterTri3.n1.setDoFTri3VP(obj.nextDoF);
                % Node 2
                obj.nextDoF = iterTri3.n2.setDoFTri3VP(obj.nextDoF);
                % Node 3
                obj.nextDoF = iterTri3.n3.setDoFTri3VP(obj.nextDoF);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Fix tri6 mesh nodes (4,5,6) coordinates
        function fixTri6Mesh(obj)
            nEle = size(obj.tri6Map,1);
            key  = keys(obj.tri6Map);
            for i = 1:nEle
                iterTri6 = obj.tri6Map(key{i});
                coor1 = iterTri6.n1.coor;
                coor2 = iterTri6.n2.coor;
                coor3 = iterTri6.n3.coor;
                % The coordinate of nodes below is defective GiD
                % Calculate updated coordinates
                coor4 = (coor1+coor2)/2;
                coor5 = (coor2+coor3)/2;
                coor6 = (coor1+coor3)/2;
                % Update object
                iterTri6.n4.updNodCoor(coor4(1,1),coor4(1,2));
                iterTri6.n5.updNodCoor(coor5(1,1),coor5(1,2));
                iterTri6.n6.updNodCoor(coor6(1,1),coor6(1,2));
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Assemble Mv global
        function M = assembleMv(obj)
            % Initialize
            nDoF  = obj.nextDoFV-1;
            M = zeros(nDoF,nDoF);
            % Loop over Tri6 elements
            nEle = size(obj.tri6Map,1);
            keyE  = keys(obj.tri6Map);
            for i = 1:nEle
                iterTri6 = obj.tri6Map(keyE{i});
                % Element contribution to M
                Me = iterTri6.getMe();
                %%% Assemble Mv using connectivity %%%
                nNod = size(iterTri6.nodMap,1);
                conn = iterTri6.conn;
                for j = 1:nNod
                    iterNod_j = iterTri6.nodMap(conn(1,j));
                    gDoF_j    = iterNod_j.gDoF;
                    lDof_j    = [2*j-1 2*j];
                    % Pupulate global M
                    for k = 1:nNod
                        iterNod_k = iterTri6.nodMap(conn(1,k));
                        gDoF_k    = iterNod_k.gDoF;
                        lDof_k    = [2*k-1 2*k];
                        % Zero-out if entry < tol
                        tol = 10e-15;
                        tempMe = Me(lDof_j,lDof_k);
                        tempMe ( abs(tempMe) < tol ) = 0;
                        % Assign values
                        M(gDoF_j,gDoF_k) = M(gDoF_j,gDoF_k)+tempMe;
                    end
                end
            end
            % Make it sparse
            M = sparse(M);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Assemble Kv global
        function K = assembleKv(obj)
            % Initialize
            nDoF  = obj.nextDoFV-1;
            K = zeros(nDoF,nDoF);
            % Loop over Tri6 elements
            nEle = size(obj.tri6Map,1);
            keyE  = keys(obj.tri6Map);
            for i = 1:nEle
                iterTri6 = obj.tri6Map(keyE{i});
                % Element contribution to K
                Ke = iterTri6.getKe();
                %%% Assemble Kv using connectivity %%%
                nNod = size(iterTri6.nodMap,1);
                conn = iterTri6.conn;
                for j = 1:nNod
                    iterNod_j = iterTri6.nodMap(conn(1,j));
                    gDoF_j    = iterNod_j.gDoF;
                    lDof_j    = [2*j-1 2*j];
                    % Pupulate global K
                    for k = 1:nNod
                        iterNod_k = iterTri6.nodMap(conn(1,k));
                        gDoF_k    = iterNod_k.gDoF;
                        lDof_k    = [2*k-1 2*k];
                        % Zero-out if entry < tol
                        tol = 10e-15;
                        tempKe = Ke(lDof_j,lDof_k);
                        tempKe ( abs(tempKe) < tol ) = 0;
                        % Assign values
                        %K(gDoF_j,gDoF_k) = K(gDoF_j,gDoF_k)+Ke(lDof_j,lDof_k);
                        K(gDoF_j,gDoF_k) = K(gDoF_j,gDoF_k)+tempKe;
                    end
                end
            end
            % Make it sparse
            K = sparse(K);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Assemble Kelas global
        function K = assembleKelas(obj)
            % Initialize
            nDoF  = obj.nextDoFV-1;
            K = zeros(nDoF,nDoF);
            % Loop over Tri6 elements
            nEle = size(obj.tri6Map,1);
            keyE  = keys(obj.tri6Map);
            for i = 1:nEle
                iterTri6 = obj.tri6Map(keyE{i});
                % Element contribution to K
                Ke = iterTri6.getKeElas();
                %%% Assemble Kv using connectivity %%%
                nNod = size(iterTri6.nodMap,1);
                keyN = keys(iterTri6.nodMap);
                for j = 1:nNod
                    iterNod_j = iterTri6.nodMap(keyN{j});
                    gDoF_j    = iterNod_j.gDoF;
                    lDof_j    = [2*j-1 2*j];
                    % Pupulate global K
                    for k = 1:nNod
                        iterNod_k = iterTri6.nodMap(keyN{k});
                        gDoF_k    = iterNod_k.gDoF;
                        lDof_k    = [2*k-1 2*k];
                        % Zero-out if entry < tol
                        tol = 10e-15;
                        tempKe = Ke(lDof_j,lDof_k);
                        tempKe ( abs(tempKe) < tol ) = 0;
                        % Assign values
                        K(gDoF_j,gDoF_k) = K(gDoF_j,gDoF_k)+tempKe;
                    end
                end
            end
            % Make it sparse
            K = sparse(K);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Assemble Q global
        function Q = assembleQ(obj,Ut)
            % Initialize
            nDoF  = obj.nextDoFV-1;
            Q = zeros(nDoF,nDoF);
            % Loop over Tri6 elements
            nEle = size(obj.tri6Map,1);
            keyE = keys(obj.tri6Map);
            for i = 1:nEle
                iterTri6 = obj.tri6Map(keyE{i});
                % Element contribution to Q
                Qe = iterTri6.getQe(Ut);
                %%% Assemble Q using connectivity %%%
                nNod = size(iterTri6.nodMap,1);
                conn = iterTri6.conn;
                for j = 1:nNod
                    iterNod_j = iterTri6.nodMap(conn(1,j));
                    gDoF_j    = iterNod_j.gDoF;
                    lDof_j    = [2*j-1 2*j];
                    % Pupulate global Q
                    for k = 1:nNod
                        iterNod_k        = iterTri6.nodMap(conn(1,k));
                        gDoF_k           = iterNod_k.gDoF;
                        lDof_k           = [2*k-1 2*k];
                        % Zero-out if entry < tol
                        tol = 10e-15;
                        tempQe = Qe(lDof_j,lDof_k);
                        tempQe ( abs(tempQe) < tol ) = 0;
                        % Assign values
                        Q(gDoF_j,gDoF_k) = Q(gDoF_j,gDoF_k)+tempQe;
                    end
                end
            end
            % Make it sparse
            Q = sparse(Q);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Use direct method for velocity BC
        function [lhs,rhs] = directMethodBCVel(obj,lhs,rhs)
            % Get number of BCs
            nBCs = size(obj.bcVMap,1);
            key  = keys(obj.bcVMap);
            % Initialize vector of global DOF with BCs
            dofBC = [];
            dofVal = [];
            % Get global DOF with BCs
            for i = 1:nBCs
                % Get BC iter
                iterBC = obj.bcVMap(key{i});
                nodID  = iterBC.nodId;
                % Get node with ID
                iterNod = obj.nodMap(nodID);
                % Get global DoF
                gDoF = iterNod.gDoF;
                % Get BC val and flag
                bcVal  = iterNod.bcVal;
                bcFlag = iterNod.bcFlag;
                for j = 1:2
                    if bcFlag(j) == 1
                        % Save global DoF with BC
                        dofBC = [dofBC; gDoF(j)];
                        % Save BC value
                        dofVal = [dofVal; bcVal(j)];
                        % Modify rhs
                        rhs = rhs-bcVal(1,j)*lhs(:,gDoF(1,j));
                    end
                end
            end
            % Set Dirichlet BCs in rhs
            for i = 1:length(dofBC)
                rhs(dofBC(i),1) = dofVal(i);
            end
            % Set to zero rhs rows with BC
            lhs(dofBC,:) = 0;
            % Set to zero rhs rows with BC
            lhs(:,dofBC) = 0;
            % Set to one diagonal entry
            for i = 1:length(dofBC)
                lhs(dofBC(i),dofBC(i)) = 1;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Use direct method for pressure BC
        function [lhs,rhs] = directMethodBCPre(obj,lhs,rhs)
            % Get number of BCs
            nBCs = size(obj.bcVMap,1);
            key  = keys(obj.bcVMap);
            % Initialize vector of global DOF with BCs
            dofBC = [];
            dofVal = [];
            % Get global DOF with BCs
            for i = 1:nBCs
                % Get BC iter
                iterBC = obj.bcVMap(key{i});
                nodID  = iterBC.nodId;
                % Get node with ID
                iterNod = obj.nodMap(nodID);
                % Get global DoF
                gDoF = iterNod.gDoF;
                % Get BC val and flag
                bcVal  = iterNod.bcVal;
                bcFlag = iterNod.bcFlag;
                % Check only the FIRST entry - IGNORE second
                for j = 1:1
                    if bcFlag(j) == 1
                        % Save global DoF with BC
                        dofBC = [dofBC; gDoF(j)];
                        % Save BC value
                        dofVal = [dofVal; bcVal(j)];
                        % Modify rhs
                        rhs = rhs-bcVal(1,j)*lhs(:,gDoF(1,j));
                    end
                end
            end
            % Set Dirichlet BCs in rhs
            for i = 1:length(dofBC)
                rhs(dofBC(i),1) = dofVal(i);
            end
            % Set to zero rhs rows with BC
            lhs(dofBC,:) = 0;
            % Set to zero rhs rows with BC
            lhs(:,dofBC) = 0;
            % Set to one diagonal entry
            for i = 1:length(dofBC)
                lhs(dofBC(i),dofBC(i)) = 1;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Enforce Dirichlet boundary conditions
        function rhs = setDirichBC(obj,rhs)
            % Get number of BCs
            nBCs = size(obj.bcVMap,1);
            key  = keys(obj.bcVMap);
            % Initialize vector of global DOF with BCs
            dofBC  = [];
            dofVal = [];
            % Get global DOF with BCs
            for i = 1:nBCs
                % Get BC iter
                iterBC = obj.bcVMap(key{i});
                nodID  = iterBC.nodId;
                % Get node with ID
                iterNod = obj.nodMap(nodID);
                % Get global DoF
                gDoF = iterNod.gDoF;
                % Get BC val and flag
                bcVal  = iterNod.bcVal;
                bcFlag = iterNod.bcFlag;
                for j = 1:2
                    if bcFlag(j) == 1
                        % Save global DoF with BC
                        dofBC = [dofBC; gDoF(j)];
                        % Save BC value
                        dofVal = [dofVal; bcVal(j)];
                    end
                end
            end
            % Enforce Dirichlet BC
            rhs(dofBC,1) = dofVal;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Assemble Kp global
        function K = assembleKp(obj)
            % Initialize
            nDoF  = obj.nextDoFP-1;
            K     = zeros(nDoF,nDoF);
            % Loop over Tri6 elements
            nEle = size(obj.tri3Map,1);
            keyE  = keys(obj.tri3Map);
            for i = 1:nEle
                iterTri3 = obj.tri3Map(keyE{i});
                % Element contribution to Kp
                Ke = iterTri3.getKe();
                %%% Assemble Kp using connectivity %%%
                nNod = size(iterTri3.nodMap,1);
                conn = iterTri3.conn;
                for j = 1:nNod
                    iterNod_j = iterTri3.nodMap(conn(1,j));
                    gDoF_j    = iterNod_j.gDoF;
                    lDof_j    = j;
                    % Pupulate global Kp
                    for k = 1:nNod
                        iterNod_k        = iterTri3.nodMap(conn(1,k));
                        gDoF_k           = iterNod_k.gDoF;
                        lDof_k           = k;
                        % Zero-out if entry < tol
                        tol = 10e-15;
                        tempKe = Ke(lDof_j,lDof_k);
                        tempKe ( abs(tempKe) < tol ) = 0;
                        % Assign values
                        K(gDoF_j,gDoF_k) = K(gDoF_j,gDoF_k)+tempKe;
                    end
                end
            end
            % Make it sparse
            K = sparse(K);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Evaluate divergence of Us
        function divUs = evalDivUs(obj,Us)
            % Initialize
            nDoF  = obj.nextDoFP-1;
            divUs = zeros(nDoF,1);
            % Loop over Tri3 elements
            nEle = size(obj.tri3Map,1);
            keyE = keys(obj.tri3Map);
            for i = 1:nEle
                % Iterator Tri3
                iterTri3 = obj.tri3Map(keyE{i});
                % Calc divergence of U using Tri6 element
                iterTri6 = obj.tri6Map(iterTri3.ID);
                divUsE   = iterTri3.evalDivUe(Us,iterTri6);
                %%% Assemble divUs using connectivity %%%
                nNod = size(iterTri3.nodMap,1);
                conn = iterTri3.conn;
                for j = 1:nNod
                    iterNod = iterTri3.nodMap(conn(1,j));
                    gDoF    = iterNod.gDoF;
                    lDof    = j;
                    % Zero-out if entry < tol
                    tol = 10e-15;
                    tempdivUsE = divUsE(lDof,1);
                    tempdivUsE ( abs(tempdivUsE) < tol ) = 0;
                    % Pupulate global gradUs
                    divUs(gDoF,1) = divUs(gDoF,1)+tempdivUsE;
                end
            end
            divUs = sparse(divUs);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Evaluate gradient of pressure
        function gradP = evalGradP(obj,p)
            % Initialize
            nDoF  = obj.nextDoFV-1;
            gradP = zeros(nDoF,1);
            % Loop over Tri6 elements
            nEle = size(obj.tri6Map,1);
            keyE = keys(obj.tri6Map);
            for i = 1:nEle
                % Iterator Tri6
                iterTri6 = obj.tri6Map(keyE{i});
                % Calc grad(p) using Tri3 element
                iterTri3 = obj.tri3Map(iterTri6.ID);
                gradPe   = iterTri6.evalGradPe(p,iterTri3);
                %%% Assemble gradP using connectivity %%%
                nNod = size(iterTri6.nodMap,1);
                conn = iterTri6.conn;
                for j = 1:nNod
                    iterNod = iterTri6.nodMap(conn(1,j));
                    gDoF    = iterNod.gDoF;
                    lDof    = j;
                    % Zero-out if entry < tol
                    tol = 10e-15;
                    tempgradPe = gradPe([2*lDof-1 2*lDof],1);
                    tempgradPe ( abs(tempgradPe) < tol ) = 0;
                    % Pupulate global gradP
                    gradP(gDoF,1) = gradP(gDoF,1)+tempgradPe;
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot normalized velocity vector field
        function plotU(obj,U)
            % Initialize
            nDoF = (2/3)*(obj.nextDoF-1);
            X    = zeros(nDoF,1);
            Y    = zeros(nDoF,1);
            Ux   = zeros(nDoF,1);
            Uy   = zeros(nDoF,1);
            % Loop over nodes
            nNod = size(obj.nodMap,1);
            keyN = keys(obj.nodMap);
            for i = 1:nNod
                iterNod = obj.nodMap(keyN{i});
                coor    = iterNod.coor;
                gDoF    = iterNod.gDoF;
                X(i,1)  = coor(1,1);
                Y(i,1)  = coor(1,2);
                Vnorm   = norm([U(gDoF(1),1) U(gDoF(2),1)]);
                Ux(i,1) = U(gDoF(1),1);%/Vnorm;
                Uy(i,1) = U(gDoF(2),1);%/Vnorm;
            end
            % Plot
            quiver(X,Y,Ux,Uy,0.5);
            xlabel('x')
            ylabel('y')
            grid minor
%             pbaspect([4 1 1])
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot curl of velocity vector field
        function plotUCurl(obj,U)
            % Initialize
            nEle  = size(obj.tri6Map,1);
            curlU = [];
            % Loop elements
            keyN = keys(obj.tri6Map);
            for i = 1:nEle
                iterTri6 = obj.tri6Map(keyN{i});
                % Get curlU at element int points
                curlUeIntPt = iterTri6.evalCurlUeIntPt(U);
                % Store data
                curlU = [curlU; curlUeIntPt];
            end
            % Plot
            X    = curlU(:,1);
            Y    = curlU(:,2);
            curl = curlU(:,3);
            c = linspace(min(curl),max(curl),length(X));
            scatter3(X,Y,curl,10,c,'filled');
            colormap("jet")
            colorbar
            view(2)
            xlabel('x')
            ylabel('y')
            pbaspect([4 1 1])
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot velocity magnitude scalar field
        function [X,Y,Unorm] = plotUMag(obj,U)
            % Initialize
            nDoF  = obj.nextDoFV-1;
            X     = zeros(nDoF,1);
            Y     = zeros(nDoF,1);
            Unorm = zeros(nDoF,1);
            % Loop over nodes
            nNod = size(obj.nodMap,1);
            keyN = keys(obj.nodMap);
            for i = 1:nNod
                iterNod = obj.nodMap(keyN{i});
                coor    = iterNod.coor;
                gDoF    = iterNod.gDoF;
                X(i,1)  = coor(1,1);
                Y(i,1)  = coor(1,2);
                Unorm(i,1) = norm([U(gDoF(1),1) U(gDoF(2),1)]);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot pressure field
        function plotP(obj,p)
            % Initialize
            nDoF = obj.nextDoFP-1;
            X    = zeros(nDoF,1);
            Y    = zeros(nDoF,1);
            P    = zeros(nDoF,1);
            % Loop over nodes
            nNod = size(obj.nodMap,1);
            keyN = keys(obj.nodMap);
            j    = 1;
            for i = 1:nNod
                iterNod = obj.nodMap(keyN{i});
                % Only use pressure nodes info
                if iterNod.isPNod == true
                    coor    = iterNod.coor;
                    gDoF    = iterNod.gDoF;
                    if gDoF(1) ~= -1
                        X(j,1)  = coor(1,1);
                        Y(j,1)  = coor(1,2);
                        P(j,1)  = p(gDoF(1),1);
                        j       = j+1;
                    end
                end
            end
            % Plot
            scatter3(X,Y,P);
            xlabel('x')
            ylabel('y')
            pbaspect([4 1 1])
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Eval body force term for difussion manufactured solution
        function F = evalMyBF(obj)
            % Initialize
            nDoF = obj.nextDoFV-1;
            F    = zeros(nDoF,1);
            % Loop over Tri6 elements
            nEle = size(obj.tri6Map,1);
            keyE  = keys(obj.tri6Map);
            for i = 1:nEle
                iterTri6 = obj.tri6Map(keyE{i});
                % Element contribution to BF
                bFe = iterTri6.getMyBF();
                %%% Assemble F using connectivity %%%
                nNod = size(iterTri6.nodMap,1);
                conn = iterTri6.conn;
                for j = 1:nNod
                    iterNod = iterTri6.nodMap(conn(1,j));
                    gDoF    = iterNod.gDoF;
                    lDof    = [2*j-1 2*j];
                    % Pupulate global F
                    F(gDoF,1) = F(gDoF,1)+bFe(lDof,1);
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Export velocity field to GiD
        function toGidV(obj,fileName,step,U)
            % GiD file name
            % fileName = 'NS_VelField';
            % Get node map size
            nNod = size(obj.nodMap,1);
            % Loop over nodes
            coor = zeros(nNod,3);
            vel  = zeros(2*nNod,1);
            for i = 1:nNod
                iterNod   = obj.nodMap(i);
                % Store coordinate
                coor(i,:) = [iterNod.iD iterNod.coor];
                % Store velocity
                gDoF = iterNod.gDoF;
                vel(i*2-1:i*2,1) = U(gDoF);
            end
            % Get element map
            nEle = size(obj.tri6Map,1);
            % Loop over elemens
            conn = zeros(nEle,7);
            for i = 1:nEle
                iterTri6 = obj.tri6Map(i);
                % Store connectivity
                conn(i,:) = [iterTri6.ID iterTri6.conn];
            end
            % Export to GiD
            ToGiD_v1_3(fileName,step,coor,conn,vel);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Export pressure field to GiD
        function toGidP(obj,fileName,step,P)
            % GiD file name
            % fileName = 'NS_PreField';
            % Get node map size
            nNod = size(obj.nodMap,1);
            % Loop over nodes
            coor = [];
            pre  = [];
            for i = 1:nNod
                iterNod   = obj.nodMap(i);
                if iterNod.isDofSet ~= false
                    % Store coordinate
                    coor = [coor; [iterNod.iD iterNod.coor]];
                    % Store pressure
                    gDoF = iterNod.gDoF;
                    pre  = [pre; P(gDoF)];
                end
            end
            % Get element map
            nEle = size(obj.tri3Map,1);
            % Loop over elemens
            conn = zeros(nEle,4);
            for i = 1:nEle
                iterTri3 = obj.tri3Map(i);
                % Store connectivity
                conn(i,:) = [iterTri3.ID iterTri3.conn];
            end
            % Export to GiD
            ToGiD_v1_3(fileName,step,coor,conn,pre);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Assemble residual
        function Res = assembleRes(obj,Uh,Uhk,iTime)
            % Initialize
            nDoF = obj.nextDoF-1;
            Res  = zeros(nDoF,1);
            % Loop over Tri3 elements
            nEle = size(obj.tri3Map,1);
            keyE = keys(obj.tri3Map);
            for i = 1:nEle
                iterTri3 = obj.tri3Map(keyE{i});
                % Element contribution to Residual
                Res_e  = iterTri3.evalResE(Uh,Uhk);
                %%% Assemble global residual using connectivity %%%
                nNod = size(iterTri3.nodMap,1);
                conn = iterTri3.conn;
                for j = 1:nNod
                    iterNod_j = iterTri3.nodMap(conn(1,j));
                    gDoF_j    = iterNod_j.gDoF;
                    lDof_j    = [j , j+3 , j+6];
                    % Populate global Res
                    Res(gDoF_j,1) = Res(gDoF_j,1)+Res_e(lDof_j,1);
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Assemble consstent Jacobian matrix
        function KT = assembleKT(obj,Uh,Uhk,iTime)
            % Initialize
            nDoF = obj.nextDoF-1;
            KT   = zeros(nDoF,nDoF);
            % Loop over Tri3 elements
            nEle = size(obj.tri3Map,1);
            keyE = keys(obj.tri3Map);
            for i = 1:nEle
                iterTri3 = obj.tri3Map(keyE{i});
                % Element contribution to KT
                KTe  = iterTri3.evalKTe(Uh,Uhk);
                %%% Assemble global consistent tangent using connectivity %%%
                nNod = size(iterTri3.nodMap,1);
                conn = iterTri3.conn;
                for j = 1:nNod
                    iterNod_j = iterTri3.nodMap(conn(1,j));
                    gDoF_j    = iterNod_j.gDoF;
                    lDof_j    = [j , j+3 , j+6];
                    % Populate global KT
                    for k = 1:nNod
                        iterNod_k = iterTri3.nodMap(conn(1,k));
                        gDoF_k    = iterNod_k.gDoF;
                        lDof_k    = [k , k+3 , k+6];
                        % Zero-out if entry < tol
                        tol = 10e-15;
                        %
                        tempKTe = KTe(lDof_j,lDof_k);
                        tempKTe ( abs(tempKTe) < tol ) = 0;
                        % Assign values
                        KT(gDoF_j,gDoF_k)  = KT(gDoF_j,gDoF_k)+tempKTe;
                    end
                end
            end
            % Make it sparse
            KT = sparse(KT);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function KT = assembleKTFN(obj,Res,Uhk,h,iTime)
            % Initialize
            nDoF = obj.nextDoF-1;
            KT   = zeros(nDoF,nDoF);
            for j = 1:nDoF
                % Unit vector
                e      = zeros(nDoF,1);
                e(j,1) = 1;
                % Finite difference
                Uhk_h   = Uhk+(h*e);
                Res_h   = obj.assembleRes(Uhk,Uhk_h,iTime);
                KT(:,j) = (1/h)*( Res_h-Res );
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Use direct method for BC
        function [lhs,rhs] = directMethodBC(obj,lhs,rhs)
            % Get number of BCs
            nBCs = size(obj.bcVMap,1);
            key  = keys(obj.bcVMap);
            % Initialize vector of global DOF with BCs
            dofBC  = [];
            dofVal = [];
            % Get global DOF with BCs
            for i = 1:nBCs
                % Get BC iter
                iterBC = obj.bcVMap(key{i});
                nodID  = iterBC.nodId;
                % Get node with ID
                iterNod = obj.nodMap(nodID);
                % Get global DoF
                gDoF = iterNod.gDoF;
                % Get BC val
                bcVal = iterNod.bcVal;
                bcVal = 0*bcVal;
                % JAL: Assumed that all nodes with constrained DoF
                %      have a zero applied velocity. In this case, 
                %      the corresponding entry in the incremental 
                %      solution vector will have a value equal to 0
                %      and the solution will be computed as Uhk+Ug
                % Get BC flag
                bcFlag = iterNod.bcFlag;
                % Check only the FIRST TWO entries (velocity BCs)
                for j = 1:2
                    if bcFlag(j) == 1
                        % Save global DoF with BC
                        dofBC = [dofBC; gDoF(j)];
                        % Save BC value
                        dofVal = [dofVal; bcVal(j)];
                        % Modify rhs
                        rhs = rhs-bcVal(1,j)*lhs(:,gDoF(1,j));
                    end
                end
            end
            % Set Dirichlet BCs in rhs
            for i = 1:length(dofBC)
                rhs(dofBC(i),1) = dofVal(i);
            end
            % Set equal to zero lhs rows with BC
            lhs(dofBC,:) = 0;
            % Set equal to zero lhs columns with BC
            lhs(:,dofBC) = 0;
            % Set equal to one diagonal entry with BC
            for i = 1:length(dofBC)
                lhs(dofBC(i),dofBC(i)) = 1;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end

