classdef analysis < handle
    %ANALYSIS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Input
        meshP;
        meshV;
        mat;
        dt;
        % Velocity field
        Ut;
        Us;
        Utp1;
        % Pressure field
        p;
        % Force term for manufactured U solution
        myF;
        myU;
        % K matrix VELOCITY - only compute once
        Kv;
        isKvComp = false;
        % Q matrix - recalculated with each solution Ut
        Q;
        % M matrix - mass-like matrix - only computed once
        Mv;
        isMvComp = false;
        %isQComp = false;
        % K matrix PRESSURE - only compute once
        Kp;
        isKpComp = false;
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = analysis(meshP,meshV,mat,dt)
            % Constructor
            obj.meshP = meshP;
            obj.meshV = meshV;
            obj.mat   = mat;
            obj.dt    = dt;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Set initial solution 0 + BCs only for velocity
        function setInitialSol(obj)
            % System size
            gDoF = obj.meshV.nextDoFV-1;
            % Initial solution
            U0 = zeros(gDoF,1);
            % Apply BCs to initial solution
            nBCs = size(obj.meshV.bcVMap,1);
            key  = keys(obj.meshV.bcVMap);
            for i = 1:nBCs
                % Get BC iter
                iterBC = obj.meshV.bcVMap(key{i});
                nodID  = iterBC.nodId;
                % Get node with ID
                iterNod = obj.meshV.nodMap(nodID);
                % Get node global DoF
                gDoF = iterNod.gDoF;
                % Get nodal BC flags
                bcFlag = iterNod.bcFlag;
                % Get nodal BC vals
                bcVal = iterNod.bcVal;
                % Loop
                for j = 1:2
                    if bcFlag(j) == 1
                        U0(gDoF(j)) = bcVal(j);
                    end
                end
            end
            obj.Ut = sparse(U0);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Solve for tentative velocity
        function solveVStar(obj)
            % Assemble Mv and set flag equal to true
            if obj.isMvComp == false                
                obj.Mv       = obj.meshV.assembleMv();
                obj.isMvComp = true;
            end
            % Assemble Kv and set flag equal to true
            if obj.isKvComp == false                
                obj.Kv       = obj.meshV.assembleKv();
                obj.isKvComp = true;
            end
            % Assemble Q(t)
            obj.Q = obj.meshV.assembleQ(obj.Ut);
            % LHS and RHS
            lhs = obj.Mv+(obj.dt)*(obj.Kv);
            rhs = (obj.Mv)*(obj.Ut)-(obj.dt)*(obj.Q)*(obj.Ut);
            % Direct method for velocity
            [lhs,rhs] = obj.meshV.directMethodBCVel(lhs,rhs);
            % Solve
            obj.Us = lhs\rhs;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Solve for pressure
        function solvePress(obj)
            % Assemble Kp and set flag equal to true
            if obj.isKpComp == false                
                obj.Kp       = obj.meshP.assembleKp();
                obj.isKpComp = true;
            end
            % Calculate div(Us) = (grad /dot Us)
            divUs = obj.meshP.evalDivUs(obj.Us);
            % LHS and RHS
            lhs = obj.Kp;
            rhs = -(obj.mat.ro/obj.dt)*(divUs);
            % Direct method for velocity
            [lhs,rhs] = obj.meshP.directMethodBCPre(lhs,rhs);
            % Solve
            obj.p = lhs\rhs;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Correct velocities for incompresibility
        function solveVelIncomp(obj)
            % Assemble Mv and set flag equal to true
            if obj.isMvComp == false                
                obj.Mv       = obj.meshV.assembleMv();
                obj.isMvComp = true;
            end
            % Calculate grad(p)
            gradP = obj.meshV.evalGradP(obj.p);
            % LHS and RHS
            lhs = obj.Mv;
            rhs = (obj.Mv)*(obj.Us)-(obj.dt/obj.mat.ro)*(gradP);
            % Direct method for velocity
            [lhs,rhs] = obj.meshV.directMethodBCVel(lhs,rhs);
            % Solve
            obj.Utp1 = lhs\rhs;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot vector field
        function plotU(obj,U)
            obj.meshV.plotU(U);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot magnitude of vector field
        function plotUMag(obj,U)
            [X,Y,Unorm] = obj.meshV.plotUMag(U);
            % Plot
            scatter3(X,Y,log(Unorm));
            xlabel('x')
            ylabel('y')
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot curl of vector field
        function plotUCurl(obj,U)
            obj.meshV.plotUCurl(U);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot pressure field
        function plotP(obj,P)
            obj.meshP.plotP(P);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Create GiD postprocess file for velocity
        function toGidV(obj,fileName,step,U)
            obj.meshV.toGidV(fileName,step,U);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Create GiD postprocess file for pressure
        function toGidP(obj,fileName,step,P)
            obj.meshP.toGidP(fileName,step,P);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function solveDifussion(obj)
            % Get rhs
            rhs = obj.meshV.evalMyBF();
            obj.myF = rhs;
            % Assemble Kv and set flag equal to true           
            obj.Kv = obj.meshV.assembleKv();
            % Set lhs
            lhs = obj.Kv;
            % Direct method for velocity
            [lhs,rhs] = obj.meshV.directMethodBCVel(lhs,rhs);
            % Solve
            obj.myU = lhs\rhs;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end

