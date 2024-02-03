classdef nonlintrans < handle
    %NONLINTRANS Summary of this class goes here
    %   Implementation of Crank-Nicolson for
    %   transient convection-difussion equation
    %   with SUPG (Streamline upwind Petrov–Galerkin)
    
    properties
        % Input
        mesh;
        mat;

        % M-mass matrix
        M;
        isMComp = false;
        % K-diffusion matrix
        K;
        isKComp = false;
        % KS-diffusion matrix
        KS;
        isKSComp = false;
        % Q-convection matrix
        Q;
        isQComp = false;
        % QS-convection matrix
        QS;
        isQSComp = false;

        % Fint := a_supg
        Fint;
        % Fext := F_supg
        Fext;
        % Residual := Fint-Fext
        Res;
        Res0;
        Res_p;

        % Tangent matrix
        KT;
        isKTComp = false;
        % Increment solution
        dUhk;
        % Increment residual
        dRes;
        % Solution at 'n'-step
        Uh;
        % Solution at 'n+1'-step & 'k'-iteration
        Uhk;

        % Essential BCs
        Ug;

        % Last solved step
        lastSlvStp = -1;
        % All solution container
        Uall;

        % Tolerance for nonlinear solver
        tolA  = 1e-04;
        tolR  = 1e-02;
        toldR = 5e-04;

        % Max number of NR iterations
        maxIter = 30;
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        function obj = nonlintrans(mesh,mat)
            % Constructor
            obj.mesh = mesh;
            obj.mat  = mat;
            
            % Set essential BC - ONLY velocity BCs
            obj.Ug = zeros(obj.mesh.nextDoF-1,1);
            nBCs   = size(obj.mesh.bcVMap,1);
            key    = keys(obj.mesh.bcVMap);
            for i = 1:nBCs
                % Get BC iter
                iterBC = obj.mesh.bcVMap(key{i});
                nodID  = iterBC.nodId;
                % Get node with ID
                iterNod     = obj.mesh.nodMap(nodID);
                iterNodGDof = iterNod.gDoF;
                % Set BC val using flag val
                bcFlag = iterBC.bcFlag;
                bcVal  = iterBC.bcVal;
                for j = 1:3
                    if bcFlag(j) == 1
                        obj.Ug(iterNodGDof(j)) = bcVal(j);
                    end
                end
            end

            % Set initial condition - step 0
            obj.Uh  = obj.Ug;
            
            % Set initial residual
            obj.Res   = 0*(obj.Ug);
            obj.Res_p = 0*(obj.Ug);
            % Set initial guess for incremental sol vector
            obj.dUhk = 0*(obj.Ug);
            % Set current NR sol vector
            obj.Uhk  = obj.Uh+obj.dUhk;
            
            % Set last solved step
            obj.lastSlvStp = 0;
            % Set solution container
            nstp          = obj.mat.nstp;
            obj.Uall      = zeros(size(obj.Ug,1),nstp+1);
            obj.Uall(:,1) = obj.Uh;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Solve
        function solve(obj)
            nSteps = obj.mat.nstp;
            for iStep = 1:nSteps
                obj.solveStep(iStep);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Solve step
        function solveStep(obj,iStep)
            % Current time
            iTime = (iStep/obj.mat.nstp)*(obj.mat.tf);
            for k = 0:obj.maxIter
                % Compute current residual
                obj.Res  = obj.mesh.assembleRes(obj.Uh,obj.Uhk,iTime);
                obj.dRes = obj.Res-obj.Res_p;
                norm(obj.Res);
                % Set initial residual
                if k == 0
                    obj.Res0 = obj.Res;
                end
                % Check residual
                if ( norm(obj.Res) <= obj.tolA || ...
                     norm(obj.Res) <= obj.tolR*norm(obj.Res0) || ...
                     norm(obj.dRes) <= obj.toldR || ...
                     k == obj.maxIter )
                    % Current solution Uhk satisfies residual or maxIter
                    % Set Uh = Uhk and proceed to next iStep
                    obj.Uh = obj.Uhk;
                    %
                    obj.Uall(:,iStep+1) = obj.Uh;
                    %
                    obj.Res_p = obj.Res;
                    %

                    % ::: START POSTPROCESSING AT EACH STEP
                    %obj.mesh.plotU(obj.Uhk)
                    path = strcat('../paraview/sol_step_',num2str(iStep),'.vtk');
                    obj.print2VTK(path,obj.Uh);
                    % ::: END START POSTPROCESSING AT EACH STEP

                    % Print convergence message
                    fprintf(strcat('\nStep Converged: Step=',num2str(iStep),'\n'));
                    fprintf(strcat('itrNR=',num2str(k),'\n'));
                    fprintf(strcat('|Res0|=',num2str(norm(obj.Res0)),'\n'));
                    fprintf(strcat('|Res|=',num2str(norm(obj.Res)),'\n'));
                    fprintf(strcat('|Res|/|Res0|=',num2str(norm(obj.Res)/norm(obj.Res0)),'\n'));
                    fprintf(strcat('|dRes|=',num2str(norm(obj.dRes)),'\n'));
                    %
                    return
                else
                    % Solve nonlinear system using NR
                    if k == 0
                        % Compute Jacobian initial approximation using 
                        % finite difference (too expensive)
                        %h = 1e-02;
                        %obj.KT = obj.mesh.assembleKTFN(obj.Res,obj.Uhk,h,iTime);
                        % Use idendity matrix
                        %obj.KT = eye(obj.mesh.nextDoF-1);
                        % Compute consistent Jacobian
                        obj.KT  = obj.mesh.assembleKT(obj.Uh,obj.Uhk,iTime);
                    else
                        % Compute Jacobian approximation using Broyden
                        %obj.KT = obj.KT + ...
                        %         ( obj.dRes - obj.KT*obj.dUhk ) * ...
                        %         ( (obj.dUhk)' / (norm(obj.dUhk))^2 );
                        % Compute consistent Jacobian
                        %obj.KT  = obj.mesh.assembleKT(obj.Uh,obj.Uhk,iTime);
                    end
                    % LHS and RHS
                    lhs = obj.KT;
                    rhs = obj.Res;
                    % Direct method for BCs
                    [lhs,rhs] = obj.mesh.directMethodBC(lhs,rhs);
                    % Solve
                    obj.dUhk = lhs\rhs;
                    % Update current NR solution
                    obj.Uhk = obj.Uhk - obj.dUhk;
                    % Update residuals
                    obj.Res_p = obj.Res;
                    % ::: Start plot :::
                    %obj.mesh.plotU(obj.Uhk)
                    % ::: End plot :::
                end
            end
            % End of NR iterations for step iStep
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Print solution to vtk format
        function print2VTK(obj,filename,Uh)
            % Get node info
            nodCor  = obj.mesh.vtkNodCor;
            nNodes  = size(nodCor,1);
            % Get element info
            tri3Con = obj.mesh.vtkTri3Con;
            nEle    = size(tri3Con,1);
            % Get nodal velocity and pressure solution
            Uvtk = zeros(nNodes,3);
            Pvtk = zeros(nNodes,1);
            % Loop over nodes
            for i = 1:nNodes
                iterNod = obj.mesh.nodMap(i);
                gDoF    = iterNod.gDoF;
                % Velocity and pressure
                Uvtk(i,1:2) = Uh(gDoF(1:2));
                Pvtk(i,1) = Uh(gDoF(3));
            end

            % Open file
            fid = fopen(filename,'w');
            % VTK DataFile Version
            fprintf(fid,'# vtk DataFile Version 2.0\n');
            % Title
            fprintf(fid,'VTK from Matlab\n');
            % Format
            fprintf(fid,'ASCII\n');
            fprintf(fid,'DATASET POLYDATA\n');
            % Print node coordinates
            fprintf(fid,['POINTS ' num2str(nNodes) ' float\n']);
            fprintf(fid,'%d %d %d\n',nodCor');
            % Print elements
            fprintf(fid,'\nPOLYGONS %d %d\n',nEle,4*nEle);
            fprintf(fid,'3 %d %d %d\n',tri3Con');
            % PRINT PRESSURE FIELD
            fprintf(fid,['\nPOINT_DATA ' num2str(nNodes)]);
            fprintf(fid,'\nSCALARS P float 1');
            fprintf(fid,'\nLOOKUP_TABLE default\n');
            fprintf(fid,'%d %d %d\n',Pvtk);
            % PRINT VELOCITY FIELD
            fprintf(fid,'\nVECTORS Velocity float\n');
            fprintf(fid,'%d %d %d\n',Uvtk');
            % Close file
            fclose(fid);

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end

