classdef tri3 < handle
    %TRI3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % ID
        ID = -1;
        % Material
        mat;
        % Element connectivity
        conn;
        % Nodes
        n1;
        n2;
        n3;
        nodMap;
        % Integration rule
        intRule;
        % Only compute Ke = B'B once since it does not change
        Ke;
        isKeComp = false;
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        function obj = tri3(iD,conn,mat,n1,n2,n3,intRule)
            % ID
            obj.ID  = iD;
            % Element connectivity
            obj.conn = conn;
            % Material
            obj.mat = mat;
            % Nodes
            obj.n1  = n1;
            obj.n2  = n2;
            obj.n3  = n3;
            % Integration rule
            obj.intRule = gaussint(intRule);
            % Element node map
            obj.nodMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
            obj.nodMap(n1.iD) = n1;
            obj.nodMap(n2.iD) = n2;
            obj.nodMap(n3.iD) = n3;
            % Indicate these are pressure nodes
            obj.n1.isPNod = true;
            obj.n2.isPNod = true;
            obj.n3.isPNod = true;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get coordinates of nodes
        function coor = getNodesCoor(obj)
            coor = zeros(3,2);
            %
            coor(1,:) = obj.n1.coor;
            coor(2,:) = obj.n2.coor;
            coor(3,:) = obj.n3.coor;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % To paramatric coordinates
        function [x,y] = doMap(obj,xi,eta)
            % Get nodal coordinates
            coor = getNodesCoor(obj);
            x1   = coor(1,1);
            x2   = coor(2,1);
            x3   = coor(3,1);
            y1   = coor(1,2);
            y2   = coor(2,2);
            y3   = coor(3,2);
            % Do mapping
            x = xi*x1+eta*x2+(1-xi-eta)*x3;
            y = xi*y1+eta*y2+(1-xi-eta)*y3;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % To physical coordinates
        function [xi,eta] = doInvMap(obj,x,y)
            % Get nodal coordinates
            coor = getNodesCoor(obj);
            x1   = coor(1,1);
            x2   = coor(2,1);
            x3   = coor(3,1);
            y1   = coor(1,2);
            y2   = coor(2,2);
            y3   = coor(3,2);
            % Do inverse mapping
            xi  = ( (y2-y3)*(x-x3) - (x2-x3)*(y-y3) )...
                / ( (x1-x3)*(y2-y3) - (y1-y3)*(x2-x3) );
            eta = ( -(y1-y3)*(x-x3) + (x1-x3)*(y-y3) )...
                / ( (x1-x3)*(y2-y3) - (y1-y3)*(x2-x3) );
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculate area
        function area = getA(obj)
            vec1 = obj.n2.coor-obj.n1.coor;
            vec2 = obj.n3.coor-obj.n1.coor;
            area = norm(cross([vec1 0],[vec2 0]))/2;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Determinant of Jacobian
        function detJ = getDetJ(obj)
            % Get nodal coordinates
            coor = getNodesCoor(obj);
            x1   = coor(1,1);
            x2   = coor(2,1);
            x3   = coor(3,1);
            y1   = coor(1,2);
            y2   = coor(2,2);
            y3   = coor(3,2);
            % Determinant of Jacobian
            detJ = (x1-x3)*(y2-y3)-(y1-y3)*(x2-x3);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Shape functions
        function N = evalShpFn(obj,x,y)
            % Get nodal coordinates
            coor = getNodesCoor(obj);
            x1   = coor(1,1);
            x2   = coor(2,1);
            x3   = coor(3,1);
            y1   = coor(1,2);
            y2   = coor(2,2);
            y3   = coor(3,2);
            % Get area
            A = getA(obj);
            % Eval N
            N1 = (1/(2*A))*(x2*y3-x3*y2+(y2-y3)*x+(x3-x2)*y);
            N2 = (1/(2*A))*(x3*y1-x1*y3+(y3-y1)*x+(x1-x3)*y);
            N3 = (1/(2*A))*(x1*y2-x2*y1+(y1-y2)*x+(x2-x1)*y);
            N = [N1 N2 N3];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Derivative of shape functions w.r.t X
        function dNdX = evaldNdX(obj,x,y)
            % Get nodal coordinates
            coor = getNodesCoor(obj);
            x1   = coor(1,1);
            x2   = coor(2,1);
            x3   = coor(3,1);
            y1   = coor(1,2);
            y2   = coor(2,2);
            y3   = coor(3,2);
            % Get area
            A = getA(obj);
            % Eval dNdX
            % X
            dN1dX = (y2-y3)/(2*A);
            dN2dX = (y3-y1)/(2*A);
            dN3dX = (y1-y2)/(2*A);
            % Y
            dN1dY = (x3-x2)/(2*A);
            dN2dY = (x1-x3)/(2*A);
            dN3dY = (x2-x1)/(2*A);
            %
            dNdX = [ dN1dX dN2dX dN3dX;...
                     dN1dY dN2dY dN3dY ];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Evaluate B matrix
        function B = evalB(obj,xi,eta)
            % Do mapping 
            [x,y] = obj.doMap(xi,eta);
            % Eval dN/dX
            dNdX = evaldNdX(obj,x,y);
            % Assemble B
            B = [ dNdX(1,1) dNdX(1,2) dNdX(1,3);...
                  dNdX(2,1) dNdX(2,2) dNdX(2,3) ];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Evaluate matrix of shape functions
        function N = evalN(obj,xi,eta)
            % Do mapping 
            [x,y] = obj.doMap(xi,eta);
            % Eval shape functions
            Ni = obj.evalShpFn(x,y);
            % Assemble N
            N = [ Ni(1) Ni(2) Ni(3) ];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Evaluate matrix of vector gradient
        function grad = evalGrad(obj,xi,eta)
            % Do mapping 
            [x,y] = obj.doMap(xi,eta);
            % Eval dN/dX
            dNdX = evaldNdX(obj,x,y);
            % Assemble N
            grad = [ dNdX(1,1) dNdX(2,1)...
                     dNdX(1,2) dNdX(2,2)...
                     dNdX(1,3) dNdX(2,3) ];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculate Ke = sum ( B'*B*detJ*w ) |_gp
        function myKe = getKe(obj)
            if obj.isKeComp == false
                % Initialize
                myKe = zeros(3,3);
                % Loop over integration points
                nGP = obj.intRule.nP;
                gP  = obj.intRule.gP;
                wP  = obj.intRule.gW;
                for i = 1:nGP
                    % Parametric constants
                    xi  = gP(i,1);
                    eta = gP(i,2);
                    w   = wP(i,1);
                    % Eval B
                    B = obj.evalB(xi,eta);
                    % Eval J
                    detJ = obj.getDetJ();
                    % Add contribution
                    myKe = myKe + B'*B*detJ*w;
                end
                % Store and set flag to true
                obj.Ke       = myKe;
                obj.isKeComp = true;
            else
                myKe = obj.Ke;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculate div of velocity field using velocity mesh
        function divUe = evalDivUe(obj,Us,tri6)
            % Initialize
            divUe = zeros(3,1);
            % Loop over integration points
            nGP = obj.intRule.nP;
            gP  = obj.intRule.gP;
            wP  = obj.intRule.gW;
            for i = 1:nGP
                % Parametric constants
                xi  = gP(i,1);
                eta = gP(i,2);
                w   = wP(i,1);
                % Eval N using Tri3
                N = evalN(obj,xi,eta);
                % Eval div of vector field using Tri6
                divU = tri6.evalDivU(Us,xi,eta);
                % Eval J
                detJ = obj.getDetJ();
                % Add contribution
                divUe = divUe + N'*divU*detJ*w;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Evaluate gradient of scalar field
        function gradPe = evalGradP(obj,P,xi,eta)
            % Do mapping 
            [x,y] = obj.doMap(xi,eta);
            % Eval dN/dX
            dNdX = evaldNdX(obj,x,y);
            % Assemble grad matrix
            grad = [ dNdX(1,1) dNdX(1,2) dNdX(1,3);...
                     dNdX(2,1) dNdX(2,2) dNdX(2,3) ];
            % Get element local solution coefficients
            Pe   = zeros(3,1);
            nNod = size(obj.nodMap,1);
            % Loop over nodes
            for i = 1:nNod
                iterNod = obj.nodMap(obj.conn(1,i));
                gDoF    = iterNod.gDoF;
                Pe(i)   = P(gDoF);
            end
            % Compute gradient of scalar
            gradPe = grad*Pe;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Evaluate element momentum residual
        function RMe = evalRMe(obj,Ue,Uke,xi,eta)
            % Get constants
            ro    = obj.mat.ro;
            mu    = obj.mat.mu;
            kappa = obj.mat.kappa;
            dt    = obj.mat.dt;
            CL    = obj.mat.CL;
            % Eval shape functions
            Ni = evalN(obj,xi,eta);
            Nx = [ Ni(1) Ni(2) Ni(3) , 0 0 0 , 0 0 0 ];
            Ny = [ 0 0 0 , Ni(1) Ni(2) Ni(3) , 0 0 0 ];
            % Eval derivative of shape functions
            Bi     = evalB(obj,xi,eta);
            dNx_dx = [ Bi(1,1) Bi(1,2) Bi(1,3) , 0 0 0 , 0 0 0 ];
            dNx_dy = [ Bi(2,1) Bi(2,2) Bi(2,3) , 0 0 0 , 0 0 0 ];
            dNy_dx = [ 0 0 0 , Bi(1,1) Bi(1,2) Bi(1,3) , 0 0 0 ];
            dNy_dy = [ 0 0 0 , Bi(2,1) Bi(2,2) Bi(2,3) , 0 0 0 ];
            dNp_dx = [ 0 0 0 , 0 0 0 , Bi(1,1) Bi(1,2) Bi(1,3) ];
            dNp_dy = [ 0 0 0 , 0 0 0 , Bi(2,1) Bi(2,2) Bi(2,3) ];
            % Term 1 - $u,t$
            T1 = (ro/dt)*[Nx ; Ny]*(Uke-Ue);
            % Term 2 - $u \cdot \nabla u$
            Uex   = Nx*Uke;
            Uey   = Ny*Uke;
            dUxdX = dNx_dx*Uke;
            dUydY = dNy_dy*Uke;
            dUxdY = dNx_dy*Uke;
            dUydX = dNy_dx*Uke;
            T2 = ro*[ Uex*dUxdX + Uey*dUxdY ;...
                      Uex*dUydX + Uey*dUydY ];
            % Term 3 - b: not implemented
            T3 = ro*zeros(2,1);
            % Term 4 - $\nabla p$
            T4 = [ dNp_dx ; dNp_dy ]*Uke;
            % Term 5 - $2\mu \Delta u$: Ignored due to linear tri3 elements
            T5 = zeros(2,1);
            % Element momentum residual at integration point
            RMe = T1+T2+T3+T4+T5;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Evaluate element continuity residual
        function RCe = evalRCe(obj,Ue,Uke,xi,eta)
            % Eval derivative of shape functions
            Bi     = evalB(obj,xi,eta);
            dNx_dx = [ Bi(1,1) Bi(1,2) Bi(1,3) , 0 0 0 , 0 0 0 ];
            dNx_dy = [ Bi(2,1) Bi(2,2) Bi(2,3) , 0 0 0 , 0 0 0 ];
            dNy_dx = [ 0 0 0 , Bi(1,1) Bi(1,2) Bi(1,3) , 0 0 0 ];
            dNy_dy = [ 0 0 0 , Bi(2,1) Bi(2,2) Bi(2,3) , 0 0 0 ];
            dNp_dx = [ 0 0 0 , 0 0 0 , Bi(1,1) Bi(1,2) Bi(1,3) ];
            dNp_dy = [ 0 0 0 , 0 0 0 , Bi(2,1) Bi(2,2) Bi(2,3) ];
            % Element continuty residual at integration point
            RCe = dNx_dx*Uke + ...
                  dNy_dy*Uke ;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Evaluate stabilitazion paramters
        function [tauM,tauC] = evalTau(obj,Uke,xi,eta)
            % Get constants
            ro    = obj.mat.ro;
            mu    = obj.mat.mu;
            kappa = obj.mat.kappa;
            dt    = obj.mat.dt;
            CL    = obj.mat.CL;
            % Eval shape functions and velocity components
            Ni = evalN(obj,xi,eta);
            Nx = [ Ni(1) Ni(2) Ni(3) , 0 0 0 , 0 0 0 ];
            Ny = [ 0 0 0 , Ni(1) Ni(2) Ni(3) , 0 0 0 ];
            u  = [ Nx ; Ny ]*Uke;
            % Get element nodal coordinates
            coor = getNodesCoor(obj);
            x1   = coor(1,1);
            x2   = coor(2,1);
            x3   = coor(3,1);
            y1   = coor(1,2);
            y2   = coor(2,2);
            y3   = coor(3,2);
            % Derivates of parametric coors wrt physical coor
            dXidX  =  (y2-y3) / ( (x1-x3)*(y2-y3) - (y1-y3)*(x2-x3) );
            dXidY  = -(x2-x3) / ( (x1-x3)*(y2-y3) - (y1-y3)*(x2-x3) );
            dEtadX = -(y1-y3) / ( (x1-x3)*(y2-y3) - (y1-y3)*(x2-x3) );
            dEtadY =  (x1-x3) / ( (x1-x3)*(y2-y3) - (y1-y3)*(x2-x3) );
            % Get mesh metric tensor
            G      = zeros(2,2);
            G(1,1) = dXidX*dXidX + dEtadX*dEtadX;
            G(1,2) = dXidX*dXidY + dEtadX*dEtadY;
            G(2,1) = dXidY*dXidX + dEtadY*dEtadX;
            G(2,2) = dXidY*dXidY + dEtadY*dEtadY;
            % TauM
            tauM = (1/ro)*...
                   ( 4/dt^2 + u'*(G*u) + CL^2*(mu/ro)^2*(sum(sum(G.*G))) )^(-1/2);
            % Tau
            tau = (1/ro)*...
                  ( u'*(G*u) + CL^2*(mu/ro)^2*(sum(sum(G.*G))) )^(-1/2);
            % TauC
            tauC = 1 / (tau*trace(G));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Evaluate element residual vector
        function resE = evalResE(obj,U,Uk)
            % Initialize
            resE = zeros(9,1);
            % Get constants
            ro    = obj.mat.ro;
            mu    = obj.mat.mu;
            kappa = obj.mat.kappa;
            dt    = obj.mat.dt;
            CL    = obj.mat.CL;
            % Get element velocity and pressure solution coefficients
            Ue  = zeros(9,1);
            Uke = zeros(9,1);
            % Number of nodes
            nNod = size(obj.nodMap,1);
            % Loop over nodes
            for i = 1:nNod
                iterNod = obj.nodMap(obj.conn(1,i));
                gDoF    = iterNod.gDoF;
                % Velocity and pressure
                Ue ([i , i+3 , i+6]) = U(gDoF);
                Uke([i , i+3 , i+6]) = Uk(gDoF);
            end
            % Average solution at n+1/2
            Use = (Ue+Uke)/2;
            % Loop over integration points
            nGP = obj.intRule.nP;
            gP  = obj.intRule.gP;
            wP  = obj.intRule.gW;
            for i = 1:nGP
                % Parametric constants
                xi  = gP(i,1);
                eta = gP(i,2);
                w   = wP(i,1);
                % Eval shape functions
                Ni = evalN(obj,xi,eta);
                Nx = [ Ni(1) Ni(2) Ni(3) , 0 0 0 , 0 0 0 ];
                Ny = [ 0 0 0 , Ni(1) Ni(2) Ni(3) , 0 0 0 ];
                Np = [ 0 0 0 , 0 0 0 , Ni(1) Ni(2) Ni(3) ];
                % Eval derivative of shape functions
                Bi     = obj.evalB(xi,eta);
                dNx_dx = [ Bi(1,1) Bi(1,2) Bi(1,3) , 0 0 0 , 0 0 0 ];
                dNx_dy = [ Bi(2,1) Bi(2,2) Bi(2,3) , 0 0 0 , 0 0 0 ];
                dNy_dx = [ 0 0 0 , Bi(1,1) Bi(1,2) Bi(1,3) , 0 0 0 ];
                dNy_dy = [ 0 0 0 , Bi(2,1) Bi(2,2) Bi(2,3) , 0 0 0 ];
                dNp_dx = [ 0 0 0 , 0 0 0 , Bi(1,1) Bi(1,2) Bi(1,3) ];
                dNp_dy = [ 0 0 0 , 0 0 0 , Bi(2,1) Bi(2,2) Bi(2,3) ];
                % Eval determinant of Jacobian
                detJ = obj.getDetJ();
                % Compute fields
                Uex  = Nx*Ue;   Ukex = Nx*Uke;   Usx = Nx*Use;
                Uey  = Ny*Ue;   Ukey = Ny*Uke;   Usy = Ny*Use;
                                                 Ps  = Np*Use;
                % Copute derivative of fields
                dUsx_dx = dNx_dx*Use;
                dUsx_dy = dNx_dy*Use;
                dUsy_dx = dNy_dx*Use;
                dUsy_dy = dNy_dy*Use;
                % Compute stress tensor components
                sig_xx = -Ps + (2*mu)*(1/2)*(dUsx_dx+dUsx_dx);
                sig_xy =   0 + (2*mu)*(1/2)*(dUsy_dx+dUsx_dy);
                sig_yx =   0 + (2*mu)*(1/2)*(dUsx_dy+dUsy_dx);
                sig_yy = -Ps + (2*mu)*(1/2)*(dUsy_dy+dUsy_dy);
                % Compute element momemtum residual                         JAL: Uke or Use ?
                RMe = obj.evalRMe(Ue,Uke,xi,eta);
                % Compute element continuity residual                       JAL: Uke or Use ?
                RCe = obj.evalRCe(Ue,Uke,xi,eta);
                % Compute element stabilization parameters                  JAL: Uke or Use ?
                [tauM,tauC] = obj.evalTau(Uke,xi,eta);
                % Residual Term 1
                resT1 = Nx' * ro*(Ukex-Uex)/dt + ...
                        Ny' * ro*(Ukey-Uey)/dt;
                % Residual Term 2
                resT2 = Nx' * ro*(Usx*dUsx_dx) + ...
                        Nx' * ro*(Usy*dUsx_dy) + ...
                        Ny' * ro*(Usx*dUsy_dx) + ...
                        Ny' * ro*(Usy*dUsy_dy) ;
                % Residual Term 3
                resT3 = dNx_dx' * sig_xx + ...
                        dNx_dy' * sig_xy + ...
                        dNy_dx' * sig_yx + ...
                        dNy_dy' * sig_yy ;
                % Residual Term 4
                resT4 = Np' * dUsx_dx + ...
                        Np' * dUsy_dy ;
                % Residual Term 5
                gamma = [Usx; Usy] * (tauM*RMe)' + ...
                        (tauM*RMe) * [Usx; Usy]' + ...
                        (tauM*RMe) * (tauM*RMe)' ;
                resT5 = dNx_dx' * (ro*gamma(1,1)) + ...
                        dNx_dy' * (ro*gamma(1,2)) + ...
                        dNy_dx' * (ro*gamma(2,1)) + ...
                        dNy_dy' * (ro*gamma(2,2)) ;
                % Residual Term 6
                resT6 = dNx_dx' * (tauC*RCe) + ...
                        dNy_dy' * (tauC*RCe) ;
                % Residual Term 7
                resT7 = dNp_dx' * (tauM*RMe(1,1)) + ... 
                        dNp_dy' * (tauM*RMe(2,1)) ;
                % Add contribution
                resE = resE + ...
                       ( resT1 + ... 
                         resT2 + ... 
                         resT3 + ...
                         resT4 + ...
                         resT5 + ...
                         resT6 + ...
                         resT7 ) * (detJ*w);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Evaluate element consistent Jacobian matrix
        function KTe = evalKTe(obj,U,Uk)
            % Initialize
            KTe = zeros(9,9);
            % Get constants
            ro    = obj.mat.ro;
            mu    = obj.mat.mu;
            kappa = obj.mat.kappa;
            dt    = obj.mat.dt;
            CL    = obj.mat.CL;
            % Get element velocity and pressure solution coefficients
            Ue  = zeros(9,1);
            Uke = zeros(9,1);
            % Number of nodes
            nNod = size(obj.nodMap,1);
            % Loop over nodes
            for i = 1:nNod
                iterNod = obj.nodMap(obj.conn(1,i));
                gDoF    = iterNod.gDoF;
                % Velocity and pressure
                Ue ([i , i+3 , i+6]) = U(gDoF);
                Uke([i , i+3 , i+6]) = Uk(gDoF);
            end
            % Average solution at n+1/2
            Use = (Ue+Uke)/2;
            % Loop over integration points
            nGP = obj.intRule.nP;
            gP  = obj.intRule.gP;
            wP  = obj.intRule.gW;
            for i = 1:nGP
                % Parametric constants
                xi  = gP(i,1);
                eta = gP(i,2);
                w   = wP(i,1);
                % Eval shape functions
                Ni = evalN(obj,xi,eta);
                Nx = [ Ni(1) Ni(2) Ni(3) , 0 0 0 , 0 0 0 ];
                Ny = [ 0 0 0 , Ni(1) Ni(2) Ni(3) , 0 0 0 ];
                Np = [ 0 0 0 , 0 0 0 , Ni(1) Ni(2) Ni(3) ];
                % Eval derivative of shape functions
                Bi     = obj.evalB(xi,eta);
                dNx_dx = [ Bi(1,1) Bi(1,2) Bi(1,3) , 0 0 0 , 0 0 0 ];
                dNx_dy = [ Bi(2,1) Bi(2,2) Bi(2,3) , 0 0 0 , 0 0 0 ];
                dNy_dx = [ 0 0 0 , Bi(1,1) Bi(1,2) Bi(1,3) , 0 0 0 ];
                dNy_dy = [ 0 0 0 , Bi(2,1) Bi(2,2) Bi(2,3) , 0 0 0 ];
                dNp_dx = [ 0 0 0 , 0 0 0 , Bi(1,1) Bi(1,2) Bi(1,3) ];
                dNp_dy = [ 0 0 0 , 0 0 0 , Bi(2,1) Bi(2,2) Bi(2,3) ];
                % Eval determinant of Jacobian
                detJ = obj.getDetJ();
                % Compute fields
                Uex  = Nx*Ue;   Ukex = Nx*Uke;   Usx = Nx*Use;
                Uey  = Ny*Ue;   Ukey = Ny*Uke;   Usy = Ny*Use;
                                                 Ps  = Np*Use;
                % Copute derivative of fields
                dUsx_dx = dNx_dx*Use;
                dUsx_dy = dNx_dy*Use;
                dUsy_dx = dNy_dx*Use;
                dUsy_dy = dNy_dy*Use;
                % Compute element momemtum residual                         JAL: Uke or Use ?
                RMe = obj.evalRMe(Ue,Uke,xi,eta);
                % Compute element continuity residual                       JAL: Uke or Use ?
                RCe = obj.evalRCe(Ue,Uke,xi,eta);
                % Compute element stabilization parameters                  JAL: Uke or Use ?
                [tauM,tauC] = obj.evalTau(Uke,xi,eta);
                % RM - Term 1 contribution to KTe
                KTe_T1v = ( Nx' * (ro/dt*Nx) ) + ...
                          ( Ny' * (ro/dt*Ny) ) ;
                KTe_T1p = zeros(9,9);
                % RM - Term 2 contribution to KTe
                KTe_T2v = ( Nx' * (ro*Nx*dUsx_dx) ) + ( Nx' * (ro*Usx*dNx_dx) ) + ...
                          ( Nx' * (ro*Ny*dUsx_dy) ) + ( Nx' * (ro*Usy*dNx_dy) ) + ...
                          ( Ny' * (ro*Nx*dUsy_dx) ) + ( Ny' * (ro*Usx*dNy_dx) ) + ...
                          ( Ny' * (ro*Ny*dUsy_dy) ) + ( Ny' * (ro*Usy*dNy_dy) ) ;
                KTe_T2p = zeros(9,9);
                % RM - Term 3 contribution to KTe
                KTe_T3v = ( dNx_dx' * (2*mu*1/2*dNx_dx) ) + ( dNx_dx' * (2*mu*1/2*dNx_dx) ) + ...
                          ( dNx_dy' * (2*mu*1/2*dNx_dy) ) + ( dNx_dy' * (2*mu*1/2*dNy_dx) ) + ...
                          ( dNy_dx' * (2*mu*1/2*dNy_dx) ) + ( dNy_dx' * (2*mu*1/2*dNx_dy) ) + ...
                          ( dNy_dy' * (2*mu*1/2*dNy_dy) ) + ( dNy_dy' * (2*mu*1/2*dNy_dy) ) ;
                KTe_T3p = ( dNx_dx' * -Np ) + ...
                          ( dNy_dy' * -Np ) ;
                % RM - Term 4 contribution to KTe
                KTe_T4v = ( dNx_dx' * (ro*Nx*tauM*RMe(1,1)) ) + ( dNx_dy' * (ro*Nx*tauM*RMe(2,1)) ) + ...
                          ( dNy_dx' * (ro*Ny*tauM*RMe(1,1)) ) + ( dNy_dy' * (ro*Ny*tauM*RMe(2,1)) ) + ...
                          ( dNx_dx' * (ro*tauM*RMe(1,1)*Nx) ) + ( dNx_dy' * (ro*tauM*RMe(1,1)*Ny) ) + ...
                          ( dNy_dx' * (ro*tauM*RMe(2,1)*Nx) ) + ( dNy_dy' * (ro*tauM*RMe(2,1)*Ny) ) ;
                KTe_T4p = zeros(9,9);
                % RM - Term 5 contribution to KTe
                KTe_T5v = ( dNx_dx' * (tauC*(dNx_dx+dNy_dy)) ) + ...
                          ( dNy_dy' * (tauC*(dNx_dx+dNy_dy)) ) ;
                KTe_T5p = zeros(9,9);
                % RC - Term 6 contribution to KTe
                KTe_T6v = ( Np' * dNx_dx ) + ...
                          ( Np' * dNy_dy );
                KTe_T6p = zeros(9,9);
                % RC - Term 7 contribution to KTe
                KTe_T7v = ( dNp_dx' * (tauM*(ro*(Nx/dt+Nx*dUsx_dx+Usx*dNx_dx+Ny*dUsx_dy+Usy*dNx_dy))) ) + ...
                          ( dNp_dy' * (tauM*(ro*(Ny/dt+Nx*dUsy_dx+Usx*dNy_dx+Ny*dUsy_dy+Usy*dNy_dy))) ) ;
                KTe_T7p = ( dNp_dx' * (tauM*dNp_dx) ) + ...
                          ( dNp_dy' * (tauM*dNp_dy) ) ;
                % Add contribution
                KTe = KTe + ...
                      ( KTe_T1v + KTe_T1p + ...
                        KTe_T2v + KTe_T2p + ...
                        KTe_T3v + KTe_T3p + ...
                        KTe_T4v + KTe_T4p + ...
                        KTe_T5v + KTe_T5p + ...
                        KTe_T6v + KTe_T6p + ...
                        KTe_T7v + KTe_T7p ) * (detJ*w);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end

