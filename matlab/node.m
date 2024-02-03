classdef node < handle
    %NODE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        iD;
        coor;
        gDoF = [-1 -1];
        % Check if node has velocity DoF
        isDofVSet = false;
        % Check if node has pressure DoF
        isDofPSet = false;
        % Check if node has velocity AND pressure DoF
        isDofVPSet = false;
        % Check if node has Dirichlet BC
        hasDirBC = false;
        bcVal;
        bcFlag;
        % Flag to check if is pressure node
        isVNod = true;
        isPNod = false;
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        function obj = node(iD,coor)
            % ID
            obj.iD   = iD;
            % Coordinate
            obj.coor = coor;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function nextDoF = setDoFTri3V(obj,dof)
            if obj.isDofVSet  == false
                obj.gDoF      = [dof dof+1];
                obj.isDofVSet = true;
                nextDoF       = dof+2;
            else
                nextDoF = dof;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function nextDoF = setDoFTri3P(obj,dof)
            if obj.isDofPSet  == false && obj.isDofVSet  == true
                obj.gDoF      = [obj.gDoF dof];
                obj.isDofPSet = true;
                nextDoF       = dof+1;
            else
                nextDoF = dof;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function nextDoF = setDoFTri3VP(obj,dof)
            if obj.isDofVPSet  == false
                obj.gDoF       = [dof dof+1 dof+2];
                obj.isDofVPSet = true;
                nextDoF        = dof+3;
            else
                nextDoF = dof;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function nextDof = setDoFTri6(obj,dof)
            if obj.isDofSet == false
                obj.gDoF     = [dof dof+1];
                obj.isDofSet = true;
                nextDof      = dof+2;
            else
                nextDof = dof;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Update coordinate
        function updNodCoor(obj,xNew,yNew)
            obj.coor = [xNew yNew];
        end
    end
end

