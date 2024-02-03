classdef gaussint < handle
    %GAUSSINT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Rule order
        nP;
        % Points and wegihts
        gP;
        gW;
    end
    
    methods
        % Constructor
        function obj = gaussint(p)
            %GAUSSINT Construct an instance of this class
            %   Detailed explanation goes here
            if p == 1
                % Order
                obj.nP = p;
                % Nodes
                obj.gP = [1/3 1/3];
                % Weight
                obj.gW = 1/2;
            elseif p == 3
                % Order
                obj.nP = p;
                % Nodes
                obj.gP = [1/6 2/3;...
                          1/6 1/6;...
                          2/3 1/6];
                % Weight
                obj.gW = [1/6; 1/6; 1/6];
            else
                disp('IntRule not supported. Use [ 1 | 3 ]');
                return
            end
        end
    end
end

