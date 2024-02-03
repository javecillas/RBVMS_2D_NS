classdef bcinfo < handle
    %BCINFO Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nodId;
        bcVal;
        bcFlag;
    end
    
    methods
        function obj = bcinfo(nodId,bcVal,bcFlag)
            obj.nodId  = nodId;
            obj.bcVal  = bcVal;
            obj.bcFlag = bcFlag;
        end
    end
end

