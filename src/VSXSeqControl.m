classdef VSXSeqControl < matlab.mixin.Copyable
    properties
        command (1,:) char
        argument (1,1) {mustBeA(argument, ["double", "VSXEvent"])}
        % condition (1,:) char
    end
    methods
        function obj = VSXSeqControl(kwargs)
            arguments
                kwargs.?VSXSeqControl
            end
            for f = string(fieldnames(kwargs))'
                obj.(f) = kwargs.(f); 
            end
        end
    end
end
