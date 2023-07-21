classdef VSXResource < matlab.mixin.Copyable
    properties
        Parameters (1,1) VSXParameters
        RcvBuffer (1,:) VSXRcvBuffer = VSXRcvBuffer.empty
        InterBuffer (1,:) VSXInterBuffer = VSXInterBuffer.empty
        ImageBuffer (1,:) VSXImageBuffer = VSXImageBuffer.empty
        DisplayWindow (1,:) VSXDisplayWindow = VSXDisplayWindow.empty
    end
    methods
        function obj = VSXResource(kwargs)
            arguments
                kwargs.?VSXResource
            end
            for f = string(fieldnames(kwargs))'
                obj.(f) = kwargs.(f); 
            end
        end
    end
end

