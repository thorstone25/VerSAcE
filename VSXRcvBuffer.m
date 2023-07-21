classdef VSXRcvBuffer < matlab.mixin.Copyable
    properties
        datatype (1,:) char = 'int16'
        rowsPerFrame (1,1) double 
        colsPerFrame (1,1) double = 128
        numFrames (1,1) double = 1
    end
    
    methods
        function obj = VSXRcvBuffer(kwargs)
            arguments, kwargs.?VSXRcvBuffer, end % arguments
            for f = string(fieldnames(kwargs))', obj.(f) = kwargs.(f); end % init

            value = obj.numFrames;
            if ~(isequal(value, 1) || rem(value, 2) == 0)
                error('Property value must be 1 or an even number.');
            end     
        end
    end
end
