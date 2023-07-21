classdef VSXImageBuffer < matlab.mixin.Copyable
    properties
% %         datatype (1,:) char = 'double'
%         rowsPerFrame (1,1) double 
%         colsPerFrame (1,1) double
%         sectionsPerFrame (1,1) double
        numFrames (1,1) double = 1
% %         firstFrame (1,1) double
% %         lastFrame (1,1) double
    end   
    
    
    methods
        function obj = VSXImageBuffer(kwargs)
            arguments
                kwargs.?VSXImageBuffer
            end
            for f = string(fieldnames(kwargs))'
                obj.(f) = kwargs.(f);
            end
        end
    end
end
