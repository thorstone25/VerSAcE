classdef VSXInterBuffer < matlab.mixin.Copyable
    properties
% %         datatype (1,:) char = 'int16'
%         rowsPerFrame (1,1) double %optional if PData defined
%         colsPerFrame (1,1) double
%         sectionsPerFrame (1,1) double
%         pagesPerFrame (1,1) double 
        numFrames (1,1) double = 1
% %         firstFrame (1,1) double
% %         lastFrame (1,1) double
    end
    
    methods
        function obj = VSXInterBuffer(kwargs)
            arguments, kwargs.?VSXInterBuffer, end
            for f = string(fieldnames(kwargs))'
                obj.(f) = kwargs.(f);
            end
        end
    end
end
