classdef VSXImageBuffer < matlab.mixin.Copyable
    properties
        datatype (1,1) string {mustBeMember(datatype, "double")} = "double"
        rowsPerFrame (1,1) double
        colsPerFrame (1,1) double
        sectionsPerFrame (1,1) double
        numFrames (1,1) double = 1
        % firstFrame (1,1) double % read-only
        % lastFrame (1,1) double % read-only
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

    methods(Static)
        function obj = fromPData(vPData, kwargs)
            arguments
                vPData (1,1) VSXPData
                kwargs.datatype (1,1) string {mustBeMember(kwargs.datatype, "double")}
                kwargs.numFrames (1,1)
            end
            args = namedargs2cell(kwargs); % forward arguments
            obj = VSXImageBuffer( ...
                "rowsPerFrame"      , vPData.Size(1), ...
                "colsPerFrame"      , vPData.Size(2), ...
                "sectionsPerFrame"  , vPData.Size(3), ...
                args{:} ...
                );

        end
    end
end
