classdef VSXRcvBuffer < matlab.mixin.Copyable
    properties
        datatype (1,1) string {mustBeMember(datatype, ["int16"])} = 'int16'
        rowsPerFrame (1,1) double 
        colsPerFrame (1,1) double = 128
        numFrames (1,1) double {mustBe1OrEven} = 1
    end

    methods
        function obj = VSXRcvBuffer(kwargs)
            arguments, kwargs.?VSXRcvBuffer, end % arguments
            for f = string(fieldnames(kwargs))', obj.(f) = kwargs.(f); end % init
        end
    end
end

% validator(s)
function mustBe1OrEven(value)
if ~(isequal(value, 1) || rem(value, 2) == 0)
    error('NumFrames must be 1 or an even number.');
end
end
