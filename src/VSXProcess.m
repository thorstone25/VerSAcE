classdef VSXProcess < matlab.mixin.Copyable
    properties
        classname (1,1) string {mustBeMember(classname, ["Image", "Doppler", "External"])} = "External"
        method (1,1) string
        Parameters (1,:) cell
    end
    
    methods
        function obj = VSXProcess(kwargs)
            arguments, kwargs.?VSXProcess, end % arguments
            for f = string(fieldnames(kwargs))', obj.(f) = kwargs.(f); end % init
        end
    end
end