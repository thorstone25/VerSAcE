classdef VSXPData < matlab.mixin.Copyable
    properties
        Coord (1,1) string {mustBeMember(Coord, ["rectangular",...
                                                 "polar",...
                                                 "spherical"])} = "rectangular"
        PDelta (1,:) double
        Size (1,:) double
        Origin (1,:) double
% %         Region struct
    end
    methods
        function obj = VSXPData(kwargs)
            arguments
                kwargs.?VSXPData
            end
            for f = string(fieldnames(kwargs))'
                obj.(f) = kwargs.(f); 
            end
        end
    end
end
