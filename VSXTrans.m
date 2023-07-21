classdef VSXTrans < matlab.mixin.Copyable
    properties
        name (1,:) char
        units (1,:) char = 'mm' %{mustBeMember(units, ["wavelengths", "mm"])} = "mm"
        frequency (1,1) double
        numelements (1,1) double 
        connType (1,1) double
    end
    methods
        function obj = VSXTrans(kwargs)
            arguments
                kwargs.?VSXTrans
            end
            for f = string(fieldnames(kwargs))'
                obj.(f) = kwargs.(f); 
            end
        end
    end
end