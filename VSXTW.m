classdef VSXTW < matlab.mixin.Copyable
    properties
        type (1,:) char = 'parametric'
%                             {mustBeMember(type,["parametric",...
%                                               "envelope",...
%                                               "pulseCode",...
%                                               "states",...
%                                               "function",...
%                                               "sampled"])} = "parametric"
        Parameters double %had (:,4)
    end
    methods
        function obj = VSXTW(kwargs)
            arguments
                kwargs.?VSXTW
            end
            for f = string(fieldnames(kwargs))'
                obj.(f) = kwargs.(f); 
            end
        end
    end
end