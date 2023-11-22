classdef VSXUI < matlab.mixin.Copyable
    properties
        Statement (:,1) string = string.empty
        Control cell = {}
        Callback cell = {}
        % handle function_handle (set by VSX)
    end
    methods
        function obj = VSXUI(kwargs)
            arguments
                kwargs.?VSXUI
            end
            for f = string(fieldnames(kwargs))'
                obj.(f) = kwargs.(f);
            end
        end
    end
end
