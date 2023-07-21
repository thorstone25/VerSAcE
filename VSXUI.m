classdef VSXUI < matlab.mixin.Copyable
    properties
        Statement (1,:) char
        Control cell
        Callback cell
        handle function_handle 
    end
end
