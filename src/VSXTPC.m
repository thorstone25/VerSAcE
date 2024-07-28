classdef VSXTPC < matlab.mixin.Copyable
    properties
        name string {mustBeScalarOrEmpty} % (optional) name identifier for profile
        maxHighVoltage double {mustBeScalarOrEmpty} % (optional) max high voltage for this profile
        highVoltageLimit double {mustBeScalarOrEmpty} % (optional) high voltage limit based on use model
        xmitDuration double {mustBeScalarOrEmpty} % (optional) longest transmit duration (usec)
        hv double{mustBeScalarOrEmpty} % (optional) initial High Voltage value at startup
        ishpt (1,1) logical = false % whether high power transmit
    end
    methods
        function obj = VSXTPC(kwargs)
            arguments
                kwargs.?VSXTPC
            end
            for f = string(fieldnames(kwargs))'
                obj.(f) = kwargs.(f); 
            end
        end
    end
end
