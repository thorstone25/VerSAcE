classdef VSXReceive < matlab.mixin.Copyable
    properties
        Apod (1,:) double
%         aperture (1,1) double
        startDepth (1,1) double
        endDepth (1,1) double
        TGC VSXTGC {mustBeScalarOrEmpty} = VSXTGC.empty
        bufnum (1,:) VSXRcvBuffer
        framenum (1,1) double = 1
        acqNum (1,1) double
        sampleMode (1,1) string {mustBeMember(sampleMode,["NS200BW","NS200BWI","BS100BW","BS67BW","BS50BW","custom"])} = "NS200BW"
        mode (1,1) double = 0
        callMediaFunc (1,1) logical = false
        decimSampleRate double {mustBeScalarOrEmpty} = []
    end
    methods
        function obj = VSXReceive(kwargs)
            arguments
                kwargs.?VSXReceive
            end
            for f = string(fieldnames(kwargs))'
                obj.(f) = kwargs.(f); 
            end
        end
    end
end