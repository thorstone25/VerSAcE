classdef VSXEvent < matlab.mixin.Copyable
    properties
        info (1,1) string = 'Event'
        tx VSXTX {mustBeScalarOrEmpty} = VSXTX.empty
        rcv VSXReceive {mustBeScalarOrEmpty} = VSXReceive.empty
        recon (1,:) VSXRecon = VSXRecon.empty
        process VSXProcess  {mustBeScalarOrEmpty} = VSXProcess.empty
        seqControl (1,:) VSXSeqControl = VSXSeqControl.empty
    end
    methods
        function obj = VSXEvent(kwargs)
            arguments
                kwargs.?VSXEvent
            end
            for f = string(fieldnames(kwargs))'
                obj.(f) = kwargs.(f); 
            end
        end
    end
end
        