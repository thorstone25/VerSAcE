classdef VSXSeqControl < matlab.mixin.Copyable
    properties
        command (1,1) string {mustBeMember(command, [ ...
            "call", "cBranch", ... "DMA", ... DMA for internal usage
            "jump", "loopCnt", "loopTst", ...
            "markTransferProcessed", "multiSysSync", "noop", "pause", "returnToMatlab", ...
            "rtn", "setTPCProfile","setRcvProfile","stop","sync", ...
            "timeToNextAcq","timeToNextEB","transferToHost","triggerIn","triggerOut",...
            "waitForTransferComplete" ...
            ])} = "noop"
        argument {mustBeScalarOrEmpty, mustBeA(argument, ["double", "VSXEvent", "VSXTPC", "VSXRcvProfile"])} = [] % TODO: VSXRcvProfile
        condition  string {mustBeScalarOrEmpty} = string.empty
        %{
        , mustBeMember(condition, [ ...
            "exitAfterJump", ... jump
            "Trigger_1_Falling","Trigger_2_Falling","Trigger_1_2_Falling", ... triggerIn
            "Trigger_1_Rising","Trigger_2_Rising","Trigger_1_2_Rising", ... triggerIn
            "syncNone", "syncADC_CLK","syncSYNC_CLK", ... triggerOut (sync-none is def)
            ])} = string.empty
        %}
    end
    methods
        function obj = VSXSeqControl(kwargs)
            arguments
                kwargs.?VSXSeqControl
            end
            for f = string(fieldnames(kwargs))'
                obj.(f) = kwargs.(f); 
            end
        end
    end
end
