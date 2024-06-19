function [vUI, vEvent] = addDataFun(fnm, vRcvBuf, vSeq, kwargs)
arguments
    fnm (1,1) string % = "RFDataStore"
    vRcvBuf (1,1) VSXRcvBuffer
    vSeq (1,:) VSXSeqControl = [VSXSeqControl('command','returnToMatlab'), VSXSeqControl('command', 'noop', 'argument', 100/0.2)];
    kwargs.name (1,1) string = "Do "+fnm;
    kwargs.UItyp (1,1) string {mustBeMember(kwargs.UItyp,["VsToggleButton", "VsPushButton", "None"])} = "None"
    kwargs.UIpos (1,1) string {mustBeMember(kwargs.UIpos, ["auto","UserB1","UserB2","UserB3","UserB4","UserB5","UserB6","UserB7","UserB8","UserC1","UserC2","UserC3","UserC4","UserC5","UserC6","UserC7","UserC8","UserA1","UserA2"])} = "auto";
end

% save all RF Data
vEvent = VSXEvent('info', kwargs.name, 'seqControl', vSeq, 'process', ... 
    VSXProcess('classname', 'External', 'method', fnm, 'Parameters',struct( ...
    'srcbuffer','receive',...
    'srcbufnum', vRcvBuf,... 
    'srcframenum', 0,...
    'dstbuffer','none' ...
)))); % Process: saving data

% optionally add UI to toggle
switch kwargs.UItyp
    case {"VsToggleButton", "VsPushButton"}
        vUI = VSXUI( ...
            'Control', {char(kwargs.UIpos),'Style',char(kwargs.UItyp),'Label', char(kwargs.name)}, ...
            'Callback', cellstr(["global TOGGLE_"+fnm+"; TOGGLE_"+fnm+" = logical(hObject.Value);", "return;"]) ...
            );
    case "None", vUI = VSXUI.empty;
    otherwise
        error("VSXOOD:addDataFun:buttonTypeNotSupported", "Button not supported.");
end

