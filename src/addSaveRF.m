function [vUI, vEvent] = addSaveRF(vbuf_rx, vSeq, kwargs)
arguments
    vbuf_rx (1,1) VSXRcvBuffer
    vSeq (1,:) VSXSeqControl = VSXSeqControl('command', 'noop', 'argument', 100/0.2);
    kwargs.frames (1,1) string {mustBeMember(kwargs.frames, ["last","all"])} = "last"
end
nm = "RFDataStore";
% select frames
switch kwargs.frames, case "all", f=0; case "last", f=-1; end

%% Process: saving data
save_rf_data = VSXProcess('classname', 'External', 'method', nm, 'Parameters',struct( ...
    'srcbuffer','receive',...
    'srcbufnum', vbuf_rx,... 
    'srcframenum',f,...
    'dstbuffer','none'));

vUI = VSXUI( ...
'Control', {'auto','Style','VsPushButton','Label', 'SAVE RFData'}, ...
'Callback', cellstr(["global TOGGLE_"+nm+"; TOGGLE_"+nm+" = true; return;"]) ...
);

% save all RF Data
vEvent = VSXEvent('info', 'Save RF Data', 'process', save_rf_data, 'seqControl', vSeq);