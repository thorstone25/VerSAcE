function [vUI, vEvent] = addReconFun(fnm, scan, vRcvBuf, vRes, vSeq, vPData, kwargs)
arguments
    fnm (1,1) string
    scan (1,1) Scan
    vRcvBuf (1,1) VSXRcvBuffer
    vRes (1,1) VSXResource
    vSeq (1,:) VSXSeqControl = VSXSeqControl('command', 'noop', 'argument', 100/0.2);
    vPData {mustBeScalarOrEmpty} = VSXPData.QUPS(scan)
    kwargs.numFrames (1,1) double = 1
    kwargs.multipage (1,1) logical = false
    kwargs.display (1,1) logical = true
    kwargs.name (1,1) string = "Do " + fnm;
    kwargs.UItyp (1,1) string {mustBeMember(kwargs.UItyp,["VsToggleButton", "VsPushButton", "None"])} = "None"
    kwargs.UIpos (1,1) string {mustBeMember(kwargs.UIpos, ["auto","UserB1","UserB2","UserB3","UserB4","UserB5","UserB6","UserB7","UserB8","UserC1","UserC2","UserC3","UserC4","UserC5","UserC6","UserC7","UserC8","UserA1","UserA2"])} = "auto";
end

% Image buffer
vbuf_im    = VSXImageBuffer.fromPData(vPData);

% Display window
vDisplayWindow = VSXDisplayWindow.QUPS(scan, ...
    'Title', char(kwargs.name), ...
    'numFrames', kwargs.numFrames, ...
    'AxesUnits', 'mm', ...
    'Colormap', gray(256) ...
    );

% Process
compute_image_process = VSXProcess('classname', 'External', 'method', fnm);
compute_image_process.Parameters = struct( ...
    'srcbuffer','receive',...
    'srcbufnum', vRcvBuf,...
    'srcframenum',0,... 0 for all frames
    'dstbuffer','image', ...
    'dstbufnum', vbuf_im, ...
    'dstframenum', -1 ... -1 for last frame
    );

display_image_process = VSXProcess('classname', 'Image', 'method', 'imageDisplay');
display_image_process.Parameters = struct( ...
    'imgbufnum', vbuf_im,...   % number of buffer to process.
    'framenum',-1,...   % (-1 => lastFrame)
    'pdatanum', vPData,...    % PData structure to use
    'pgain',1.0,...            % pgain is image processing gain
    'reject',2,...      % reject level
    'persistMethod','simple',...
    'persistLevel',20,...
    'interpMethod','4pt',...
    'grainRemoval','none',...
    'processMethod','none',...
    'averageMethod','none',...
    'compressMethod','power',...
    'compressFactor',40,...
    'mappingMethod','full',...
    'display',1,...      % display image after processing
    'displayWindow', vDisplayWindow ...
    );

% Event
vEvent = VSXEvent(...
    'info', [char(kwargs.name) ' Recon'], ...
    'process', compute_image_process, ...
    'seqControl', vSeq ...
    );

vEvent(2) = VSXEvent(...
    'info', [char(kwargs.name) ' Disp'], ...
    'process', display_image_process, ...
    'seqControl', vSeq ...
    );

% optionally add UI to toggle
switch kwargs.UItyp
    case {"VsToggleButton", "VsPushButton"}
        vUI = VSXUI( ...
            'Control', {char(kwargs.UIpos),'Style',char(kwargs.UItyp),'Label', char(kwargs.name)}, ...
            'Callback', cellstr(["global TOGGLE_"+fnm+"; TOGGLE_"+fnm+" = logical(hObject.Value);", "return;"]) ...
            );
    case "None", vUI = VSXUI.empty;
    otherwise
        error("VerSAcE:addReconFun:buttonTypeNotSupported", "Button not supported.");
end

%% Add to required Resource buffer
if ~isempty(vbuf_im       ), vRes.ImageBuffer(end+1)    = vbuf_im       ; end
if ~isempty(vDisplayWindow), vRes.DisplayWindow(end+1)  = vDisplayWindow; end
