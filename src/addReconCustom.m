function [vEvent, vUI] = addReconCustom(scan, vbuf_rx, vres, vSeq, vPData, kwargs)
arguments
    scan (1,1) Scan
    vbuf_rx (1,1) VSXRcvBuffer
    vres VSXResource {mustBeScalarOrEmpty} = VSXResource.empty
    vSeq (1,:) VSXSeqControl = VSXSeqControl.empty
    vPData (1,1) = VSXPData.QUPS(scan)
    kwargs.numFrames (1,1) double = 1
    kwargs.multipage (1,1) logical = false
    kwargs.display (1,1) logical = true
end

nm = "RFDataImg"; % function name

% Image buffer
vbuf_inter = VSXInterBuffer.fromPData(vPData, 'numFrames', kwargs.numFrames); 
vbuf_im    = VSXImageBuffer.fromPData(vPData);

% Display window
vDisplayWindow = VSXDisplayWindow.QUPS(scan, ...
    'Title', 'QUPS Recon', ...
    'numFrames', kwargs.numFrames, ...
    'AxesUnits', 'mm', ...
    'Colormap', gray(256) ...
    );

% Process
compute_image_process = VSXProcess('classname', 'External', 'method', nm);
compute_image_process.Parameters = struct( ...
    'srcbuffer','receive',...
    'srcbufnum', vbuf_rx,...
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
    'info', 'QUPS Recon', ...
    'process', compute_image_process, ...
    'seqControl', vSeq ...
    );

vEvent(2) = VSXEvent(...
    'info', 'QUPS Disp', ...
    'process', display_image_process, ...
    'seqControl', vSeq ...
    );

% UI
vUI = VSXUI( ...
    'Control', {'auto','Style','VsToggleButton','Label', 'Image RFData', 'Callback', str2func("do"+nm)}, ...
    'Statement', cellstr(["global TOGGLE_"+nm+"; TOGGLE_"+nm+" = false; return;"]), ... init
    'Callback', cellstr(["do"+nm+"(varargin)", "global TOGGLE_"+nm+"RFDataImg; TOGGLE_"+nm+" = logical(UIState); return;"]) ...
    );

%% Add to required Resource buffer
if isempty(vres)
    warning("VerSAcE:addReconCustom:NoResource", "No resource given; allocated buffers may not be properly linked.");
else
    if ~isempty(vbuf_im       ), vres.ImageBuffer(end+1)    = vbuf_im       ; end
    if ~isempty(vbuf_inter    ), vres.InterBuffer(end+1)    = vbuf_inter    ; end
    if ~isempty(vDisplayWindow), vres.DisplayWindow(end+1)  = vDisplayWindow; end
end
