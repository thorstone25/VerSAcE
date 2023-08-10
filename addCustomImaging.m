function [vUI, vEvent] = addCustomImaging(scan, vResource, vbuf_rx, vPData, vSeq, kwargs)
arguments
    scan Scan
    vResource VSXResource
    vbuf_rx (1,1) VSXRcvBuffer
    vPData {mustBeScalarOrEmpty} = VSXPData.empty
    vSeq (1,:) VSXSeqControl = VSXSeqControl.empty
    kwargs.numFrames (1,1) double = 1
    kwargs.multipage (1,1) logical = false
    kwargs.display (1,1) logical = true
end

% PData
if isempty(vPData), vPData = VSXPData.QUPS(scan); end

% Image buffer
vbuf_inter = VSXInterBuffer.fromPData(vPData, 'numFrames', kwargs.numFrames);
vbuf_im = VSXImageBuffer.fromPData(vPData);

% Display window
vDisplayWindow = VSXDisplayWindow.QUPS(scan, ...
    'Title', 'Custom Imaging', ...
    'numFrames', kwargs.numFrames, ...
    'AxesUnits', 'mm', ...
    'Colormap', gray(256) ...
    );

% Process
% nm = "imagingProc"; % function name
nm = "RFDataCImage";
compute_image_process = VSXProcess('classname', 'External', 'method', nm);
compute_image_process.Parameters = {
    'srcbuffer','receive',...
    'srcbufnum', vbuf_rx,...
    'srcframenum',0,... 0 for all frames
    'dstbuffer','image', ...
    'dstbufnum', vbuf_im, ...
    'dstframenum', -1,... -1 for last frame
    };

display_image_process = VSXProcess('classname', 'Image', 'method', 'imageDisplay');
display_image_process.Parameters = { ...
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
    };

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
    'Control', {'UserB4','Style','VsToggleButton','Label', 'RealTimeImg', 'Callback', str2func("do"+nm)}, ...
    'Statement', cellstr(["global TOGGLE_"+nm+"; TOGGLE_"+nm+" = false; return;"]), ... init
    'Callback', cellstr(["do"+nm+"(varargin)", "global TOGGLE_"+nm+"; TOGGLE_"+nm+" = logical(UIState); return;"]) ...
    );

%% Add to required Resource buffer
if ~isempty(vbuf_im       ), vResource.ImageBuffer(end+1)    = vbuf_im       ; end
if ~isempty(vbuf_inter    ), vResource.InterBuffer(end+1)    = vbuf_inter    ; end
if ~isempty(vDisplayWindow), vResource.DisplayWindow(end+1)  = vDisplayWindow; end
end