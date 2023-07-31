function [vBlock, chd, vTrans] = QUPS2VSX(us, xdc, vResource, kwargs)
% QUPS2VSX - Verasonics structure converter
%
% [vBlock, vPData, vTrans] = QUPS2VSX(us) converts the UltrasoundSystem us
% into a VSXBlock vBlock, with a default VSXPData vPData and a Trans
% structure vTrans. These can be used with VSXBlock.link to generate
% Verasonics compatible configuration structure.
%
% [...] = QUPS2VSX(us, xdc) where xdc is a string uses the named transducer
% xdc and the `computeTrans` utility to generate the transducer struct.
% If xdc is a struct, it is assumed to be the verasonics Trans struct and
% the `computeTrans` utility is used to complete the definition.
% If xdc is a Transducer, a custom transducer definition is generated. The 
% default is us.xdc.
%
% [...] = QUPS2VSX(us, xdc, vResource) uses the VSXResource vResource
% instead of creating a new VSXResource. You must specify this if
arguments
    us (1,1) UltrasoundSystem
    xdc (1,1) {mustBeA(xdc, ["string", "struct", "Transducer"])} = us.xdc % transducer name
    vResource (1,1) VSXResource = VSXResource()
    kwargs.units (1,1) string {mustBeMember(kwargs.units, ["mm", "wavelengths"])} = "mm"
    kwargs.vTW (1,1) VSXTW = VSXTW('type', 'parametric', 'Parameters', [us.xdc.fc/1e6, 0.67, 1, 1]); %-
    kwargs.vTGC VSXTGC {mustBeScalarOrEmpty} = VSXTGC('CntrlPts', [0,297,424,515,627,764,871,1000],...
        'rangeMax', hypot(us.scan.zb(2), us.scan.xb(2)) ./ us.lambda); %-
    kwargs.frames (1,1) {mustBeInteger, mustBePositive} = 1;
    kwargs.sample_mode (1,1) string {mustBeMember(kwargs.sample_mode, ["NS200BW", "NS200BWI", "BS100BW",  "BS67BW",  "BS50BW",  "custom"])} = "NS200BW"
    kwargs.custom_fs (1,1) double
    kwargs.recon_VSX (1,1) logical = true
    kwargs.recon_custom (1,1) logical = true
    kwargs.recon_custom_delays (1,1) logical = false
    kwargs.saver_custom (1,1) logical = true
end

% squash obj to struct warning
warning_state = warning('off', 'MATLAB:structOnObject');

% init
vUI     = reshape(VSXUI.empty   , [1 0]);
% vPData  = reshape(VSXPData.empty, [1 0]);

%% Trans
% set sound speed
c0 = us.seq.c0;
vResource.Parameters.speedOfSound = c0;

if isa(xdc, 'string') % interpret as name
    vTrans = struct('name', char(xdc), 'units', char(kwargs.units));
elseif isa(xdc, 'struct') % interpret as the Trans struct
    vTrans = xdc;
elseif isa(xdc, "Transducer") % make custom
    vTrans = xdc.QUPS2VSX();
else
    error("Unrecognized input for xdc.")
end

% have VSX compute remaining properties
vTrans = computeTrans(vTrans);

% convert to QUPS
xdc = Transducer.Verasonics(vTrans);

% set global params
lambda = c0 / (1e6*vTrans.frequency); % wavelengths

% get the scan region in units of wavelengths so we know the ROI
assert(isa(us.scan, 'ScanCartesian')); % should be ScanCartesian
scan = scale(us.scan, 'dist', 1./lambda);

% get temporal sampling region
dnear = floor(2 * scan.zb(1)); % nearest distance (2-way)
dfar  =  ceil(2 * hypot(range(scan.xb), scan.zb(end))); % furthest distance (2-way)

%% TW
% TODO: include TPC
% TODO: handle return of Trans
vTW = computeTWWaveform(kwargs.vTW, vTrans, vResource); % fill out the waveform

%% Allocate buffers and Set Parameters
fs_available = 250 ./ (100:-1:4); % all supported sampling frequencies

% decimation frequency as per sampleMode
if kwargs.sample_mode == "NS200BW"
    fs_decim = fs_available(find(fs_available >= 4 * vTrans.frequency, 1)); 
elseif kwargs.sample_mode == "BS100BW"
    fs_decim = fs_available(find(fs_available >= 2 * vTrans.frequency, 1));
elseif kwargs.sample_mode == "BS67BW"
    fs_decim = fs_available(find(fs_available >= (2*0.67) * vTrans.frequency, 1));
elseif kwargs.sample_mode == "BS50BW"
    fs_decim = fs_available(find(fs_available >= 1 * vTrans.frequency, 1));
elseif kwargs.sample_mode == "custom"
    fs_decim = fs_available(find(fs_available >= kwargs.custom_fs, 1));
end

% get the output data buffer length
spw = fs_decim / vTrans.frequency; % samples per wave
bufLen = 128 * ceil(2 * (dfar - dnear - 1/spw) * spw / 128); % only modulus 128 sample buffer length supported
T = bufLen; % alias

% make new rx buffer
vbuf_rx    = VSXRcvBuffer(  'numFrames', kwargs.frames, ...
    'rowsPerFrame', T * us.seq.numPulse,...  %-
    'colsPerFrame', vResource.Parameters.numRcvChannels...
    );

% add buffers to the resources
vResource.RcvBuffer(end+1)   = vbuf_rx;

%% TX
vTX = copy(repmat(VSXTX('waveform', vTW), [1, us.seq.numPulse]));

% get delay and apodization matrices
delay = - us.seq.delays(xdc) * xdc.fc;
apod  = us.seq.apodization(xdc);
t0    = min(delay, [], 1); % min over elements
delay = delay - t0;

% - Set event specific TX attributes.
for i = 1:us.seq.numPulse
    % set beamforming geometry
    switch us.seq.type
        case "VS"
            vTX(i).FocalPt = us.seq.focus(:,i)/lambda;
            % vTX(i).Origin(1) = us.seq.focus(1,i)/lambda;
            % vTX(i).focus     = us.seq.focus(3,i)/lambda;
        case "PW"
            vTX(i).Steer = deg2rad([us.seq.angles(i), 0]);
        case "FSA"
            %do nothing
    end

    % set delays and apodization
    vTX(i).Apod  =  apod(:,i)';
    vTX(i).Delay = delay(:,i)';

    % TODO: use computeTXDelays instead
    % vTX(i).Delay = computeTXDelays(struct(vTX(i))); % requires base
end

% TODO: hack compute TXPD
for i = 1:us.seq.numPulse
    %         vTX(i).TXPD = computeTXPD(struct(vTX(i)), struct(vPData));
end

%% TGC
kwargs.vTGC.Waveform = computeTGCWaveform(kwargs.vTGC, 1e6*vTrans.frequency);

%% Rcv
% default
vRcv = VSXReceive('startDepth', dnear, 'endDepth', dfar, 'bufnum', vbuf_rx, 'sampleMode', kwargs.sample_mode);
if kwargs.sample_mode == "custom", vRcv.decimSampleRate = fs_decim; end
vRcv.Apod = ones([1, vResource.Parameters.numRcvChannels]);
vRcv.TGC = kwargs.vTGC;

% replicate
vRcv = copy(repmat(vRcv,[us.seq.numPulse, kwargs.frames]));

% - Set event specific Receive attributes.
for f = 1:kwargs.frames
    % move points before (or after?) first receive of the frame
    vRcv(1,f).callMediaFunc = true;
    for i = 1:us.seq.numPulse
        vRcv(i,f).framenum = f;
        vRcv(i,f).acqNum   = i;
    end
end

%% SeqControl
t_puls = round(T / fs_decim) + 50; % pulse wait time in usec
t_frm = t_puls*us.seq.numPulse*1.2; % frame wait time in usec
wait_for_tx_pulse        = VSXSeqControl('command', 'timeToNextAcq', 'argument', t_puls);
wait_for_pulse_sequence  = VSXSeqControl('command', 'timeToNextAcq', 'argument', t_frm ); % max TTNA is 4190000
transfer_to_host         = VSXSeqControl('command', 'transferToHost');
no_operation             = VSXSeqControl('command', 'noop', 'argument', 100/0.2); % 'condition', 'Hw&Sw');

%% Event loop
% loop through all events and frames
% ---------- Events ------------- %
vEvent = copy(repmat(VSXEvent('seqControl', wait_for_tx_pulse), [us.seq.numPulse, kwargs.frames]));

for f = 1:kwargs.frames
    for i = 1:us.seq.numPulse % each transmit
        vEvent(i,f) = VSXEvent('info',"Tx "+i,'tx',vTX(i),'rcv',vRcv(i,f));
    end

    % transfer data to host using the last event of the frame
    vEvent(i,f).seqControl = [wait_for_pulse_sequence, transfer_to_host]; % modify last acquisition vEvent's seqControl
end
    
%% Add Events per frame 
% make recon image and return to MATLAB
if kwargs.recon_VSX
    return_to_matlab = VSXSeqControl('command', 'returnToMatlab');
    [recon_event] = addVSXRecon(scan, vResource, vTX, vRcv, return_to_matlab, 'multipage', false, 'numFrames', kwargs.frames);
    recon_event.recon.newFrameTimeout = 50 * 1e-3*wait_for_pulse_sequence.argument; % wait time in msec ~ heuristic is approx 50 x acquisition time
    vEvent(end+1,:) = recon_event; % assign to end of TX-loop
    vEvent(end,:) = copy(vEvent(end,:)); % copy to make unique Events for each frame
else
    recon_event = VSXEvent.empty;
end

%% Add Events at the end of the loop
% vectorize
vev_post = reshape(VSXEvent.empty, [1,0]);

% custom delays reconstruction process and ui
if kwargs.recon_custom_delays
    [vUI(end+1), vev_post(end+1)] = addCustomDelayRecon(vbuf_rx, no_operation);
end

% custom beamformer reconstruction process and ui
if kwargs.recon_custom
    if isempty(recon_event), vPData = VSXPData.empty; 
    else, vPData = recon_event.recon.pdatanum; % reuse PData
    end
    [vUI(end+1), vev_post(end+1:2)] = addCustomRecon(scan, vResource, vbuf_rx, vPData, no_operation);
end

% custom saving process and ui
if kwargs.saver_custom
    [vUI(end+1), vev_post(end+1)] = addCustomSaver(vbuf_rx, no_operation);
end

% return to start of block after all Events
% TODO: add option for what to do after the end of all events
% vEvent(end+1) = VSXEvent('info', 'Jump back', 'seqControl', ...
%     VSXSeqControl('command', 'jump', 'argument', vEvent(1)) ...
%     );

% ---------- Events ------------- %

%% Block
vBlock = VSXBlock('capture', vEvent, 'post', vev_post, 'next', vEvent(1), 'vUI', vUI);

%% Create a template ChannelData object
% TODO: call computeTWWaveform to get TW.peak correction
t0l = 2 * vTrans.lensCorrection; % lens correction in wavelengths
t0p = - vTW.peak; % peak correction in wavelengths
x = zeros([T us.seq.numPulse vTrans.numelements kwargs.frames, 0], 'single');
chd = ChannelData('data', x, 'fs', 1e6*fs_decim, 't0', (t0 + t0l + t0p)./xdc.fc, 'order', 'TMNF');

%% added External Functions/Callback
% restore warning state
warning(warning_state);

% done!
return; 

function [vEvent, vPData] = addVSXRecon(scan, vResource, vTX, vRcv, vSeq, kwargs)
arguments
    scan Scan
    vResource VSXResource
    vTX VSXTX
    vRcv VSXReceive
    vSeq (1,:) VSXSeqControl = VSXSeqControl.empty
    kwargs.numFrames (1,1) double = 1
    kwargs.multipage (1,1) logical = false
    kwargs.display (1,1) logical = true
end
%% ReconInfo
% We need 1 ReconInfo structures for each transmit
vReconInfo = copy(repmat(VSXReconInfo('mode', 'accumIQ'), size(vTX)));  % default is to accumulate IQ data.

% - Set specific ReconInfo attributes.
for j = 1:numel(vReconInfo) % For each reconinfo object
    vReconInfo(j).txnum  = vTX( j); % tx
    vReconInfo(j).rcvnum = vRcv(j); % rx
    vReconInfo(j).regionnum = j; % region index TODO: VSXRegion
    if kwargs.multipage
        vReconInfo(j).pagenum = j; % page number
    end
end
vReconInfo( 1 ).mode = 'replaceIQ'; % on first tx, replace IQ data
vReconInfo(end).mode = 'accumIQ_replaceIntensity'; % on last tx, accum and "detect"

if isscalar(vReconInfo), vReconInfo.mode = 'replaceIntensity'; end % 1 tx

%% PData
% create the pixel data
vPData = VSXPData.QUPS(scan);

% create a default region for each transmit
% TODO: VSXRegion
vPData.Region = repmat(struct('Shape', struct('Name','Custom'), 'PixelsLA', int32(0:scan.nPix-1), 'numPixels', scan.nPix), size(vReconInfo));

% fill out the regions
% TODO: VSXRegion
vPData.Region = computeRegions(struct(vPData));

%% Make Buffers
% if inter(mediate) buffers needed, create one
if any(contains([vReconInfo.mode], "IQ"))
    % get required number of pages
    if any(contains([vReconInfo.mode], "pages"))
        P = vResource.Parameters.numRcvChannels; % for page acquisitions, each page is a receive channel
    else 
        P = max([1, vReconInfo.pagenum]); % use max page number in reconinfo
    end
    
    % make a buffer
    vbuf_inter = VSXInterBuffer.fromPData(vPData, 'numFrames', kwargs.numFrames, 'pagesPerFrame', P); 
else
    vbuf_inter = VSXInterBuffer.empty; % no buffer
end

% if image buffers needed, create one
if any(contains([vReconInfo.mode], "Intensity"))
    vbuf_im = VSXImageBuffer.fromPData(vPData, 'numFrames', kwargs.numFrames);
else
    vbuf_im = VSXImageBuffer.empty; % no buffer
end

%% Recon
vRecon = VSXRecon('pdatanum', vPData, 'IntBufDest', vbuf_inter, 'ImgBufDest', vbuf_im, "RINums", vReconInfo);

%% Display window
if kwargs.display
    vDisplayWindow = VSXDisplayWindow.QUPS(scan, ...
        'Title', 'VSX Beamformer', ...
        'numFrames', kwargs.numFrames, ...
        'AxesUnits', 'mm', ...
        'Colormap', gray(256) ...
        );

    %% Process
    display_image_process = VSXProcess('classname', 'Image', 'method', 'imageDisplay');
    display_image_process.Parameters = {
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
        'displayWindow', vDisplayWindow, ...
        };

    %% Event
    vEvent = VSXEvent(...
        'info', 'VSX Recon', ...
        'recon', vRecon, ...
        'process', display_image_process, ...
        'seqControl', vSeq ...
        );
else
    vEvent = VSXEvent.empty;
end

%% Add to required Resource buffer
if ~isempty(vbuf_inter    ), vResource.InterBuffer(end+1)    = vbuf_inter    ; end
if ~isempty(vbuf_im       ), vResource.ImageBuffer(end+1)    = vbuf_im       ; end
if ~isempty(vDisplayWindow), vResource.DisplayWindow(end+1)  = vDisplayWindow; end

function [vUI, vEvent] = addCustomRecon(scan, vResource, vbuf_rx, vPData, vSeq, kwargs)
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
vbuf_im = VSXImageBuffer.fromPData(vPData);

% Display window
vDisplayWindow = VSXDisplayWindow.QUPS(scan, ...
    'Title', 'QUPS Recon', ...
    'numFrames', kwargs.numFrames, ...
    'AxesUnits', 'mm', ...
    'Colormap', gray(256) ...
    );

% Process
nm = "RFDataImg"; % function name
compute_image_process = VSXProcess('classname', 'External', 'method', nm);
compute_image_process.Parameters = {
    'srcbuffer','receive',...
    'srcbufnum', vbuf_rx,...
    'srcframenum',0,...
    'dstbuffer','image', ...
    'dstbufnum', vbuf_im, ...
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
    'process', compute_image_process ...
    );

vEvent(2) = VSXEvent(...
    'info', 'QUPS Disp', ...
    'process', display_image_process, ...
    'seqControl', vSeq ...
    );

% UI
vUI = VSXUI( ...
    'Control', {'UserB3','Style','VsToggleButton','Label', 'Image RFData', 'Callback', str2func("do"+nm)}, ...
    'Statement', cellstr(["global TOGGLE_"+nm+"; TOGGLE_"+nm+" = false; return;"]), ... init
    'Callback', cellstr(["do"+nm+"(varargin)", "global TOGGLE_"+nm+"RFDataImg; TOGGLE_"+nm+" = logical(UIState); return;"]) ...
    );

%% Add to required Resource buffer
if ~isempty(vbuf_im       ), vResource.ImageBuffer(end+1)    = vbuf_im       ; end
if ~isempty(vDisplayWindow), vResource.DisplayWindow(end+1)  = vDisplayWindow; end

function [vUI, vEvent] = addCustomDelayRecon(vbuf_rx, vSeq)
arguments
    vbuf_rx (1,1) VSXRcvBuffer
    vSeq (1,:) VSXSeqControl = VSXSeqControl.empty
end
%% Add custom data processing
proc_rf_data = VSXProcess('classname', 'External', 'method', 'RFDataProc');
proc_rf_data.Parameters = {                                                     
    'srcbuffer','receive',...                                                   
    'srcbufnum', vbuf_rx,...                                                    
    'srcframenum',0,...                                                         
    'dstbuffer','none'};                                                        

vUI = VSXUI( ...
'Control', {'UserB2','Style','VsToggleButton','Label', 'Process RFData', 'Callback', @doRFDataProc}, ...
'Statement', cellstr(["global TOGGLE_RFDataProc; TOGGLE_RFDataProc = false; return;"]), ... init
'Callback', cellstr(["doRFDataProc(varargin)", "global TOGGLE_RFDataProc; TOGGLE_RFDataProc = logical(UIState); return;"]) ...
);

% process RF Data
vEvent = VSXEvent('info', 'Proc RF Data', 'process', proc_rf_data, 'seqControl', vSeq);                               

function [vUI, vEvent] = addCustomSaver(vbuf_rx, vSeq)
arguments
    vbuf_rx (1,1) VSXRcvBuffer
    vSeq (1,:) VSXSeqControl = VSXSeqControl.empty
end
%% Process: saving data
save_rf_data = VSXProcess('classname', 'External', 'method', 'RFDataStore');
save_rf_data.Parameters = {
    'srcbuffer','receive',...
    'srcbufnum', vbuf_rx,... 
    'srcframenum',0,...
    'dstbuffer','none'};

vUI = VSXUI( ...
'Control', {'UserB1','Style','VsPushButton','Label', 'SAVE RFData', 'Callback', @doRFDataStore}, ...
'Callback', cellstr(["doRFDataStore(varargin)", "global TOGGLE_RFDataStore; TOGGLE_RFDataStore = true; return;"]) ...
);

% save all RF Data
vEvent = VSXEvent('info', 'Save RF Data', 'process', save_rf_data, 'seqControl', vSeq);

