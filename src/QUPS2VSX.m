function [vBlock, chd, Trans] = QUPS2VSX(us, xdc, vResource, kwargs)
% QUPS2VSX - Construct a VSXBlock 
%
% vBlock = QUPS2VSX(us) converts the UltrasoundSystem us into a 
% VSXBlock vBlock. This can be used to generate a Verasonics compatible 
% configuration struct.
% 
% [vBlock, chd] = QUPS2VSX(...) additionally returns a template ChannelData
% chd with data sizing matching the output data in dimensions 1:4 and 
% size 0 in dimension 5. The data can be preallocated with 
% `chd.data(:,:,:,:,1) = 0;`
% 
% [vBlock, chd, Trans] = QUPS2VSX(...) additionally returns a Verasonics Trans
% struct Trans.
%
% [...] = QUPS2VSX(us, xdc) where xdc is a string uses the named transducer
% xdc and the `computeTrans` utility to generate the transducer struct.
% If xdc is a struct, it is assumed to be the verasonics Trans struct and
% the `computeTrans` utility is used to complete the definition.
% If xdc is a Transducer, a custom transducer definition is generated. The 
% default is us.xdc.
%
% [...] = QUPS2VSX(us, xdc, vResource) uses the VSXResource vResource
% instead of creating a new VSXResource. You must specify and retain a 
% single vResource whenever buffers are necessary or when combining 
% multiple blocks.
%
% [...] = QUPS2VSX(..., 'frames', F) creates F frames per block. The default 
% is 1.
%
% [...] = QUPS2VSX(..., 'units', 'wavelengths') or 
% [...] = QUPS2VSX(..., 'units', 'mm') sets the units to wavelengths or
% millimeters respectively. The default is 'mm'.
%
% [...] = QUPS2VSX(..., 'sample_mode', samp) uses the sampling mode given by 
% the string samp. The default is "NS200BW".
%
% [...] = QUPS2VSX(..., 'custom_fs', fs) sets a custom sampling frequency fs. 
% This overrides the sampling frequency chosen by Verasonics.
%
% [...] = QUPS2VSX(..., 'recon_VSX', true) adds basic image reconstruction 
% via Verasonics' Recon structures. The default is false.
%
% [...] = QUPS2VSX(..., 'recon_custom', true) adds basic image reconstruction 
% via QUPS. The default is false.
%
% [...] = QUPS2VSX(..., 'range', [lb ub]) sets the data acquisition range
% from lb to ub in meters. The default is determined by the
% closest/furthest pixel of the scan.
%
% [...] = QUPS2VSX(..., 'set_foci', true) sets the Verasonics 'Origin',
% 'Steer', and 'focus' properties for each 'TX' struct. This appears to
% *override* the delays from QUPS, making them innacurate! However, these
% properties are required for the 'computeTXPD' utility. The default is
% false.
%
% Example:
% % Create an UltrasoundSystem
% xdc = TransducerArray.L11_5v();
% seq = Sequence('type', 'FSA', 'numPulse', xdc.numel);
% scan = ScanCartesian('x', -20e-3 : 0.1e-3 : 20e-3, 'z', 0 : 0.1e-3 : 40e-3);
% us = UltrasoundSystem('scan', scan, 'seq', seq, 'xdc', xdc);
%
% % Create the VSXBlock
% vRes = VSXResource();
% [vBlock, chd, Trans] = QUPS2VSX(us, "L11-5v", vRes, 'recon_VSX', true);
%
% % Link to create the struct
% s = vBlock.link(vRes);
% s.Trans = Trans; % include the Trans struct
%
% % Save to a MAT-file
% filename = fullfile('MatFiles','qups-vsx.mat');
% save(filename, '-struct', 's'); % save the fields of 's'
%
% % Launch VSX
% VSX;
% 
% See also VSXBLOCK/LINK
arguments
    us (1,1) UltrasoundSystem
    xdc (1,1) {mustBeA(xdc, ["string", "struct", "Transducer"])} = us.xdc % transducer name
    vResource (1,1) VSXResource = VSXResource()
    kwargs.apod (:,:,:,:,1) = 1; % image apodization per tx
    kwargs.units (1,1) string {mustBeMember(kwargs.units, ["mm", "wavelengths"])} = "mm"
    kwargs.vTW (1,:) VSXTW = VSXTW('type', 'parametric', 'Parameters', [us.xdc.fc/1e6, 0.67, 1, 1]);
    kwargs.vTGC VSXTGC {mustBeScalarOrEmpty} = VSXTGC( ...
        'CntrlPts', [0,297,424,515,627,764,871,1000],...
        'rangeMax', hypot(us.scan.zb(2), us.scan.xb(2)) ./ us.lambda ...
        ); % default Time Gain Compensation (TGC)
    kwargs.vTPC VSXTPC {mustBeScalarOrEmpty} = VSXTPC.empty; 
    kwargs.frames (1,1) {mustBeInteger, mustBePositive} = 1;
    kwargs.sample_mode (1,1) string {mustBeMember(kwargs.sample_mode, ["NS200BW", "NS200BWI", "BS100BW",  "BS67BW",  "BS50BW",  "custom"])} = "NS200BW"
    kwargs.custom_fs (1,1) double {mustBeInRange(kwargs.custom_fs, 2.5e6, 62.5e6)}
    kwargs.range (1,2) {mustBeReal, mustBeNonnegative} = [us.scan.zb(1), hypot(max(abs(us.scan.xb)), max(abs(us.scan.zb)))]
    kwargs.recon_VSX (1,1) logical = false
    kwargs.recon_custom (1,1) logical = false
    kwargs.recon_custom_delays (1,1) logical = false
    kwargs.saver_custom (1,1) logical = false
    kwargs.multi_rx (1,1) logical = (isa(xdc, 'struct') && (xdc.numelements > vResource.Parameters.numRcvChannels)) ...
        || (isa(xdc, 'Transducer') && (xdc.numel > vResource.Parameters.numRcvChannels))
    kwargs.set_foci (1,1) logical = false
end

% squash obj to struct warning
warning_state = warning('off', 'MATLAB:structOnObject');

% init
vUI     = reshape(VSXUI.empty   , [1 0]);

% tx sequence
seq = us.seq;

%% Trans
% set sound speed
c0 = seq.c0;
vResource.Parameters.speedOfSound = c0;

if isa(xdc, 'string') % interpret as name
    Trans = struct('name', char(xdc), 'units', char(kwargs.units));
elseif isa(xdc, 'struct') % interpret as the Trans struct
    Trans = xdc;
elseif isa(xdc, "Transducer") % make custom
    Trans = xdc.QUPS2VSX();
else
    error("Unrecognized input for xdc.")
end

% have VSX compute remaining properties
Trans = computeTrans(Trans);

% convert to QUPS
xdc = Transducer.Verasonics(Trans);

% transmit multiplexing - override the requested sequence if (clearly) 
% unsatisfiable
% TODO: validate with Trans.connector
apd = seq.apodization(xdc);
M = max(sum(apd,1)); % max number of simultaneously transmitting elements
if M > vResource.Parameters.numTransmit
    warning("QUPS2VSX:multiplexTX", ...
        "Multiplexing the transmit aperture to satisfy the number of transmit channels." ...
        );
    seq = multiplex(seq, xdc, vResource.Parameters.numTransmit);
end

% receive multiplexing factor
if kwargs.multi_rx
    Mx = ceil(xdc.numel / vResource.Parameters.numRcvChannels);
else
    Mx = 1;
end

% set global params
lambda = c0 / (1e6*Trans.frequency); % wavelengths

% get the scan region in units of wavelengths so we know the ROI
assert(isa(us.scan, 'ScanCartesian')); % should be ScanCartesian
scan = scale(us.scan, 'dist', 1./lambda); % convert to wavelengths

% get temporal sampling region (wavelengths)
% dnear = floor(2 * scan.zb(1)); % nearest distance (2-way)
% dfar  =  ceil(2 * hypot(range(scan.xb), scan.zb(end))); % furthest distance (2-way)
[dnear, dfar] = deal(kwargs.range(1)./lambda, kwargs.range(2)./lambda);

%% TW
vTW = arrayfun(@(tw) tw.computeTWWaveform(Trans, vResource, kwargs.vTPC), kwargs.vTW); % fill out the waveform (and TPC)

%% Allocate buffers and Set Parameters
fs_available = 250 ./ (100:-1:4); % all supported sampling frequencies (MHz)

% decimation frequency as per sampleMode
switch kwargs.sample_mode
    case "NS200BW", fs_target = Trans.frequency * 4;
    case "BS100BW", fs_target = Trans.frequency * 2;
    case "BS67BW",  fs_target = Trans.frequency * 2 * 0.67;
    case "BS50BW",  fs_target = Trans.frequency * 1;
    case "custom",  fs_target = kwargs.custom_fs;
    otherwise, error("Sampling mode '" + kwargs.sample_mode + "' unsupported.");
end
% fs_decim = fs_available(find(fs_available >= fs_target, 1, 'first')); % get first frequency above limit
fs_decim = fs_available(argmin(abs(fs_available - fs_target))); % round to nearest supported frequency

% get the output data buffer length
spw = fs_decim / Trans.frequency; % samples per wave
T = 128 * ceil(2 * (dfar - dnear - 1/spw) * spw / 128); % only modulus 128 sample buffer length supported

% make new rx buffer
vbuf_rx    = VSXRcvBuffer(  'numFrames', kwargs.frames, ...
    'rowsPerFrame', T * Mx * seq.numPulse, ... 
    'colsPerFrame', vResource.Parameters.numRcvChannels...
    );

% add buffers to the resources
vResource.RcvBuffer(end+1)   = vbuf_rx;

%% TX
vTX = copy(repmat(VSXTX(), [1, seq.numPulse]));
b = ~isscalar(vTW); % whether vTW is an array
for i = 1:numel(vTX), vTX(i).waveform = vTW(i*b+~b); end % index=i or index=1

% get delay and apodization matrices
posn  = xdc.positions();
delay = - seq.delays(xdc) * xdc.fc;
apod  = seq.apodization(xdc);
% delay(~apod) = nan; % deactivate non-active elements, TODO: debug this
t0    = min(delay, [], 1, 'omitnan'); % min over elements
delay = delay - t0; % don't include 0
if max(delay,[],'all','omitnan') / xdc.fc > 45.5e-6
    error("QUPS2VSX:unsupportedSequence", ...
        "The delays exceed the maximum supported delay of " + 45.5 + " us." ...
        );
end

% - Set event specific TX attributes.
for i = 1:seq.numPulse
    % NOTE: Verasonics overrides users delays when the {'focus','Steer','Origin'} properties are deleted!!!!!
    % These should be left at 0 to ensure that they exist, and the QUPS delays are used.
    % There may be conflicts with TXPD and the way that VSX computes minimum delays
    if kwargs.set_foci % set beamforming geometry
        switch seq.type
            case {"FC","DV","VS"} % focused / diverging / virtual source (assumed focused)
                % TODO: set for arb. Transducer types
                if ~any(arrayfun(@(t) isa(xdc,t), ["TransducerArray", "TransducerMatrix"]))
                    warning("Focal sequence 'Origin' property assumes a planar transducers.");
                end
                vTX(i).FocalPt = seq.focus(:,i) ./ lambda;
                vTX(i).focus   = seq.focus(3,i) ./ lambda; % set z value
                vTX(i).Origin  = seq.focus(:,i) ./ lambda .* [1 1 0]'; % set to z=0
            case "PW"
                vTX(i).Steer = deg2rad([seq.angles(i), 0]);
                % vTX(i).focus = 0; % focal distance == 0 for PW
            case "FSA"
                % vTX(i).focus = 0; % focal distance == 0 for compliance
                %do nothing
        end
    end

    % set delays and apodization
    vTX(i).Apod  =           apod(:,i)' ;
    vTX(i).Delay = nan2zero(delay(:,i)');

    % attempt to find the virtual origin as center of active aperture
    % approximate for non-planar transducers
    im = median(find(apod(:,i)), 2); % middle element
    vTX(i).Origin = mean(posn(:,[floor(im), ceil(im)]),2); % 

    % TODO: use computeTXDelays instead?
    % vTX(i).Delay = computeTXDelays(struct(vTX(i))); % requires base workspace variables
end

%% TGC
kwargs.vTGC.Waveform = computeTGCWaveform(kwargs.vTGC, 1e6*Trans.frequency);

%% Rcv
% default
vRcv = VSXReceive('startDepth', dnear, 'endDepth', dfar, 'bufnum', vbuf_rx, 'sampleMode', kwargs.sample_mode);
if kwargs.sample_mode == "custom", vRcv.decimSampleRate = fs_decim; end
vRcv.TGC = kwargs.vTGC;

% replicate
vRcv = copy(repmat(vRcv,[Mx, seq.numPulse, kwargs.frames]));

% set receive apodization - use "Dynamic HVMux Apertures"
N = Trans.numelements; % total elements
M = N / Mx; % active aperture size (M == N if Mx == 1 for no multiplexing)
ap = [ones([1, M]), zeros([1, N - M])]; % active aperture
for i = 1:Mx % for each multiplexed aperture
    [vRcv(i,:).Apod] = deal(circshift(ap, M*(i-1))); % shift to that portion
end

% - Set event specific Receive attributes.
for f = 1:size(vRcv,3)
    % move points before (or after?) first receive of the frame
    vRcv(1,1,f).callMediaFunc = true;
    for i = 1:size(vRcv,2)
        for m = 1:Mx
            vRcv(m,i,f).framenum = f;
            vRcv(m,i,f).acqNum   = sub2ind(size(vRcv,1:2), m, i);
        end
    end
end

%% SeqControl
t_puls = round(T / fs_decim) + 50; % pulse wait time in usec
t_frm = t_puls*Mx*seq.numPulse*1.2; % frame wait time in usec
wait_for_tx_pulse        = VSXSeqControl('command', 'timeToNextAcq', 'argument', t_puls);
wait_for_pulse_sequence  = VSXSeqControl('command', 'timeToNextAcq', 'argument', t_frm ); % max TTNA is 4190000
transfer_to_host         = VSXSeqControl('command', 'transferToHost');
no_operation             = VSXSeqControl('command', 'noop', 'argument', 100/0.2); % 'condition', 'Hw&Sw');

%% Event loop
% loop through all events and frames
% ---------- Events ------------- %
vEvent = copy(repmat(VSXEvent('seqControl', wait_for_tx_pulse), size(vRcv)));

for f = 1:kwargs.frames
    for i = 1:seq.numPulse % each transmit
        for m = 1:Mx % each receive aperture
            vEvent(m,i,f) = VSXEvent('info',"Tx "+i+" - Ap "+m+" - Frame "+f,'tx',vTX(i),'rcv',vRcv(m,i,f));
        end
    end

    % transfer data to host using the last event of the frame
    vEvent(m,i,f).seqControl = [wait_for_pulse_sequence, copy(transfer_to_host)]; % modify last acquisition vEvent's seqControl
end
vEvent = reshape(vEvent, [prod(size(vEvent,1:2)), size(vEvent,3)]); % combine rx multiplexing
    
%% Add Events per frame 
% make recon image and return to MATLAB
if kwargs.recon_VSX
    return_to_matlab = VSXSeqControl('command', 'returnToMatlab');
    [recon_event] = addVSXRecon(scan, vResource, vTX, vRcv, return_to_matlab, 'multipage', false, 'numFrames', 1+0*kwargs.frames, 'apod', kwargs.apod);
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
    [vUI(end+1), vev_post(end+(1:2))] = addCustomRecon(scan, vResource, vbuf_rx, vPData, no_operation);
end

% custom saving process and ui
if kwargs.saver_custom
    [vUI(end+1), vev_post(end+1)] = addCustomSaver(vbuf_rx, no_operation);
end

% ---------- Events ------------- %

%% Block
vBlock = VSXBlock('capture', vEvent, 'post', vev_post, 'next', vEvent(1), 'vUI', vUI);
if ~isempty(kwargs.vTPC), vBlock.vTPC = kwargs.vTPC; end

%% Create a template ChannelData object
% TODO: call computeTWWaveform to get TW.peak correction
t0l = - 2 * Trans.lensCorrection; % lens correction in wavelengths
t0p = - vTW.peak; % peak correction in wavelengths
x = zeros([T Mx*seq.numPulse Trans.numelements kwargs.frames, 0], 'single'); % pre-allocated array
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
    kwargs.apod (:,:,:,:,:) {mustBeNumericOrLogical} = 1
end

%% Validate sizing
Tx = prod(size(vRcv,1:2));
assert(all(any(size(kwargs.apod,1:5) == [[scan.size,Tx,kwargs.numFrames]; ones(1,5)],1)), ...
    "The size of the apodization(" +...
    join(string(size(kwargs.apod,1:5)),", ") +...
    ") must be compatible with the scan (" + ...
    join(string(scan.size),", ") +...
    ") and the number of transmits (" + Tx + ")." ...
);

%% ReconInfo
% We need 1 ReconInfo structures for each transmit
vReconInfo = copy(repmat(VSXReconInfo('mode', 'accumIQ'), size(vRcv,1:2)));  % default is to accumulate IQ data.

% - Set specific ReconInfo attributes.
for i = 1:size(vReconInfo,1) % for each rx of the reconinfo object (multiplexing)
for j = 1:size(vReconInfo,2) % for each tx of the reconinfo object
    [vReconInfo(i,j,:).txnum    ] = deal(vTX(   j)); % tx
    [vReconInfo(i,j,:).rcvnum   ] = deal(vRcv(i,j)); % rx
    [vReconInfo(i,j,:).regionnum] = deal(       j );% region index TODO: VSXRegion
    if kwargs.multipage
        [vReconInfo(i,j,:).pagenum] = deal(sub2ind(size(vReconInfo),i,j)); % page number
    end
end
end
vReconInfo( 1 ).mode = 'replaceIQ'; % on first tx, replace IQ data
vReconInfo(end).mode = 'accumIQ_replaceIntensity'; % on last tx, accum and "detect"

if isscalar(vReconInfo), vReconInfo.mode = 'replaceIntensity'; end % 1 tx

%% PData
% create the pixel data
vPData = VSXPData.QUPS(scan);

% get the apodization in it's full size
ap = logical(kwargs.apod); % mask
Psz = [scan.size prod(size(vRcv,1:2)), size(vRcv,3)]; % full sizs (I x Tx x F)
ap = repmat(ap, Psz ./ size(ap,1:5)); % explicit brodcast

% convert to address / count format
aps = num2cell(ap, 1:3); % pack pixels per cell
cnt  = cellfun(@nnz,   aps,  "UniformOutput", false); % number of pixels
addr = cellfun(@find,  aps,  "UniformOutput", false); % linear address (1-based)
addr = cellfun(@int32, addr, "UniformOutput", false); % to int
addr = cellfun(@(x)x-1,addr, "UniformOutput", false); % to 0-based

% create a default region for each transmit
% TODO: VSXRegion
% vPData.Region = repmat(struct('Shape', struct('Name','Custom'), 'PixelsLA', int32(0:scan.nPix-1), 'numPixels', scan.nPix), [1 numel(vReconInfo)]);
vPData.Region = cellfun(@(LA, NP) struct('Shape', struct('Name','Custom'), 'PixelsLA', LA, 'numPixels', NP), addr(1:Tx), cnt(1:Tx));

% fill out the regions
% TODO: VSXRegion
vPData.Region = computeRegions(struct(vPData));

%% Make Buffers
% if inter(mediate) buffers needed, create one
if true || any(contains([vReconInfo.mode], "IQ"))
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
vRecon = VSXRecon('pdatanum', vPData, 'IntBufDest', vbuf_inter, 'ImgBufDest', vbuf_im, "RINums", vReconInfo(:));

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
vbuf_inter = VSXInterBuffer.fromPData(vPData, 'numFrames', kwargs.numFrames); 
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
    'Control', {'UserB3','Style','VsToggleButton','Label', 'Image RFData', 'Callback', str2func("do"+nm)}, ...
    'Statement', cellstr(["global TOGGLE_"+nm+"; TOGGLE_"+nm+" = false; return;"]), ... init
    'Callback', cellstr(["do"+nm+"(varargin)", "global TOGGLE_"+nm+"RFDataImg; TOGGLE_"+nm+" = logical(UIState); return;"]) ...
    );

%% Add to required Resource buffer
if ~isempty(vbuf_im       ), vResource.ImageBuffer(end+1)    = vbuf_im       ; end
if ~isempty(vbuf_inter    ), vResource.InterBuffer(end+1)    = vbuf_inter    ; end
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

