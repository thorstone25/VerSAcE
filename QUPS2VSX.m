function [vBlock, vPData, vTrans, vUI, chd] = QUPS2VSX(us, xdc, vResource, kwargs)
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
% [...] = QUPS2VSX(us, xdc_name, vResource) uses the VSXResource vResource
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
end

% squash obj to struct warning
warning_state = warning('off', 'MATLAB:structOnObject');

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

% set global params
lambda = c0 / (1e6*vTrans.frequency); % wavelengths

% get the scan region in units of wavelengths so we know the ROI
assert(isa(us.scan, 'ScanCartesian')); % should be ScanCartesian
scan = scale(us.scan, 'dist', 1./lambda);

% get temporal sampling region
dnear = floor(2 * scan.zb(1)); % nearest distance (2-way)
dfar  =  ceil(2 * hypot(range(scan.xb), scan.zb(end))); % furthest distance (2-way)

%% PData
vPData = VSXPData();
% vPData.PDelta = [0.5, 0, 0.5];
vPData.PDelta = [scan.dx, 0, scan.dz];
vPData.Size   = [scan.nz, scan.nx, scan.ny]; %-
vPData.Origin = [scan.xb(1), scan.yb(1), scan.zb(1)];

% TODO: compute pixel regions
% vPData.Region = computeRegions(struct(vPData));

% %     vResource.DisplayWindow(end+1) = VSXDisplayWindow('ReferencePt', vPData.Origin);
vDisplayWindow = VSXDisplayWindow( ...
    'ReferencePt', vPData.Origin, ...
    'Title', 'L11-5vFlashAngles', ...
    'pdelta', scan.dx, ...; % 0.35
    'Position', [250, 89.5, scan.nx, scan.nz], ...
    'ReferencePt', [vPData(1).Origin(1),0,vPData(1).Origin(3)], ... ;   % 2D imaging is in the X,Z plan
    'numFrames', kwargs.frames, ...
    'AxesUnits', 'mm', ...
    'Colormap', gray(256) ...
);
vResource.DisplayWindow(end+1) = vDisplayWindow;

%% Allocate buffers and Set Parameters

% infer the selected sampling frequency for NS200BW
fs_available = 250 ./ (100:-1:4); % all supported sampling frequencies
fs_decim = fs_available(find(fs_available >= 4 * vTrans.frequency, 1, 'first')); % decimation frequency for NS200BW

% get the output data buffer length
spw = fs_decim / vTrans.frequency; % samples per wave
bufLen = 128 * ceil(2 * (dfar - dnear) * spw / 128); % only modulus 128 sample buffer length supported

% T = 256 * 20; %% bufLen;
T = bufLen;

% make new buffers
vbuf_inter = VSXInterBuffer('numFrames', kwargs.frames);
vbuf_im    = VSXImageBuffer('numFrames', kwargs.frames);
vbuf_rx    = VSXRcvBuffer(  'numFrames', kwargs.frames, ...
    'rowsPerFrame', T * us.seq.numPulse,...  %-
    'colsPerFrame', vResource.Parameters.numRcvChannels...
    );

% add buffers to the resources
vResource.InterBuffer(end+1) = vbuf_inter;
vResource.ImageBuffer(end+1) = vbuf_im;
vResource.RcvBuffer(end+1)   = vbuf_rx;

%% TW
% vTW = kwargs.vTW;

%% TX
vTX = copy(repmat(VSXTX('waveform', kwargs.vTW), [1, us.seq.numPulse]));

% get delay and apodization matrices
delay = - us.seq.delays(us.xdc) * us.xdc.fc;
apod  = us.seq.apodization(us.xdc);
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
vRcv = VSXReceive('startDepth', dnear, 'endDepth', dfar, 'bufnum', vbuf_rx);
vRcv.Apod = ones([1, vResource.Parameters.numRcvChannels]);
vRcv.TGC = kwargs.vTGC;

% replicate
vRcv = copy(repmat(vRcv,[us.seq.numPulse, kwargs.frames]));

% - Set event specific Receive attributes.
for i = 1:kwargs.frames
    % move points before (or after?) first receive of the frame
    vRcv(1,i).callMediaFunc = true;
    for j = 1:us.seq.numPulse
        vRcv(j,i).framenum = i;
        vRcv(j,i).acqNum = j;
    end
end

%% Recon
vRecon = VSXRecon('pdatanum', vPData, 'IntBufDest', vbuf_inter, 'ImgBufDest', vbuf_im);
%{
vRecon.senscutoff = 0.6;
vRecon.pdatanum = vPData;
vRecon.rcvBufFrame = -1;
vRecon.IntBufDest = vbuf_inter;
vRecon.IntBufDestFrm = 1;
vRecon.ImgBufDest = vbuf_im;
vRecon.ImgBufDestFrm = -1;
%}

% Define ReconInfo structures.
% We need 1 ReconInfo structures for each transmit
vReconInfo = copy(repmat(VSXReconInfo('mode', 'accumIQ'), [1,us.seq.numPulse]));  % default is to accumulate IQ data.

% - Set specific ReconInfo attributes.
if us.seq.numPulse <= 1
    vReconInfo(1).mode = 'replaceIntensity';
else
    vReconInfo(1).mode = 'replaceIQ'; % replace IQ data
    for j = 1:us.seq.numPulse  % For each row in the column
        vReconInfo(j).txnum = vTX(j);
        vReconInfo(j).rcvnum = vRcv(j);
    end
    vReconInfo(us.seq.numPulse).mode = 'accumIQ_replaceIntensity'; % accum and detect
end

% associate the reconinfo
vRecon.RINums = vReconInfo;

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

save_rf_data = VSXProcess('classname', 'External', 'method', 'RFDataStore');
save_rf_data.Parameters = {
    'srcbuffer','receive',...
    'srcbufnum', vbuf_rx,... 
    'srcframenum',0,...
    'dstbuffer','none'};


%% SeqControl
jump_to_image_start      = VSXSeqControl('command', 'jump', 'argument', 1);
wait_for_tx_pulse        = VSXSeqControl('command', 'timeToNextAcq', 'argument', 160);
wait_for_pulse_sequence  = VSXSeqControl('command', 'timeToNextAcq', 'argument', 19040);
return_to_matlab         = VSXSeqControl('command', 'returnToMatlab');
transfer_to_host         = VSXSeqControl('command', 'transferToHost');
no_operation             = VSXSeqControl('command', 'noop', 'argument', 100/0.2); % 'condition', 'Hw&Sw');

%% Event

% loop through all events and frames
% ---------- Events ------------- %
vEvent = copy(repmat(VSXEvent('seqControl', wait_for_tx_pulse), [us.seq.numPulse, kwargs.frames]));

for f = 1:kwargs.frames
    for i = 1:us.seq.numPulse % each transmit
        vEvent(i,f).info = 'Full aperture.';
        vEvent(i,f).tx  = vTX(i);
        vEvent(i,f).rcv = vRcv(i,f);
        vEvent(i,f).rcv.acqNum = i;
    end

    % transfer data to host using the last event
    vEvent(i,f).seqControl = [wait_for_pulse_sequence, transfer_to_host]; % modify last acquisition vEvent's seqControl

    % actions per frame go here
    % post-processing events and return to MATLAB
    vEvent(i+1,f) = VSXEvent(...
        'info', 'recon and process', ...
        'recon', vRecon, ...
        'process', display_image_process, ...
        'seqControl', return_to_matlab ...
        );
end

% first event 
acq_start_evnt = vEvent(1,1);

% vectorize
vEvent = vEvent(:);

% save all RF Data
vEvent(end+1) = VSXEvent(...
    'info', 'Save RF Data', ...
    'process', save_rf_data,...
    'seqControl', no_operation...
    );

% return to start of block
vEvent(end+1) = VSXEvent(...
    'info', 'Jump back',...
    'seqControl', jump_to_image_start);

% ------------ Events ------------ %
jump_to_image_start.argument = acq_start_evnt;

%% ADDED UI
vUI = VSXUI();
vUI.Control =  {'UserB1','Style','VsPushButton','Label', 'SAVE RFData', 'Callback', @doRFDataStore};
vUI.Callback = cellstr(["doRFDataStore(varargin)", "global toggle; toggle true; return;"]);

%% Block
vBlock = VSXBlock('vsxevent', vEvent);

%% Create a template ChannelData object
x = zeros([T us.seq.numPulse vTrans.numelements kwargs.frames, 0], 'int16');
chd = ChannelData('data', x, 'fs', fs_decim, 't0', t0, 'order', 'TMNF');

%% added External Functions/Callback
% restore warning state
warning(warning_state);

end
