function [vBlock, vPData, vTrans, vUI, t0, vTW, vTX, vRcv, vRecon, display_image_process, vSeqControl, vTGC, vReconInfo] = QUPS2VSX(us, xdc, vResource, kwargs)
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
    kwargs.vTGC (1,1) VSXTGC = VSXTGC('CntrlPts', [0,297,424,515,627,764,871,1000],...
        'rangeMax', hypot(us.scan.zb(2), us.scan.xb(2)) ./ us.lambda); %-
    kwargs.frames (1,1) {mustBeInteger, mustBePositive} = 1;
end

% squash obj to struct warning
warning_state = warning('off', 'MATLAB:structOnObject');

%% Handle inputs
F = kwargs.frames;

%% Trans
% set sound speed
c0 = us.seq.c0;

if isa(xdc, 'string') % interpret as name
    vTrans.name = char(xdc);
    vTrans.units = char(kwargs.units);
elseif isa(xdc, 'struct') % interpret as the Trans struct
    vTrans = xdc;
elseif isa(xdc, "Transducer") % make custom
    warning("Untested.");
    if     isa(xdc, 'TransducerArray'),     xdctype = 0;
    elseif isa(xdc, 'TransducerConvex'),    xdctype = 1;
    elseif isa(xdc, 'TransducerMatrix'),    xdctype = 2;
    else,                                   xdctype = 2;
        % treat all others as matrix, since it is the most flexible
    end
    imp = copy(us.xdc.impulse);
    imp.fs = 250e6; % sample at 250MHz
    [th, phi] = orientations(xdc); % element orientations (deg)
    pn = positions(xdc); % element positions (m)
    lambda = c0 ./ xdc.fc; % wavelength
    vTrans = struct( ...
        strip('name            '), "custom", ...
        strip('units           '), kwargs.units, ...
        strip('frequency       '), us.xdc.fc / 1e6, ...
        strip('Bandwidth       '), us.xdc.bw / 1e6, ...
        strip('type            '), xdctype, ...
        strip('numelements     '), us.xdc.numel, ...
        strip('ElementPos      '), [pn * 1e3; deg2rad(th); deg2rad(phi)]', ...
        strip('elementWidth    '), us.xdc.width, ...
        strip('IR1wy           '), wv.samples, ...
        strip('connType        '), -1, ... -1 specifies a UTA compatible connecter
        strip('elevationFocusMm'), us.xdc.el_focus * 1e3, ...
        strip('elevationApertureMm'), us.xdc.height * 1e3 ...
        );
    switch xdctype
        case {0, 2, 4}, vTrans.spacing = us.xdc.pitch ./ lambda;
        case {1}, vTrans.radius = us.xdc.radius ./ lambda;
    end
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
dnear = 2 * scan.zb(1); % nearest distance (2-way)
dfar  = 2 * hypot(range(scan.xb), scan.zb(end)); % furthest distance (2-way)

%% PData
vPData = VSXPData();
% vPData.PDelta = [0.5, 0, 0.5];
vPData.PDelta = [scan.dx, 0, scan.dz];
vPData.Size(1) = scan.nz; %-
vPData.Size(2) = scan.nx; %-
vPData.Size(3) = scan.ny; %-
vPData.Origin = [scan.xb(1), scan.yb(1), scan.zb(1)];

% TODO: compute pixel regions
% vPData.Region = computeRegions(struct(vPData));

% %     vResource.DisplayWindow(end+1) = VSXDisplayWindow('ReferencePt', vPData.Origin);
vDisplayWindow = VSXDisplayWindow('ReferencePt', vPData.Origin);
vDisplayWindow.Title = 'L11-5vFlashAngles';
vDisplayWindow.pdelta = scan.dx; % 0.35;
vDisplayWindow.Position = [250, 89.5, scan.nx, scan.nz];
vDisplayWindow.ReferencePt = [vPData(1).Origin(1),0,vPData(1).Origin(3)];   % 2D imaging is in the X,Z plane
vDisplayWindow.numFrames = F;
vDisplayWindow.AxesUnits = 'mm';
vDisplayWindow.Colormap = gray(256);
vResource.DisplayWindow(end+1) = vDisplayWindow;

%% Allocate buffers and Set Parameters

% infer the selected sampling frequency for NS200BW
fs_available = 250 ./ (100:-1:4); % all supported sampling frequencies
fs_decim = fs_available(find(fs_available >= 4 * vTrans.frequency, 1, 'first')); % decimation frequency for NS200BW

% get the output data buffer length
spw = fs_decim / vTrans.frequency; % samples per wave
bufLen = 256 * ceil((dfar - dnear) * spw / 256); % only modulus 256 sample buffer length supported

% T = 256 * 20; %% bufLen;
T = bufLen;

% make new buffers
vbuf_inter = VSXInterBuffer('numFrames', F);
vbuf_im = VSXImageBuffer('numFrames', F);
vbuf_rx   = VSXRcvBuffer('rowsPerFrame', T * us.seq.numPulse,...  %-
    'colsPerFrame', vResource.Parameters.numRcvChannels,...
    'numFrames', F);
vResource.InterBuffer(end+1) = vbuf_inter;
vResource.ImageBuffer(end+1) = vbuf_im;
vResource.RcvBuffer(end+1)   = vbuf_rx;

% set other resource params
vResource.Parameters.speedOfSound = c0;
vResource.Parameters.verbose = 2; % default verbosity


%% TW
vTW = kwargs.vTW;

%% TX
vTX.Origin = [0.0,0.0,0.0]; %%
vTX = copy(repmat(VSXTX(), [1, us.seq.numPulse]));
[vTX.waveform] = deal(vTW);

% get delay and apodization matrices
delay  = - us.seq.delays(us.xdc) * us.xdc.fc;
apod = us.seq.apodization(us.xdc);
t0 = min(delay, [], 1); % min over elements
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
vTGC = kwargs.vTGC;
vTGC.Waveform = computeTGCWaveform(vTGC, 1e6*vTrans.frequency);

%% Rcv
vRcv = VSXReceive();
%     vRcv.aperture = 1;
vRcv.Apod = ones(1, [vResource.Parameters.numRcvChannels]);
vRcv.startDepth = floor(dnear); %% 2
vRcv.endDepth = ceil(dfar); %% 256
vRcv.TGC = kwargs.vTGC(1);
vRcv.bufnum = vResource.RcvBuffer;
vRcv.framenum = 1;
vRcv.acqNum = 1;
vRcv.callMediaFunc = 0;

vRcv = copy(repmat(vRcv,1,us.seq.numPulse));

vRcv(1).callMediaFunc = 1; %%
% - Set event specific Receive attributes.
for i = 1:vResource.RcvBuffer(1).numFrames
    for j = 1:us.seq.numPulse
        vRcv(j).framenum = i;
        vRcv(j).acqNum = j;
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
% We need na ReconInfo structures for na steering angles.
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
display_image_process = VSXProcess();
display_image_process.classname = 'Image';
display_image_process.method = 'imageDisplay';
display_image_process.Parameters = {
    'imgbufnum', vbuf_im,...   % number of buffer to process.
    'framenum',-1,...   % (-1 => lastFrame)
    'pdatanum', vPData,...    % number of PData structure to use v%TODO: handle in linking stage.
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

save_rf_data = VSXProcess();
save_rf_data.classname = 'External';
save_rf_data.method = 'RFDataStore';
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

% loop through all events
% ---------- Events ------------- %
vEvent = copy(repmat(VSXEvent('seqControl', wait_for_tx_pulse), [1 us.seq.numPulse]));
% vEvent = sort(vEvent);

for i = 1:us.seq.numPulse % each transmit
    vEvent(i).info = 'Full aperture.';
    vEvent(i).tx  = vTX(i);
    vEvent(i).rcv = vRcv(i);
    vEvent(i).rcv.acqNum = i;
end

% transfer data to host using the last event
vEvent(i).seqControl = [wait_for_pulse_sequence, transfer_to_host]; % modify last acquisition vEvent's seqControl

% post-processing events and return to MATLAB
vEvent(end+1) = VSXEvent(...
    'info', 'recon and process', ...
    'recon', vRecon, ...
    'process', display_image_process, ...
    'seqControl', return_to_matlab ...
    );

% save RF Data
vEvent(end+1) = VSXEvent(...
    'info', 'Save RF Data', ...
    'process', save_rf_data,...
    'seqControl', no_operation...
    );

% return to start of block
vEvent(end+1) = VSXEvent(...
    'info', 'Jump back',...
    'seqControl', [jump_to_image_start]);

% ------------ Events ------------ %
jump_to_image_start.argument = vEvent(1);

%% ADDED UI
vUI = VSXUI();
vUI.Control =  {'UserB1','Style','VsPushButton','Label', 'SAVE RFData', 'Callback', @doRFDataStore};

%% Block
vBlock = VSXBlock();
vBlock.vsxevent = vEvent;

%% added External Functions/Callback
% restore warning state
warning(warning_state);

end