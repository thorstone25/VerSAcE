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
% [...] = QUPS2VSX(..., 'frames', F) creates F frames per block. This must
% be 1 or an even number. The default is 1.
%
% [...] = QUPS2VSX(..., 'units', 'wavelengths') or 
% [...] = QUPS2VSX(..., 'units', 'mm') sets the units to wavelengths or
% millimeters respectively. The default is 'mm'.
%
% [...] = QUPS2VSX(..., 'sample_mode', samp) uses the sampling mode given by 
% the string samp. The default is "NS200BW".
%
% [...] = QUPS2VSX(..., 'custom_fs', fs) sets a custom sampling frequency
% fs in Hz. This overrides the sampling frequency chosen by Verasonics.
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
% potentially *override* the delays from QUPS, making them innacurate! 
% However, these properties are required for the 'computeTXPD' utility. The
% default is true.
%
% [...] = QUPS2VSX(..., 'multi_rx', false) disables receive multiplexing.
% If there are more elements than channels, an aperture must be selected
% for each `VSXReceive`. For example, to use aperture 33 for all events:
% ```
% vb = QUPS2VSX(...); % create a block
% rxs = [vb.capture.rcv]; % extract all VSXReceives
% [rxs.aperture] = deal(33); % set the aperture on Receive
% vs = vb.link(); % link
% ```
% The default is true when there are more elements than channels.
% 
% Example:
% % Create an UltrasoundSystem
% [xdc, Trans] = TransducerVerasonics('P4-2v');
% seq = Sequence('type', 'FSA', 'numPulse', xdc.numel);
% scan = ScanCartesian('x', -30e-3 : 0.1e-3 : 30e-3, 'z', 0 : 0.1e-3 : 100e-3);
% us = UltrasoundSystem('scan', scan, 'seq', seq, 'xdc', xdc);
%
% % Create the VSXBlock
% vRes = VSXResource();
% [vBlock, chd] = QUPS2VSX(us, Trans, vRes, 'recon_VSX', true);
%
% % Link to create the struct
% vs = vBlock.link(vRes, Trans);
%
% % Save to a MAT-file
% filename = fullfile('MatFiles','qups-vsx.mat');
% save(filename, '-struct', 'vs'); % save the fields of 'vs'
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
    kwargs.vTW (1,:) VSXTW = VSXTW('type', 'parametric', 'Parameters', [us.xdc.fc/1e6, 0.67, 2, 1]);
    kwargs.vTGC VSXTGC {mustBeScalarOrEmpty} = VSXTGC( ...
        'CntrlPts', [0,297,424,515,627,764,871,1000],...
        'rangeMax', max(vecnorm(us.scan.positions,2,1),[],'all') ./ us.lambda ...
        ); % default Time Gain Compensation (TGC)
    kwargs.vTPC VSXTPC {mustBeScalarOrEmpty} = VSXTPC.empty; 
    kwargs.frames (1,1) {mustBeInteger, mustBePositive} = 1;
    kwargs.sample_mode (1,1) string {mustBeMember(kwargs.sample_mode, ["NS200BW", "NS200BWI", "BS100BW",  "BS67BW",  "BS50BW",  "custom"])} = "NS200BW"
    kwargs.custom_fs (1,1) double {mustBeInRange(kwargs.custom_fs, 2.5e6, 62.5e6)}
    kwargs.range (1,2) {mustBeReal, mustBeNonnegative} = [us.scan.zb(1), hypot(max(abs(us.scan.xb)), max(abs(us.scan.zb)))]
    kwargs.recon_VSX (1,1) logical = false
    kwargs.saver_custom (1,1) logical = false
    kwargs.multi_rx (1,1) logical = (isa(xdc, 'struct') && (xdc.numelements > vResource.Parameters.numRcvChannels)) ...
        || (isa(xdc, 'Transducer') && (xdc.numel > vResource.Parameters.numRcvChannels))
    kwargs.set_foci (1,1) logical = true
    kwargs.prf (1,1) double {mustBeNonnegative} = 0 % PRF in Hz
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
M = max(sum(logical(apd),1)); % max number of simultaneously transmitting elements
if M > vResource.Parameters.numTransmit
    warning("QUPS2VSX:multiplexTX", ...
        "Multiplexing the transmit aperture to satisfy the number of transmit channels." ...
        );
    [seq, kwargs.apod] = multiplex(seq, xdc, vResource.Parameters.numTransmit, 'apod', kwargs.apod);
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
if isfield(kwargs, "custom_fs")
    fs_target = 1e-6*kwargs.custom_fs;
else
    fs_target = Trans.frequency * 4;
end
fs_decim = fs_available(argmin(abs(fs_available - fs_target))); % round to nearest supported frequency

switch kwargs.sample_mode
    case "NS200BW", fs_div = 4;
    case "BS100BW", fs_div = 2;
    case "BS67BW",  fs_div = 2 * 0.67;
    case "BS50BW",  fs_div = 1;
    case "custom",  fs_div = 1;
    otherwise, error("Sampling mode '" + kwargs.sample_mode + "' unsupported.");
end
% effective sampling frequency
fs_eff = fs_decim / fs_div; 

% get the output data buffer length
spw = fs_decim / fs_div / Trans.frequency; % samples per wave
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
delay = - seq.delays(xdc) * xdc.fc;
apod  = seq.apodization(xdc);
% delay(~apod) = nan; % deactivate non-active elements, TODO: debug this
t0    = min(delay, [], 1, 'omitnan'); % min over elements
if all(t0(1) == t0, 2), t0 = t0(1); end % scalarize if possible
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
                vTX(i).FocalPt = seq.focus(:,i) ./ lambda; %  when set, {focus,Steer,Origin} are ignored
                if isa(xdc, "TransducerConvex"), % depth is to the radial origin
                    vTX(i).focus = (norm(seq.focus(:,i)-xdc.center) - xdc.radius) ./ lambda; 
                else % array, matrix | TODO: skip on generic
                    vTX(i).focus = seq.focus(3,i) ./ lambda; % depth is to the z == 0 plane
                end
                vTX(i).Origin = xdc.focActive(apod(:,i), 0) ./ lambda; % get the beam origin using this apodization | TODO: skip on generic
            case "PW"
                vTX(i).Steer = deg2rad([seq.angles(i), 0]);
                vTX(i).focus = 0; % focal distance == 0 for PW
                vTX(i).Origin = 0; % = [0 0 0] for PW
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
    % im = median(find(apod(:,i)), 2); % middle element

    % TODO: use computeTXDelays instead?
    % vTX(i).Delay = computeTXDelays(struct(vTX(i))); % requires base workspace variables
end

%% TGC
kwargs.vTGC.Waveform = computeTGCWaveform(kwargs.vTGC, 1e6*Trans.frequency);

%% Rcv
% default
vRcv = VSXReceive('startDepth', dnear, 'endDepth', dfar, 'bufnum', vbuf_rx, 'sampleMode', kwargs.sample_mode);
if isfield(kwargs,'custom_fs'), vRcv.decimSampleRate = fs_decim; end
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
t_puls = round(2 * dfar / Trans.frequency) + 5; % pulse wait time in usec (+5us buffer)
t_frm  = double(t_puls*Mx*seq.numPulse); % total frame time in usec 
wait_for_tx_pulse = VSXSeqControl('command', 'timeToNextAcq', 'argument', t_puls);
if kwargs.prf,
    t_last = round((1e6/kwargs.prf) - (t_frm - t_puls));
    assert(t_last >= t_puls, "VERSACE:QUPS2VSX:UnsatisfiablePRF", ... 
        "A PRF of "+kwargs.prf+" Hz at a range of "+1e3*kwargs.range(end)+" mm for "+Mx+" x "+seq.numPulse+" transmits requires a final wait time of "+t_last+" usec, which subceeds the minimum required pulse wait time of "+t_puls+" usec."...
        +newline+"The maximum satisfiable PRF is "+(1e6/(t_frm - t_puls))+" Hz."...
        );
    wait_for_pulse_sequence  = VSXSeqControl('command', 'timeToNextAcq', 'argument', t_last); % max TTNA is 4190000
else
    wait_for_pulse_sequence = wait_for_tx_pulse; % no PRF requested - use the same wait time
end    
transfer_to_host         = VSXSeqControl('command', 'transferToHost');
no_operation             = VSXSeqControl('command', 'noop', 'argument', 100/0.2); % 'condition', 'Hw&Sw');

%% Event loop
% loop through all events and frames
% ---------- Events ------------- %
vEvent = repmat(VSXEvent(), size(vRcv));

for f = 1:kwargs.frames
    for i = 1:seq.numPulse % each transmit
        for m = 1:Mx % each receive aperture
            vEvent(m,i,f) = VSXEvent( ...
                'tx',vTX(i), 'rcv',vRcv(m,i,f), ...
                'seqControl', wait_for_tx_pulse,...
                'info', "Tx "+i+" - Ap "+m+" - Frame "+f ...
                );
        end
    end

    % transfer data to host using the last event of the frame
    vEvent(m,i,f).seqControl = [wait_for_pulse_sequence, copy(transfer_to_host)]; % modify last acquisition vEvent's seqControl
end
    
%% Add Events per frame 
% make recon image and return to MATLAB
if kwargs.recon_VSX
    [recon_event, vPData] = addReconVSX(scan, vTX, vRcv, vResource, 'multipage', false, 'numFrames', kwargs.frames, 'apod', kwargs.apod);
    recon_event.recon.newFrameTimeout = 50 * t_puls*Mx*seq.numPulse; % wait time in usec ~ heuristic is approx 50 x acquisition time

    vEvent = reshape(vEvent, [1 Mx*seq.numPulse kwargs.frames]); % combine rx multiplexing
    vEvent(:,end+1,:) = recon_event; % assign to end of TX-loop
    vEvent(:,end,:) = copy(vEvent(:,end,:)); % copy to make unique Events for each frame
else
    vPData = VSXPData.empty;
end

%% Add Events at the end of the loop
% vectorize
vev_post = reshape(VSXEvent.empty, [1,0]);

% custom saving process and ui
if kwargs.saver_custom
    [vUI(end+1), vev_post(end+1)] = addSaveRF(vbuf_rx, no_operation);
end

% ---------- Events ------------- %

%% Block
vBlock = VSXBlock('capture', vEvent, 'post', vev_post, 'next', vEvent(1), 'vUI', vUI);
if ~isempty(kwargs.vTPC), vBlock.vTPC = kwargs.vTPC; end

%% Create a template ChannelData object
t0l = - 2 * Trans.lensCorrection; % lens correction in wavelengths
t0p = - vTW.peak; % peak correction in wavelengths
x = zeros([T Mx*seq.numPulse Trans.numelements kwargs.frames, 0], 'single'); % pre-allocated array
chd = ChannelData('data', x, 'fs', 1e6*fs_eff, 't0', (t0 + t0l + t0p)./xdc.fc, 'order', 'TMNF');

%% added External Functions/Callback
% restore warning state
warning(warning_state);


