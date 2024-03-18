 %#ok<*UNRCH> 
 
 %% Create different system configurations
c0 = 1500; % sound speed
% xdc_name = "L12-3v"; % choose transducer
% xdc_name = "L12-5 50mm"; % choose transducer
xdc_name = "P4-2v";
Trans = computeTrans(struct("name", char(xdc_name), 'units', 'mm')); % transducer
xdc = Transducer.Verasonics(Trans);
uss = UltrasoundSystem('xdc', xdc); 
uss.seq.c0 = c0;
uss.scan.zb(2) = 150 * uss.lambda;
uss.scan.xb = 1e-3 * [-40 40];
[uss.scan.dx, uss.scan.dz] = deal(uss.lambda / 2);

apd_center = false; % use only ceneter aperture

% sequences
pf = [0;0;50e-3] + [1e-3;0;0] .* (-20 : 5 : 20);
seqfsa = Sequence(     'type', 'FSA', 'c0', c0, 'numPulse', xdc.numel);
seqpw = SequenceRadial('type', 'PW' , 'c0', c0, 'angles', -25 : 0.5 : 25); 
seqfc = Sequence(      'type', 'FC' , 'c0', c0, 'focus', pf);

% spatial downsampling
Ds = 4; 
seqfsadv = Sequence(  'type', 'FSA', 'c0', c0, 'numPulse', xdc.numel / Ds);
seqfsadv.apd = zeros([xdc.numel, seqfsadv.numPulse]);
seqfsadv.apd(1:Ds:end,:) = eye(xdc.numel / Ds);

uss = copy(repmat(uss, [1,4]));
[uss.seq] = deal(seqfsa, seqpw, seqfc, seqfsadv);

% apodization schemes
apod0 = cell([1, numel(uss)]);
for i = 1:numel(uss)
    switch uss(i).seq.type
        case "FSA", apod0{i} = uss(i).apAcceptanceAngle(45);
        case "PW",  apod0{i} = swapdim(uss(i).apTxParallelogram(uss(i).seq.angles),4,5);
        case "FC",  apod0{i} = swapdim(uss(i).apMultiline(),4,5);
        case "DV",  apod0{i} = uss(i).apAcceptanceAngle(45);
    end
end

%%
% selection
seq_ind = [2 4]; % sequence index
[us, apod] = deal(copy(uss(seq_ind)), apod0(seq_ind)); % choose pulse sequence template
if false && apd_center
    N = min(128, us.xdc.numel); % active aperture size
    apd = circshift([ones(1,N), zeros(1,us.xdc.numel-N)], (us.xdc.numel-N)/2)';
    us.seq.apodization_ = repmat(apd, [1 us.seq.numPulse]);
end

% constant resources
vres = VSXResource(); % system-wide resource
vTW = VSXTW('type','parametric', 'Parameters', [Trans.frequency, 0.67, 1, 1]); % tx waveform
vTPC = VSXTPC('name','Default', 'hv', Trans.maxHighVoltage); % max power

% make blocks
for i = 1:numel(us)
[vb(i), chd(i)] = QUPS2VSX(us(i), Trans, vres, "apod", apod{i} ...
    ,"frames", 4 ...
    ,'vTW', vTW, 'vTPC', vTPC ...
    ,'recon_VSX', any(i == 1) ...
    ,'saver_custom', false ...
    ,'set_foci', true ...
    ,'range', [0 us(i).scan.zb]*1e-3 ... in m
    ,'range', [0 200]*1e-3 ... in m
); % make VSX block
end

% connect blocks
vb(1).next = vb(2).capture(1);
vb(2).next = vb(1).capture(1);

% set peak cut-off for VSX imaging
vec = @(x)x(:);
evs = arrayfun(@(x) {vec(x.capture(:))'}, vb);
evs = [evs{:}];
txs = [evs.tx];
[txs.peakCutOff] = deal(2);

% Terminate when done acquiring (not working ...)
% vb.post(end+1) = VSXEvent('info', 'Stop', 'seqControl', [...
%     VSXSeqControl('command', 'stop'), VSXSeqControl('command', 'returnToMatlab')...
% ]);
% vb.next = VSXEvent.empty;

%{
% [vb(1), chd(1)] = QUPS2VSX(uss(1), Trans, vres, "frames", 1, 'vTW', vTW); % make VSX block
% [vb(2), chd(2)] = QUPS2VSX(uss(2), Trans, vres, "frames", 4, 'vTW', vTW); % make VSX block
% [vb.next] = deal(vb(2).capture(1), vb(1).capture(1)); % start at beginning of alternate sequence
%}
% DEBUG: test the manual receive delays
%{
for i = 1:numel(seq_ind)
    [~, tau_rx, tau_tx] = bfDAS(uss(seq_ind(i)), chd(i), 'delay_only', true);
    vRecon = unique([vb(i).capture.recon]); % find Recon (exactly 1 exists)
    if ~isscalar(vRecon), continue; end
    setVSXLUT(vRecon, tau_rx, tau_tx - swapdim(chd(i).t0,chd(i).mdim,5), uss(seq_ind(i)).xdc.fc);% broken for Vantage 4.3
end
%}

%% convert to VSX structures
global vs; % make global to access in within function
vs = link(vb, vres, Trans, 'TXPD', true); % link
pt1; vs.Media = Media; % add simulation media

% DEBUG: test the manual receive delays
% [~, tau_rx, tau_tx] = bfDAS(us, chd, 'delay_only', true);
% [vs.Recon, vs.ReconInfo] = setVSXLUT(vs.Recon, vs.ReconInfo, vs.PData, tau_rx, tau_tx + swapdim(chd.t0,chd.mdim,5), us.xdc.fc);

% force in simulation mode for testing
vs.Resource.Parameters.simulateMode = 1; % 1 to force simulate mode, 0 for hardware

%% save 
filename = char(fullfile("MatFiles","qups-vsx.mat")); 
save(filename, '-struct', 'vs');

% save
save(fullfile("MatFiles","qups-conf.mat"), "us", "chd");

% set output save directory
global VSXOOD_SAVE_DIR;
VSXOOD_SAVE_DIR = fullfile(pwd, 'tmp');
if ~exist(VSXOOD_SAVE_DIR, 'dir'), mkdir(VSXOOD_SAVE_DIR); end

% clear external functions
clear RFDataImg RFDataProc RFDataStore RFDataCImage imagingProc cEstFSA_RT;

% VSX;

%% 
