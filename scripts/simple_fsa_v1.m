%% Originally modified from test_sequences.m
% Description: FSA sequence for Verasonics L12-3V
% To run, you will need to have QUPS, VSX-OOD, and Vantage added to path!
% To run, you will need to be in the Vantage working directory.
 
%% Create different system configurations
c0 = 1540; % sound speed
xdc_name = "L12-3v"; % choose transducer
% xdc_name = "L12-5 50mm"; % choose transducer
% xdc_name = "P4-2v";
Trans = computeTrans(struct("name", char(xdc_name), 'units', 'mm')); % transducer
xdc = Transducer.Verasonics(Trans);
uss = UltrasoundSystem('xdc', xdc); 
uss.seq.c0 = c0;
uss.scan.zb(2) = 250 * uss.lambda;
uss.scan.xb = 1e-3 * [-40 40];
[uss.scan.dx, uss.scan.dz] = deal(uss.lambda / 2);

apd_center = false; % use only ceneter aperture

% sequences
seqfsa = Sequence('type', 'FSA', 'c0', c0, 'numPulse', xdc.numel);
seqpw = SequenceRadial('type', 'PW' , 'c0', c0, 'angles', -25 : 0.5 : 25); 

% spatial downsampling
Ds = 4; 
seqfsadv = Sequence(  'type', 'FSA', 'c0', c0, 'numPulse', xdc.numel / Ds);
seqfsadv.apd = zeros([xdc.numel, seqfsadv.numPulse]);
seqfsadv.apd(1:Ds:end,:) = eye(xdc.numel / Ds);

uss = copy(repmat(uss, [1,1]));
% [uss.seq] = deal(seqfsa, seqpw, seqfsadv);
% [uss.seq] = deal(seqfsa, seqpw, seqfsadv);
uss.seq = seqfsa;

% apodization schemes
apod0 = cell([1, numel(uss)]);
for i = numel(uss):-1:1
     apod0{i} = uss(i).apAcceptanceAngle(45);
%     switch i
%       case 1, apod0{i} = uss(i).apAcceptanceAngle(45);
      %case 2, apod0{i} = swapdim(apTxParallelogram(uss(i),uss(i).seq.angles, [-15 15]),4,5);
      %case 3, apod0{i} = swapdim(uss(i).apMultiline(),4,5);
      %case 3, apod0{i} = sub(uss(i).apAcceptanceAngle(45),1:Ds:uss(i).xdc.numel,4);
%     end
end

%%
% selection
seq_ind = [1]; % sequence index
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

% % make blocks
% for i = 1:numel(us)
% [vb(i), chd(i)] = QUPS2VSX(us(i), Trans, vres ... , 'apod', apod{i} ...
%     ,"frames", 1 ...
%     ,'vTW', vTW, 'vTPC', vTPC ...
%     ,'recon_VSX', true ...
%     ,'saver_custom', true ...
%     ,'set_foci', true ...
%     ,'range', [0 us(i).scan.zb]*1e-3 ... in m
%     ,'range', [0 200]*1e-3 ...
% ); % make VSX block
% end

% make blocks
for i = 1:numel(us)
[vb(i), chd(i)] = QUPS2VSX(us(i), Trans, vres ... , 'apod', apod{i} ...
    ,"frames", 1 ...
    ,'vTW', vTW, 'vTPC', vTPC ...
    ,'recon_VSX', true ...
    ,'saver_custom', true ...
    ,'range', [0 200]*1e-3 ... in m
); % make VSX block
end

% connect blocks
% vb(1).next = vb(2).capture(1);
% vb(2).next = vb(1).capture(1);

% set peak cut-off for VSX imaging
vec = @(x)x(:);
evs = arrayfun(@(x) {vec(x.capture(:))'}, vb);
evs = [evs{:}];
txs = [evs.tx];
[txs.peakCutOff] = deal(2 / 128);

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

% %% Callback function pre-processing
% No pilot pulse indices (just want to beamform the data later)
ppi = [];

% whether tx multiplexed
tx_multi = max(sum(logical(us.seq.apodization(us.xdc)),1)) > 128;

% pre-processing indexing
evi = find(startsWith({vs.Event.info}, "Tx ")); % events with transmits
txi = double(string(extractBetween({vs.Event(evi).info}, "Tx ", " - Ap"))); % physical transmit indices
if tx_multi, txi = ceil(txi/2); end % to match duplication
bfi = find(~ismember(txi, ppi)); % beamforming event indices

% remove pilot pulses from Sequence/ChannelData definitions
[usv, chdv, us, chd] = deal(us, chd, copy(us), copy(chd)); % store old, save new
us.seq = uss.seq(seq_ind); % final sequence

sz = size(chd.data);
sz(chd.mdim) = nnz(~ismember(txi, ppi));
chd.data = zeros(sz, 'like', chd.data);
if tx_multi, chd.t0 = chd.t0(1:2:end); end % every other is true tx
chd.t0 = chd.t0(~ismember(1:numel(chd.t0), ppi)); % skip pilot pulses

% pass beamforming params for bfQUPS
global QUPS_BF_PARAMS; 
QUPS_BF_PARAMS.ppi = [];
QUPS_BF_PARAMS.rx_multi = us.xdc.numel > 128; % more elems than channels
QUPS_BF_PARAMS.tx_multi = tx_multi; % more active tx elems than channels
QUPS_BF_PARAMS.rcvbuf = 1; % matching Receive.bufnum
QUPS_BF_PARAMS.evi = evi;
QUPS_BF_PARAMS.bfi = bfi;
QUPS_BF_PARAMS.D = chd.getPassbandFilter(us.xdc.bw);
QUPS_BF_PARAMS.fig = 1; % figure number

% Save verasonics pre-process variables 
filename = char(fullfile("MatFiles","qups-vsx.mat")); 
save(filename, '-struct', 'vs');

% get a copy of this file
setup_file = mfilename('fullpath') + ".m";
if ~isempty(setup_file) && ~startsWith(setup_file, tempdir)
    code = readlines(setup_file); 
else 
    code = string.empty; 
end

% save
conf_file = fullfile("MatFiles","qups-conf.mat"); % configuration
save(conf_file, "us", "chd", "QUPS_BF_PARAMS", "code");
%save(conf_file, "us", "chd");

% Save QUPS configuration variables
%save(fullfile("MatFiles","qups-conf.mat"), "us", "chd");

% set output save directory
global VSXOOD_SAVE_DIR;
VSXOOD_SAVE_DIR = fullfile(pwd, 'tmp');
if ~exist(VSXOOD_SAVE_DIR, 'dir'), mkdir(VSXOOD_SAVE_DIR); end

% clear external functions
clear RFDataImg RFDataProc RFDataStore RFDataCImage imagingProc cEstFSA_RT;

VSX;

%% Save verasonics post-process variables 
% This is relevant for future channel data reconstruction by using the
% verasonics structs
global VSXOOD_SAVE_DIR; % because VSX clears the workspace

vs = update_vstruct();
save(fullfile(VSXOOD_SAVE_DIR, 'qups-vsx-post.mat'), '-struct', 'vs');
