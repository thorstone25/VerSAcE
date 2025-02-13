%% Originally modified from test_sequences.m
% Description: Simple FSA sequence for Verasonics L12-3V
% Saves everything you would need in order to post-beamform RF data
% You may use the outputs of this script with bfRData_detailed.m
% to produce a beamformed image.
% To run, you will need to have QUPS, VSX-OOD, and Vantage added to path!
% OPtional: see the tutorial in the README.md file to setup a matlab proj.
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

% spatial downsampling
Ds = 4; 
seqfsadv = Sequence(  'type', 'FSA', 'c0', c0, 'numPulse', xdc.numel / Ds);
seqfsadv.apd = zeros([xdc.numel, seqfsadv.numPulse]);
seqfsadv.apd(1:Ds:end,:) = eye(xdc.numel / Ds);

uss = copy(repmat(uss, [1,1]));
uss.seq = seqfsa;

% apodization schemes
apod0 = cell([1, numel(uss)]);
apod0{1} = uss(1).apAcceptanceAngle(45);

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

% set peak cut-off for VSX imaging
vec = @(x)x(:);
evs = arrayfun(@(x) {vec(x.capture(:))'}, vb);
evs = [evs{:}];
txs = [evs.tx];
[txs.peakCutOff] = deal(2 / 128);

%% convert to VSX structures
global vs; % make global to access in within function
vs = link(vb, vres, Trans, 'TXPD', true); % link
pt1; vs.Media = Media; % add simulation media

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
filename = char(fullfile(vantageroot, "MatFiles","qups-vsx.mat")); 
save(filename, '-struct', 'vs');

% get a copy of this file
setup_file = mfilename('fullpath') + ".m";
if ~isempty(setup_file) && ~startsWith(setup_file, tempdir)
    code = readlines(setup_file); 
else 
    code = string.empty; 
end

% save
conf_file = fullfile(vantageroot, "MatFiles","qups-conf.mat"); % configuration
save(conf_file, "us", "chd", "QUPS_BF_PARAMS", "code");

% set output save directory
global VERSACE_PARAMS;
VERSACE_PARAMS.save_dir = fullfile(pwd, 'tmp');
if ~exist(VERSACE_PARAMS.save_dir, 'dir'), mkdir(VERSACE_PARAMS.save_dir); end

% clear external functions
clear RFDataImg RFDataProc RFDataStore RFDataCImage imagingProc cEstFSA_RT;

run VSX;

%% Save verasonics post-process variables 
% This is relevant for future channel data reconstruction by using the
% verasonics structs
global VERSACE_PARAMS; % because VSX clears the workspace

vs = update_vstruct();
save(fullfile(VERSACE_PARAMS.save_dir, 'qups-vsx-post.mat'), '-struct', 'vs');
