%% Construct a system
% transducer
% xdc_name = "L12-5 50mm"; % choose transducer
xdc_name = "L12-3v"; % choose transducer
c0 = 1500; % sound speed 
Trans = computeTrans(struct("name", char(xdc_name), 'units', 'mm')); % transducer
xdc = Transducer.Verasonics(Trans); % VSX -> QUPS

% imaging region
scan = ScanCartesian('x', 1e-3*[-20 20], 'z', 1e-3*[0 50], 'y', 0);

% pulse sequences
seq0 = Sequence('type','FSA','c0', c0, 'numPulse', xdc.numel); % FSA acquisition
% seq0 = SequenceRadial('type','PW','c0', c0, 'angles',-23 : 2 : 23); % PW acquisition
P = seq0.numPulse/8; % pilot every P pulses
seqp = Sequence('type','VS' ,'c0', c0, 'focus', [0 0 50e-3]' * ones([1,1+seq0.numPulse/P])); % pilot pulses
M = max(0, xdc.numel - 128) / 2; % buffer
seqp.apodization_ = repmat([zeros(M,1);ones(xdc.numel-2*M,1);zeros(M,1)],[1,seqp.numPulse]); % restrict pilot pulses to within 128 elements to avoid tx multiplexing

% systems
us0 = UltrasoundSystem('xdc', xdc, 'seq', seq0, 'scan', scan);
usp = UltrasoundSystem('xdc', xdc, 'seq', seqp, 'scan', scan);
[scan.dx, scan.dz] = deal(us0.lambda / 2);

% combine FSA and pilot sequences
txapd0 = reshape(seq0.apodization(xdc), xdc.numel, P, []);
txapdp = reshape(seqp.apodization(xdc), xdc.numel, 1, []);
txapd  = cat(2, txapd0, txapdp(:,:,2:end));
txapd  = [txapdp(:,:,1), txapd(:,:)];

txtau0 = reshape(seq0.delays(xdc), xdc.numel, P, []);
txtaup = reshape(seqp.delays(xdc), xdc.numel, 1, []);
txtau  = cat(2, txtau0, txtaup(:,:,2:end));
txtau  = [txtaup(:,:,1), txtau(:,:)];

% pilot pulse indices
M = size(txtau,2); % # physical transmits (pulses)
ppi = 1 : (1+P) : M;

%% create VSX objects
% number of frames
F = 1;

% constant resources
vres = VSXResource(); % system-wide resource
vres.Parameters.simulateMode = 1; % 1 to force simulate mode, 0 for hardware

% excitation pulse
vTW = VSXTW('type','parametric', 'Parameters', [Trans.frequency, 0.67, 1, 1]); % tx waveform
vTPC = VSXTPC('name','Default', 'hv', Trans.maxHighVoltage/2); % half max power

us = copy(us0);
us.seq = SequenceGeneric('apod', txapd, 'del', txtau, 'c0', seq0.c0, 'numPulse', size(txapd,2));

% make VSX blocks: reference and pilot
[vb, chd] = QUPS2VSX(us, Trans, vres ...
    ,'range', [0 50]*1e-3  ...
    ... ,'range', [0 us.scan.zb(2)] ...
    ,'vTW', vTW, 'vTPC', vTPC ...
    ,"frames", F,'set_foci', false ...
    ,'recon_VSX', false ...
    ,'recon_custom', false, 'saver_custom', true ...
    ); % ref
rxbuf = unique([cat(2,vb.capture.rcv).bufnum]);

% add QUPS processing
%{
[ui, ev] = addReconFun("bfQUPS", rxbuf, us.scan, vres, 'UItyp', 'VsToggleButton', 'UIpos', "UserB4");
[vb.vUI(end+(1:numel(ui))), vb.post(end+(1:numel(ev)))] = deal(ui, ev);
iopts = struct(ev(2).process.Parameters{:}); % imaging process params
% iopts.displayWindow.Colormap; % display window colormap
iopts.reject = 0; % percent of maxV/4 that is rejected
iopts.persistMethod = 'none';
iopts.pgain = 4; % sqrt(seq0.numPulse);
ev(2).process.Parameters = namedargs2cell(iopts); % set modified valus
%}

% add sound speed estimation
[ui, ev] = addDataFun("ceQUPS", rxbuf, 'UItyp', 'VsToggleButton', 'UIpos', "UserB5");
[vb.vUI(end+(1:numel(ui))), vb.post(end+(1:numel(ev)))] = deal(ui, ev);


% remove recon for pilot pulses
recon = unique([vb.capture.recon]); % Recon
if ~isempty(recon)
bound_mode = {recon.RINums([1 end]).mode}; % begin/end mode
recon.RINums(ppi) = []; % delete pilot pulses from reconinfo
[recon.RINums([1 end]).mode] = bound_mode{:}; % replace
end

%% convert to VSX structures
vs = link(vb, vres, Trans, 'TXPD', false); % link
pt1; vs.Media = Media; % add simulation media
% vs.Media.function = 'movePoints'; % make points move

%% Callback function pre-processing
% whether tx multiplexed
tx_multi = max(sum(logical(us.seq.apodization(us.xdc)),1)) > 128;

% pre-processing indexing
evi = find(startsWith({vs.Event.info}, "Tx ")); % events with transmits
txi = double(string(extractBetween({vs.Event(evi).info}, "Tx ", " - Ap"))); % physical transmit indices
if tx_multi, txi = ceil(txi/2); end % to match duplication
bfi = find(~ismember(txi, ppi)); % beamforming event indices

% remove pilot pulses from Sequence/ChannelData definitions
[usv, chdv, us, chd] = deal(us, chd, copy(us), copy(chd)); % store old, save new
us.seq = seq0; % final sequence

sz = size(chd.data);
sz(chd.mdim) = nnz(~ismember(txi, ppi));
chd.data = zeros(sz, 'like', chd.data);
if tx_multi, chd.t0 = chd.t0(1:2:end); end % every other is true tx
chd.t0 = chd.t0(~ismember(1:numel(chd.t0), ppi)); % skip pilot pulses

% pass beamforming params for bfQUPS
global QUPS_BF_PARAMS; 
QUPS_BF_PARAMS.ppi = ppi;
QUPS_BF_PARAMS.rx_multi = us.xdc.numel > 128; % more elems than channels
QUPS_BF_PARAMS.tx_multi = tx_multi; % more active tx elems than channels
QUPS_BF_PARAMS.rcvbuf = 1; % matching Receive.bufnum
QUPS_BF_PARAMS.evi = evi;
QUPS_BF_PARAMS.bfi = bfi;
QUPS_BF_PARAMS.D = chd.getPassbandFilter(us.xdc.bw);
QUPS_BF_PARAMS.fig = 1; % figure number


%% save 
conf_file = fullfile("MatFiles","qups-conf.mat"); % configuration
filename = char(fullfile("MatFiles","qups-vsx.mat")); % Vantage
save(filename, '-struct', 'vs');

% get a copy of this file
setup_file = mfilename('fullpath');
if ~isempty(setup_file) && ~startsWith(setup_file, tempdir)
    code = readlines(setup_file); 
else 
    code = string.empty; 
end

% save
save(conf_file, "us", "chd", "QUPS_BF_PARAMS", "code");

% set save directory for data store (TODO: rename and document)
global VSXOOD_SAVE_DIR;
VSXOOD_SAVE_DIR = fullfile(pwd, "data/pilot-pulse-test",string(datetime('now'),'yyMMddHHmmSS')); % make a path relative to the current location
if ~exist(VSXOOD_SAVE_DIR, 'dir'), mkdir(VSXOOD_SAVE_DIR); end
copyfile(conf_file, VSXOOD_SAVE_DIR); % save a copy of og files too
copyfile(filename , VSXOOD_SAVE_DIR);

% clear external functions
clear bfQUPS ceQUPS;% cEstFSA_RT;

% VSX;

%% 
