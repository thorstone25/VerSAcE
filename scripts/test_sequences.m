 %#ok<*UNRCH> unreachable code due to const coded flags
 %#ok<*GVMIS> using global variables
 
 %% Create different system configurations
c0 = 1500; % sound speed
xdc_name = "L12-3v"; % choose transducer
% xdc_name = "L12-5 50mm"; % choose transducer
% xdc_name = "P4-2v";
Trans = computeTrans(struct("name", char(xdc_name), 'units', 'mm')); % transducer
xdc = Transducer.Verasonics(Trans);
uss = UltrasoundSystem('xdc', xdc); 
uss.seq.c0 = c0;
uss.scan.zb(2) = 150 * uss.lambda;
uss.scan.xb = 1e-3 * [-20 20];
[uss.scan.dx, uss.scan.dz] = deal(uss.lambda / 2);

apd_center = false; % use only ceneter aperture

% sequences
apd   = Sequence.apWalking(xdc.numel,xdc.numel/4,2); % active aperture
pf    = xdc.focActive(apd, 50e-3);
seqfc  = Sequence(      'type', 'FC' , 'c0', c0, 'focus', pf, 'apd', apd);
seqfsa = Sequence(      'type', 'FSA', 'c0', c0, 'numPulse', xdc.numel);
seqpw  = SequenceRadial('type', 'PW' , 'c0', c0, 'angles', -25 : 0.5 : 25);

% spatial downsampling
Ds = 4; 
seqfsadv = Sequence(  'type', 'FSA', 'c0', c0, 'numPulse', xdc.numel / Ds);
seqfsadv.apd = zeros([xdc.numel, seqfsadv.numPulse]);
seqfsadv.apd(1:Ds:end,:) = eye(xdc.numel / Ds);

uss = copy(repmat(uss, [1,4]));
[uss.seq] = deal(seqfsa, seqpw, seqfc, seqfsadv);

% transmit-aperture image apodization schemes
apod0 = cell([1, numel(uss)]);
for i = numel(uss):-1:1
    switch i
      case 1, apod0{i} = uss(i).apAcceptanceAngle(45);
      case 2, apod0{i} = swapdim(apTxParallelogram(uss(i),uss(i).seq.angles, [-15 15]),4,5);
      case 3, apod0{i} = swapdim(uss(i).apMultiline(),4,5);
      case 4, apod0{i} = sub(uss(i).apAcceptanceAngle(45),1:Ds:uss(i).xdc.numel,4);
    end
end

%%
% selection
seq_ind = [3 1]; % sequence index
[us, apod] = deal(copy(uss(seq_ind)), apod0(seq_ind)); % choose pulse sequence template
if apd_center && false
    N = min(128, us.xdc.numel); % active aperture size
    apd = circshift([ones(1,N), zeros(1,us.xdc.numel-N)], (us.xdc.numel-N)/2)';
    us.seq.apodization_ = repmat(apd, [1 us.seq.numPulse]);
end

% constant resources
vres = VSXResource(); % system-wide resource
vTW  = VSXTW('type','parametric', 'Parameters', [Trans.frequency, 0.67, 1, 1]); % tx waveform
vTPC = VSXTPC('name','Default', 'hv', Trans.maxHighVoltage); % max power

% make blocks
for i = 1:numel(us)
    [vb(i), chd(i)] = QUPS2VSX(us(i), Trans, vres ...
        ... ,"apod", apod{i} ... TODO: adapt for tx-apodization
        ,"frames", 4 ...
        ,'vTW', vTW, 'vTPC', vTPC ...
        ,'recon_VSX', i == 1 ... imaging
        ,'saver_custom', i == 2 ... data collection
        ,'set_foci', true ...
        ,'range', [0 us(i).scan.zb] ... in m
        ,'range', [0 50]*1e-3 ... in m
        ); %#ok<SAGROW> % make VSX block

    % tag the block index
    for j = 1:numel(vb(i).capture)
        vb(i).capture(j).info = join([vb(i).capture(j).info,"Blk "+i], " - ");
    end
end

% connect multiple blocks (set the 'next' prop to first 'capture')
[vb(1:end-1).next] = dealfun(@(x)x(1), vb(2:end).capture); % each to next
vb(end).next = vb(1).capture(1); % last to first

% set peak cut-off for VSX imaging (heuristic)
for i = 1:numel(us)
    txs = [cat(2,cat(2,vb(i).capture).tx)];
    txpow = 2 * (sum(us(i).seq.apodization(us(i).xdc),1) ./ us(i).xdc.numel); % total weight per tx
    txpow = max(txpow); % use a single max power value
    [txs.peakCutOff] = deal(txpow);
end

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
filename = char(fullfile(vantageroot, "MatFiles","qups-vsx.mat")); 
save(filename, '-struct', 'vs');

% save
save(fullfile(vantageroot, "MatFiles","qups-conf.mat"), "us", "chd");

% set output save directory
global VERSACE_PARAMS;
VERSACE_PARAMS.save_dir = fullfile(pwd, 'tmp');
if ~exist(VERSACE_PARAMS.save_dir, 'dir'), mkdir(VERSACE_PARAMS.save_dir); end

% clear external functions
clear RFDataImg RFDataProc RFDataStore;

%% Run
run VSX;

%% Post processing - parse and beamform data from the 2nd block
blk = 2; % vsxblock index

% grab most recent dataset
global VERSACE_PARAMS; %#ok<REDEFGG> % cleared by VSX
fls = dir(fullfile(VERSACE_PARAMS.save_dir, '*_*_*.mat')); % access mat-files in or folder
dates = reshape(datetime([fls.datenum], 'ConvertFrom', 'datenum'), size(fls)); % file dates
i = argmax(dates); % most recent file
if isempty(i), disp("No files found"); return; end
vs = load(fullfile(fls(i).folder, fls(i).name)); % data

% extract block 2 (the data block)
evinf = string({vs.Event.info}); % event descriptions
evi = contains(evinf, "Tx") & contains(evinf, "Blk "+blk) & ~contains(evinf, "Jump"); % filter
[evir, evit] = deal(evi, evi & contains(evinf, "Frame 1") & contains(evinf, "Ap 1")); % corresponding receive and transmit events
rxi = [vs.Event(evir).rcv]; assert(all(rxi > 0), "Detected a non-receive event!" ); % rx indices
txi = [vs.Event(evit).tx ]; assert(all(txi > 0), "Detected a non-transmit event!"); % tx indices
txs = vs.TX(txi); % TX
rxs = vs.Receive(rxi); % Receive
tw  = vs.TW(unique([txs.waveform])); % TW

% convert to QUPS
c0 = vs.Resource.Parameters.speedOfSound;
[us, chd] = UltrasoundSystem.Verasonics(vs.Trans, txs, tw, ...
    'c0', c0, 'Receive', rxs, 'RcvData', {vs.RData}, "PData", vs.PData ...
    );
scat = Scatterers.Verasonics(vs.Media,'scale',us.lambda,'c0',c0);

% de-multiplexing
Mx = 2; % 1 if none, 2 if Tx or Rx, 4 if both Tx and Rx
x = 0; for i = 1:Mx, x = x + sub(chd.data,i:Mx:chd.M,chd.mdim); end % sum over physical transmits
chd.data = x;
chd = join(chd, chd.mdim); % HACK: make t0 scalar over tx if possible

% beamform (sanity check)
chd = hilbert(singleT(chd));
D = chd.getPassbandFilter(us.xdc.bw);
b = DAS(us, filter(chd,D));

%% display
figure("Name", "Post Processing Sanity Check"); tiledlayout('flow');
nexttile(); imagesc(us.scan, b); dbr b-mode 60; 
hold on; plot(us.xdc); plot(scat, 'r.');
ttls = "Tx "+(1:chd.M)'+" : Frm "+(1:size(chd.data,4));
nexttile(); h = imagesc(chd); 
animate(chd.data, h, 'title', ttls(:)', 'loop',false);




