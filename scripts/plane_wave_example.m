 %#ok<*UNRCH> unreachable code due to const coded flags
 %#ok<*GVMIS> using global variables
 
 %% Create a system configuration
xdc_name = "L12-3v"; % choose transducer
% xdc_name = "L12-5 50mm"; % choose transducer
% xdc_name = "P4-2v";
c0 = 1500; % sound speed
Trans = computeTrans(struct("name", char(xdc_name), 'units', 'mm')); % transducer
xdc = Transducer.Verasonics(Trans); % VSX -> QUPS
seq = SequenceRadial('type', 'PW' , 'c0', c0, 'angles', -8 : 0.5 : 8); % plane waves
seq.apd = [zeros(1,32), ones(1,128), zeros(1,32)]' .* ones(1,seq.numPulse); % restrict to center channels
scan = ScanCartesian('x', -20e-3 : 10*20e-6 : 20e-3, 'z', 0 : 20e-6 : 40e-3); % image domain
us  = UltrasoundSystem('seq', seq, 'xdc', xdc, 'scan', scan, 'fs', 4*xdc.fc);
% apod = swapdim(apTxParallelogram(us, seq.angles, [-15 15]),4,5);

%%
fc = 1e6*Trans.frequency;

% constant resources
vres = VSXResource(); % system-wide resource
vTW  = VSXTW('type','parametric', 'Parameters', [fc/1e6, 0.67, 1, 1]); % tx waveform
vTPC = VSXTPC('name','Default', 'hv', Trans.maxHighVoltage); % max power
vTGC = VSXTGC( ...
        'CntrlPts', [0,297,424,515,627,764,871,1000],...
        'rangeMax', 50*1e-3 ./ us.lambda ...
        );

% make a block
[vb, chd] = QUPS2VSX(us, Trans, vres ...
    ... ,"apod", apod ...
    ,"frames", 2 ...
    ,"custom_fs", 25e6 ...
    ,'vTW', vTW, 'vTPC', vTPC, 'vTGC', vTGC ...
    ,'recon_VSX', true ... imaging
    ,'saver_custom', true ... data collection
    ,'set_foci', true ...
    ,'range', [0 50]*1e-3 ... in m
    ); % make VSX block

% TODO: quit when done
% vb.next = VSXEvent("info", "quit", "seqControl", VSXSeqControl("command","returnToMatlab"));

% set peak cut-off for VSX imaging (heuristic)
% TODO: move to linking? Move to addVSXRecon?
for i = 1:numel(us)
    txs = [cat(2,cat(2,vb.capture).tx)];
    txpow = 2 * (sum(us.seq.apodization(us.xdc),1) ./ us.xdc.numel); % total weight per tx
    txpow = max(txpow); % use a single max power value
    [txs.peakCutOff] = deal(txpow);
end

%% convert to VSX structures
global vs; % make global to access in within function
vs = link(vb, vres, Trans, 'TXPD', false); % link
pt1; vs.Media = Media; % add simulation media

% force in simulation mode for testing
vs.Resource.Parameters.simulateMode = 1; % 1 to force simulate mode, 0 for hardware

%% save 
filename = char(fullfile(vantageroot, "MatFiles","qups-vsx.mat")); 
save(filename, '-struct', 'vs');
save(replace(filename,"qups-vsx.mat","qups-conf.mat"), "us", "chd");

% set output save directory
global VERSACE_PARAMS;
svnm = "f"+round(fc)+"-"+string(datetime("now"), "yyMMdd-HHmmss");
VERSACE_PARAMS.save_dir = fullfile(pwd, "data", "freq-est", svnm);
if ~exist(VERSACE_PARAMS.save_dir, 'dir'), mkdir(VERSACE_PARAMS.save_dir); end

% clear external functions
clear RFDataImg RFDataProc RFDataStore;

%% Run
VSX;

%% Post-processing - parse and save the raw data
global VERSACE_PARAMS hf; %#ok<REDEFGG> % cleared by VSX
fl = dir(fullfile(VERSACE_PARAMS.save_dir, "*.mat")); % get the files
fl = string(fullfile({fl.folder}, {fl.name})); % all mat-files
if isempty(fl), disp("No files found."); return; % nothing saved
else, fl = fl(1); % choose first
end

% load
load(fl); % everything

% parse
c0 = Resource.Parameters.speedOfSound;
[us, chd0] = UltrasoundSystem.Verasonics(Trans, TX, TW, 'c0', c0, 'PData', PData, 'Receive', Receive, 'RcvData', {RData});

% pre-processing
chd0.data = chd0.data(:,1:2:end,:,:) + chd0.data(:,2:2:end,:,:);
chd = filter(hilbert(singleT(chd0)), chd0.getPassbandFilter(us.xdc.bw));

% image
[us.scan.dx, us.scan.dz] = deal(us.lambda / 4);
apod = us.apApertureGrowth(1.5);
b = DAS(us, chd, 'apod', apod);

% display
hf = figure(1); 
imagesc(us.scan, b); 
dbr b-mode 60;

% convert to native structs
us   = obj2struct(us );
chd  = obj2struct(chd);
chd0 = obj2struct(chd0);
b    = gather(b);

% save
save(fl, '-append', 'us', 'chd', 'chd0', 'b', 'apod');
saveas( hf, replace(fl,".mat",".png"));
savefig(hf, replace(fl,".mat",".fig"), 'compact');

return;
%% Post processing for multiple blocks: parse and beamform data from the 2nd block
blk = 2; % vsxblock index

% grab most recent dataset
global VERSACE_PARAMS; %#ok<REDEFGG> % cleared by VSX
flds = dir(fullfile(VERSACE_PARAMS.save_dir, '*_*_*.mat')); % access mat-files in or folder
dates = reshape(datetime([flds.datenum], 'ConvertFrom', 'datenum'), size(flds)); % file dates
i = argmax(dates); % most recent file
vs = load(fullfile(flds(i).folder, flds(i).name)); % data

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
nexttile(); h = imagesc(chd); animate(chd.data, h, title=ttls(:)');




