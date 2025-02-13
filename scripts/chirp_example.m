%% Construct the system description
c0 = 1500; 
[xdc, Trans] = TransducerVerasonics("P4-2v");
seq = SequenceRadial('type', 'PW', 'angles', -25 : 2.5 : 25, 'c0', c0);
scn = ScanCartesian('x', 80e-3 * [-1 1]./2, 'z', [0, 100e-3]);
us = UltrasoundSystem('xdc', xdc, 'seq', seq, 'scan', scn);
[scn.dx, scn.dz] = deal(us.lambda / 4);

%% Create necessary VSX objects (outside of QUPS)
hashpt = true; % set true if system has HIFU or Extended Transmit support
vRes = VSXResource(); % global resource definition
vTPC = VSXTPC('name','High', 'hv', Trans.maxHighVoltage, 'ishpt', hashpt); % max power

% Define a stepped LFM chirp
df = 1.0; % step size in MHz
cyc_per_frq = 2; % cycles per frequency
frqs = 250./2./flip(6:197); % all supported tx frequencies
frqs = repmat(frqs(xdc.bw(1)/1e6 < frqs & frqs < xdc.bw(end)/1e6), [cyc_per_frq, 1]); % within bandwidth
frqs = frqs(:,arrayfun(@(f)argmin(abs(frqs(1,:) - f)), xdc.bw(1)/1e6 : df : xdc.bw(end)/1e6)); % 1MHz step size
vTW = VSXTW('type','envelope', 'envFrequency', frqs(:), 'envPulseWidth', repmat(0.67, [1 numel(frqs)]), 'envNumCycles', numel(frqs)); % make waveform

% Create a VSXBlock
vb = QUPS2VSX(us, Trans, vRes, 'vTW', vTW, 'vTPC', vTPC, 'recon_VSX', true);

% Create the structs
vs = link(vb, vRes, Trans, 'TXPD', true); % link
pt1; vs.Media = Media; % add simulation media

% Save
filename = fullfile(vantageroot, 'MatFiles','qups-vsx.mat');
save(filename, '-struct', 'vs');

%% Run
run VSX;

%% Post-process
% import to QUPS
if ~exist('RcvData', 'var'), RcvData = {RData}; end % alias
c0   = Resource.Parameters.speedOfSound;
[us, chd] = UltrasoundSystem.Verasonics(Trans, TX, TW, 'c0', c0, 'Receive', Receive, 'RcvData', RcvData, 'PData', PData); % import to QUPS
sct = Scatterers.Verasonics(Media, 'c0', c0, 'scale', c0./us.xdc.fc);
[v, ~, wv] = Waveform.Verasonics(TW, us.xdc.fc); % voltage, 2wy waveform

% pre-processing
chd = hilbert(singleT(chd));

% pulse compression
wv = Waveform('t', wv.time, 'samples', wv.samples .* hamming(numel(wv.samples))); % apply windowing
chdp = convt(chd, reverse(wv)); % matched filter
chdp.t0 = chdp.t0 - 2*Trans.lensCorrection / us.xdc.fc; % residual phase adjustment?

% image
b = DAS(us, chdp, us.apAcceptanceAngle(30), us.apTxParallelogram());

% plot
figure; plot(us); hold on; plot(sct, '.'); % geometry
figure; imagesc(chd); dbr echo 80; animate(chd.data, 'loop', false); % Channel Data
figure; imagesc(us.scan, b); dbr b-mode 60; hold on; plot(sct, 'r.'); % B-mode
