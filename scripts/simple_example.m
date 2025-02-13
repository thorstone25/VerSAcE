%% Construct the system description
c0 = 1500;
[xdc, Trans] = TransducerVerasonics("P4-2v");
seq = SequenceRadial('type', 'PW', 'angles', -25 : 2.5 : 25, 'c0', c0);
scn = ScanCartesian('x', 80e-3 * [-1 1]./2, 'z', [0, 100e-3]);
us = UltrasoundSystem('xdc', xdc, 'seq', seq, 'scan', scn);
[scn.dx, scn.dz] = deal(us.lambda / 4);

%% Create necessary VSX objects (outside of QUPS)
vRes = VSXResource(); % global resource definition
vTW  = VSXTW('type','parametric', 'Parameters', [Trans.frequency, 0.67, 2, 1]); % tx waveform
vTPC = VSXTPC('name','Default', 'hv', Trans.maxHighVoltage); % max power

% Create a VSXBlock
vb = QUPS2VSX(us, Trans, vRes, 'vTW', vTW, 'vTPC', vTPC, 'recon_VSX', true, 'saver_custom',true);
global VERSACE_PARAMS;
VERSACE_PARAMS.save_dir = fullfile(pwd, 'data', 'test')

% Create the structs
vs = link(vb, vRes, Trans, 'TXPD', true); % link
pt1; vs.Media = Media; % add simulation media

% Save
filename = fullfile(vantageroot, 'MatFiles', 'qups-vsx.mat');
save(filename, '-struct', 'vs');

%% Run
run VSX;

%% Post-process
if ~exist('RcvData', 'var'), RcvData = {RData}; end
c0   = Resource.Parameters.speedOfSound;
[us, chd] = UltrasoundSystem.Verasonics(Trans, TX, TW, 'c0', c0, 'Receive', Receive, 'RcvData', RcvData, 'PData', PData); % import to QUPS                                [us, chd] = UltrasoundSystem.Verasonics(Trans, TX, TW, 'Receive', Receive, 'RcvData', RcvData); % import to QUPS
sct = Scatterers.Verasonics(Media, 'c0', c0, 'scale', c0./us.xdc.fc);
[v, ~, wv] = Waveform.Verasonics(TW, us.xdc.fc); % excitation, 2wy waveform

% plot the import
chd = hilbert(singleT(chd));
figure; plot(us); hold on; plot(sct, '.');
figure; imagesc(chd); dbr echo 80; animate(chd.data, 'loop', false);

