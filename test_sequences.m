%% Create different system configurations
c0 = 1500; % sound speed
xdc = TransducerArray.L11_5v();
uss = UltrasoundSystem('xdc', xdc); 
uss.seq.c0 = c0;
[uss.scan.dx, uss.scan.dz] = deal(uss.lambda / 2);


% sequences
pf = [0;0;50e-3] + [1e-3;0;0] .* (-20 : 5 : 20);
seqfsa = Sequence(     'type', 'FSA', 'c0', c0, 'numPulse', xdc.numel);
seqpw = SequenceRadial('type', 'PW' , 'c0', c0, 'angles', -25 : 0.5 : 25); 
seqfc = Sequence(      'type', 'VS' , 'c0', c0, 'focus', pf);
seqfsadv2 = Sequence(     'type', 'FSA', 'c0', c0, 'numPulse', xdc.numel / 2);
seqfsadv2.apodization_ = zeros([xdc.numel, xdc.numel / 2]);
seqfsadv2.apodization_(1:2:end,:) = eye(xdc.numel / 2);

uss = copy(repmat(uss, [1,4]));
[uss.seq] = deal(seqfsa, seqpw, seqfc, seqfsadv2);

%%
% selection
seq_ind = 4;
us = copy(uss(seq_ind)); % choose pulse sequence template
xdc_name = "L7-4"; % choose transducer

% create VSX objects
% constant resources
vres = VSXResource(); % system-wide resource
Trans = computeTrans(struct("name", char(xdc_name), 'units', 'mm')); % transducer
vTW = VSXTW('type','parametric', 'Parameters', [Trans.frequency, 0.67, 1, 1]); % tx waveform

% make blocks
[vb, chd] = QUPS2VSX(us, Trans, vres, "frames", 1, 'vTW', vTW, ...
    'recon_VSX', false, 'recon_custom', true, ...
    'custom_imaging', true, ... % sound speed imaging
'recon_custom_delays', false, 'saver_custom', false ...
    ); % make VSX block
% [vb(1), chd(1)] = QUPS2VSX(uss(1), Trans, vres, "frames", 1, 'vTW', vTW); % make VSX block
% [vb(2), chd(2)] = QUPS2VSX(uss(2), Trans, vres, "frames", 4, 'vTW', vTW); % make VSX block
% [vb.next] = deal(vb(2).capture(1), vb(1).capture(1)); % start at beginning of alternate sequence

% DEBUG: test the manual receive delays
%{
for i = 1:numel(seq_ind)
    [~, tau_rx, tau_tx] = bfDAS(uss(seq_ind(i)), chd(i), 'delay_only', true);
    vRecon = unique([vb(i).capture.recon]); % find Recon (exactly 1 exists)
    if ~isscalar(vRecon), continue; end
    setVSXLUT(vRecon, tau_rx, tau_tx - swapdim(chd(i).t0,chd(i).mdim,5), uss(seq_ind(i)).xdc.fc);% broken for Vantage 4.3
end
%}

% convert to VSX structures
vs = link(vb, vres); % link
vs.Trans = Trans; % add Trans
pt1; vs.Media = Media; % add simulation media

% DEBUG: test the manual receive delays
% [~, tau_rx, tau_tx] = bfDAS(us, chd, 'delay_only', true);
% [vs.Recon, vs.ReconInfo] = setVSXLUT(vs.Recon, vs.ReconInfo, vs.PData, tau_rx, tau_tx + swapdim(chd.t0,chd.mdim,5), us.xdc.fc);

% force in simulation mode for testing
vs.Resource.Parameters.simulateMode = 0; % 1 to force simulate mode, 0 for hardware

% save 
filename = char(fullfile("MatFiles","qups-vsx.mat")); 
save(filename, '-struct', 'vs');

% fix Transducer to match VSX
us.xdc = Transducer.Verasonics(Trans);

% save
save(fullfile("MatFiles","qups-conf.mat"), "us", "chd");

% clear external functions
clear RFDataImg RFDataProc RFDataStore RFDataCImage;

% VSX;

%% 
