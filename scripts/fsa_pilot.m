%% Construct a system
% transducer
c0 = 1500; % sound speed 
Trans = computeTrans(struct("name", char("L7-4"), 'units', 'mm')); % transducer
xdc = Transducer.Verasonics(Trans); % VSX -> QUPS

% pulse sequences
P = 16; % pilot every P pulses
seq0 = Sequence('type','FSA','c0', c0, 'numPulse',xdc.numel); % FSA acquisition
seqp = Sequence('type','FC' ,'c0', c0, 'focus',[0;0;50e-3]*ones([1,1+xdc.numel/P])); % pilot pulses

% systems
us0 = UltrasoundSystem('xdc', xdc, 'seq', seq0);
usp = UltrasoundSystem('xdc', xdc, 'seq', seqp);

%% create VSX objects
% constant resources
vres = VSXResource(); % system-wide resource

% excitation pulse
vTW = VSXTW('type','parametric', 'Parameters', [Trans.frequency, 0.67, 1, 1]); % tx waveform

% make VSX blocks: reference and pilot
F = 10; % number of frames
[vb0, chd0] = QUPS2VSX(us0, Trans, vres, "frames", F, 'vTW', vTW, 'recon_VSX',true); % ref
[vbp, chdp] = QUPS2VSX(usp, Trans, vres, "frames", F, 'vTW', vTW); % pilot

% find transmit events before which we will insert the pilot pulses
[vbp.capture.info] = dealfun(@(x)"Pilot " + x, vbp.capture.info); % prepend name for clarity
istr   = "Tx "+(0:P:xdc.numel)'+" - Ap 1 - Frame "+1;  % matching info strings
[~, i] = ismember(istr, [vb0.capture.info]); % find matching (previous) transmits
es     = [VSXEvent; vb0.capture(i(2:end,1))]; % insert pulse after these events (dummy event to cause matching index to be 0 for frame 1)

% merge the blocks to insert the pilot pulses
vb  = copy(vb0); % make the final block the modified fsa block
for j = 1:numel(es) % for each remaining pilot pulse
    [~,k] = ismember(es(j), vb.capture); % find matching event (row)
    vb.capture = cat(2, vb.capture(:,1:k,:), vbp.capture(:,j,:), vb.capture(:,k+1:end,:)); % insert pilot pulse
end
vb.next = vb.capture(1); % return to beginning of the block at end

%% convert to VSX structures
vs = link(vb, vres); % link
vs.Trans = Trans; % add Trans
pt1; vs.Media = Media; % add simulation media
vs.Media.function = 'movePoints'; % make points move

% force in simulation mode for testing
vs.Resource.Parameters.simulateMode = 1; % 1 to force simulate mode, 0 for hardware

% save 
filename = char(fullfile("MatFiles","qups-vsx.mat")); 
save(filename, '-struct', 'vs');

% save
[us, chd] = deal(us0, chd0);
save(fullfile("MatFiles","qups-conf.mat"), "us", "chd");

% clear external functions (reset)
clear RFDataImg RFDataProc RFDataStore imagingProc;

%% Launch
% run VSX;
