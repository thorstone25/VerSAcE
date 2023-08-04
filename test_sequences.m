%% Create different system configurations
c0 = 1500; % sound speed
xdc = TransducerArray.L11_5v();
uss = UltrasoundSystem('xdc', xdc); 
uss.seq.c0 = c0;
[uss.scan.dx, uss.scan.dz] = deal(uss.lambda / 2);


% sequences
pf = [0;0;50e-3] + [1e-3;0;0] .* (-20 : 5 : 20);
seqfsa = Sequence(     'type', 'FSA', 'c0', c0, 'numPulse', xdc.numel);
seqpw = SequenceRadial('type', 'PW' , 'c0', c0, 'angles', -25 : 5 : 25); 
seqfc = Sequence(      'type', 'VS' , 'c0', c0, 'focus', pf);

uss = copy(repmat(uss, [1,3]));
[uss.seq] = deal(seqfsa, seqpw, seqfc);

%%
vres = VSXResource();
xdc_name = "L11-5v";
us = copy(uss(2)); % just 1 selection for now
[vb, vp, Trans, vu, chd] = QUPS2VSX(us, xdc_name, vres, "frames", 1); % make block
vs = link(vb, vres, vp, vu); % link
vs.Trans = Trans; % add Trans
pt1; vs.Media = Media; % add simulation media

% force in simulation mode for testing
vs.Resource.Parameters.simulateMode = 1; % 1 to force simulate mode, 0 for hardware

% save 
filename = char(fullfile("MatFiles","qups-vsx.mat")); 
save(filename, '-struct', 'vs');

% fix Transducer to match VSX
us.xdc = Transducer.Verasonics(Trans);

% save
save(fullfile("MatFiles","qups-conf.mat"), "us", "chd");

VSX;

%% 
