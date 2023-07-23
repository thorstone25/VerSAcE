%% Create different system configurations
c0 = 1500; % sound speed
xdc = TransducerArray.L11_5v();
uss = UltrasoundSystem('xdc', xdc); 
uss.seq.c0 = c0;
[uss.scan.dx, uss.scan.dz] = deal(uss.lambda / 4);


% sequences
pf = [0;0;50e-3] + [1e-3;0;0] .* (-20 : 5 : 20);
seqfsa = Sequence(     'type', 'FSA', 'c0', c0, 'numPulse', xdc.numel);
seqpw = SequenceRadial('type', 'PW' , 'c0', c0, 'angles', -15 : 5 : 15); 
seqfc = Sequence(      'type', 'VS' , 'c0', c0, 'focus', pf);

uss = copy(repmat(uss, [1,3]));
[uss.seq] = deal(seqfsa, seqpw, seqfc);

%%
vres = VSXResource();
for i = 3:-1:1
    [vb(i), vp(i), Trans, vu(i)] = QUPS2VSX(uss(i), "L11-5v", vres, "frames", 2);
end
vs = link(vb, vres, vp, vu);

%% 