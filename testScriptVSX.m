block = VSXBlock();
tx1 = VSXTX();
tx2 = VSXTX(); 
wv1 = 5; % VSXTW(); % DEBUG


event1 = VSXEvent();
event1.tx = tx1;
event1.tx.waveform = wv1;
% event1.rcv = 5;
% event1.process = 5;
% event1.seqControl = 5;

event2 = VSXEvent();
event2.tx = tx2;
event2.tx.waveform = wv1;
% event2.rcv = 6;
% event2.process = 6;
% event2.seqControl = 6;

block.vsxevent = [event1, event2];

[Event, TX, Receive, Process, SeqControl] = block.link();


%% TESTING
TXexp = [struct(tx1), struct(tx2)];
for i = 1:numel(TX)
   for f = string(fieldnames(TX(i)))'
       if (TX(i).(f) == TXexp(i).(f))
           disp("pass");
       end
   end
end

assert(Event(1).tx == 1)
assert(Event(2).tx == 2)



disp(Event);
disp(TX);
disp(Receive);
disp(Process);
disp(SeqControl);
