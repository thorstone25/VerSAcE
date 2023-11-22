% Create an instance of the VSXBlock class
block = VSXBlock();

% Create some sample VSXEvent objects and assign values to their properties
event1 = VSXEvent();
event1.tx = 1;
event1.rcv = 2;
event1.process = 1;
event1.seqControl = 1;

event2 = VSXEvent();
event2.tx = 2;
event2.rcv = 1;
event2.process = 2;
event2.seqControl = 2;

% Assign the VSXEvent objects to the vsxevent property of the block
block.vsxevent = [event1, event2];

% Call the link method to retrieve the linked components
[Event, TX, Receive, Process, SeqControl] = block.link();

% Display the properties of the linked components
disp(Event);
disp(TX);
disp(Receive);
disp(Process);
disp(SeqControl);
