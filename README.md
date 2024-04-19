# VerSAcE: Verasonics Simple Acquisition Exporter
VerSAcE is an object oriented design wrapper and programming interface for Verasonics Vantage software that allows users to rapidly modify and combine acquisition sequences.

## How to use this Project
You can utilize this project at 3 levels:
* High-level - export a [QUPS](https://github.com/thorstone25/qups) `UltrasoundSystem` configuration to VSX structs
* Mid-level - use the `VSXStruct` wrappers and linking stage to program with a tree structure and avoid indexing errors 
* Low-level - use the `VSXStruct` wrappers for argument validation to avoid typos, then use `struct` to convert properties.

## Quick Start
1. Initialize Vantage
```
cd ~/Projects/Vantage;
activate;
```
2. Open the VerSAcE project
```
openProject ~/Projects/versace;
```
3. Return to the Vantage directory
```
cd ~/Projects/Vantage;
```
4. Open and edit the example test script
```
edit test_sequences;
... modify 'xdc_name' to the current transducer ...
... modify/select the desired Sequences via 'seq*' / 'seq_ind' variables ...
```
5. Run the example script
```
test_sequences;
VSX;
```

## Structure
VerSAcE is designed to simplify acquisition scripts by replacing indexing references with handle classes that can be 'linked' by index in a final processing stage. Each handle class a simple wrapper around the corresponding native struct e.g. `VSXTW` is a wrapper for the `TW` struct, and `VSXTX` is a wrapper for the `TX` struct, and so on. Each wrapper additionally performs argument validation to prevent typos, and in some cases maps enumerated options by name to their corresponding numeric value. Typically, you can use the 'tab' key when only a predefined list of options is acceptable. 

Numeric references to other classes are replaced by instances of the other class. For example, `tx` and `rcv` properties of the `VSXEvent` class must be a `VSXTX` and a `VSXReceive` respectively, and the `waveform` property of the `VSXTX` class must be a `VSXTW`. Linking is handled by the `link` method of the [`VSXBlock`](#the-vsxblock) class.

### The VSXBlock
The `VSXBlock` class provides an abstraction of an acquisition sequence, and performs the linking method. The block separates events into 3 consecutive sets:
* `capture` - the set of events for capturing data
* `post` - the set of post-processing events for the captured data
* `next` - the next event _after_ the block is complete

Typical usage would be to place all events acquiring data and processing data per frame into the `capture` array. The `post` events would contain final image display or data saving events that operate on the entire set of data acquisitions. Finally the `next` event references which event to jump to after the acquisition, and would typically be the first `capture` event of the same block to repeat acquisition. If `next` is an independent event, then no 'jump to X' event is created.

In addition to numeric indexing, the `link` method performs other final pre-processing steps:
* `string`s and `logical`s are converted to their equivalent `char` or `double` types corresponding to the Vantage documentation
* if #elements > #channels, a `TX.aperture` is selected from the available apertures within the `Trans` struct that satisfies `TX.Apod`
* UI controls with a position of 'auto' are assigned a user position
* (optional) transmit power density (TXPD) is computed

### Modifying the VSXBlock
With this structure, one can therefore use traditional indexing to modify and replace objects. Use the `copy` method to represent a different instance. For example, you can implement tissue harmonic imaging by duplicating transmit events, inverting a parametric waveform, and accumulating the negative pulse receive data within the same acquisition:

```
% create a VSXBlock ...
Trans = computeTrans(struct('name', 'P4-2v', 'units', 'mm')); % transducer
us = UltrasoundSystem('xdc', Transducer.Verasonics(Trans)); % defaults
vres = VSXResource(); % global 
vb = QUPS2VSX(us, Trans, vres); % convert to a VSXBlock

% create the positive and negative pulses
vTWp = VSXTW('type','parametric', 'Parameters', [Trans.frequency, 0.67, 2, +1]); % positive pulse
vTWm = VSXTW('type','parametric', 'Parameters', [Trans.frequency, 0.67, 2, -1]); % negative pulse

% create the positive and negative acquisitions
vEvents = vb.capture; % (acq x frames) extract original
vEventsTHI = repmat(vEvents, [2 1]); % copy 2 acq. per acquitision
vEventsTHI(1:2:end, :) =      vEvents ; % reassign original
vEventsTHI(2:2:end, :) = copy(vEvents); %   assign copy

% reference transmits
vTXp = [vEventsTHI(1:2:end, :).tx]; % positive TXs
vTXm = [vEventsTHI(2:2:end, :).tx]; % negative TXs

% assign waveforms
[vTXp.waveform] = deal(vTWp); % assign positive waveform to positive TX
[vTXm.waveform] = deal(vTWm); % assign negative waveform to negative TX

% set the negative pulse acquisitions to in-place accumulation
vRcvm = [vEventsTHI(2:2:end, :).rcv]; % negative Receives
[vRcvm.mode] = deal(1); % set to 'accumulate' mode
```

### Multi-block Sequences
To combine multiple acquisitions sequences, simply set the `next` property to the first `capture` of a corresponding block, e.g.

```
% setup UltrasoundSystem objects via QUPS
xdc = TransducerVerasonics("P4-2v");
seq = [ ...
    SequenceRadial('type', 'PW', 'angles', -25 : 2.5 : 25), ...
    Sequence(      'type', 'FSA', 'numPulse', xdc.numel) ...
];
us = [
    UltrasoundSystem('xdc', xdc, 'seq', seq(1)),
    UltrasoundSystem('xdc', xdc, 'seq', seq(2)) ...
];

% create VSXBlocks
vb    = QUPS2VSX(us(1), ...); % 1st acquisition sequence
vb(2) = QUPS2VSX(us(2), ...); % 2nd acquisition sequence

% link
vb(1).next = vb(2).capture(1); % after 1st, go to 2nd
vb(2).next = vb(1).capture(1); % after 2nd, go to 1st

% export
vs = vb.link();
filename = 'MatFiles/qups-vsx.mat';
save(filename, '-struct', 'vs');

% run
VSX;
```

## Setting up a MATLAB project for path management
It may be helpful to create a MATLAB project to manage your QUPS, VerSAcE, and Vantage paths.
1. Change directory to your versace repository
2. Open a matlab editor, navigate to the "Home" tab and select New->Project->From Folder->Create (were vsx-ood path is chosen)
3. In the "Project" tab, select References->Browse->*navigate to qups repo*->Add
4. In the "Project" tab, select References->Browse->*navigate to Vantage repo (e.g., Vantage-4.8.4)*->Add
5. In your "Current Folder" MATLAB pane, double-click the "Versace.prj" file. You will need to open the file twice due to some MATLAB bug. Now, you should have all necessary QUPS, VerSAcE, and Vantage paths added to your MATLAB $PATH environment variable.
6. Change directory to your Vantage directory.
7. In the "Home" tab, select Open->*navigate to versace/scripts dir*->*Choose custom script*. You should be able to run this custom script while your current working dir is still on Vantage.

