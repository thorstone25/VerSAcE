function RFDataProc(RData)
global TOGGLE_RFDataProc;
persistent us chd0 D; % UltrasoundSystem, ChannelData, filter
persistent him a; % image handle, apodization
persistent b k PRE_ARGS POST_ARGS; % DAS arguments
persistent isHadamard; % is a hadamard encoded transmit?

% default to false, just in case
if isempty(TOGGLE_RFDataProc), TOGGLE_RFDataProc = false; end

% load data and pre-allocate outputs the first time
if TOGGLE_RFDataProc && (isempty(us) || isempty(chd0) || isempty(him) || ~isvalid(him))
    disp("Initializing ... ");
    
    % init system config and data
    dat = load("qups-conf.mat", "us", "chd"); % load classes
    [us, chd] = deal(dat.us, dat.chd); % assign
    
    % HACK: fix the start time for the pulse
    TW = evalin('base', 'TW');
    chd.t0 = chd.t0 - TW.peak ./ us.xdc.fc;
    
    % HACK: fix the buffer length
    Receive = evalin('base', 'Receive');
    sz = size(chd.data);
    sz(1) = mode([Receive.endSample] - [Receive.startSample] + 1);
    chd.data = zeros(sz, 'like', chd.data);
    
    % HACK: check if hadamard transmit
    isHadamard = isequal(us.seq.apodization(us.xdc), hadamard(us.xdc.numel));

    % adjust scan to start at 2mm
    us.scan = copy(us.scan); 
    us.scan.z = us.scan.z + 2e-3;
    
    % allocate data
    chd = (singleT(chd)); % typing
    chd.data(:,:,:,:,1) = 0; % pre-allocation hack
    
    % filtering
    D = chd.getPassbandFilter(us.xdc.bw,31);
    
    % init beamforming args
    a = us.apApertureGrowth(1); % apodization
    [b, k, PRE_ARGS, POST_ARGS] = DAS(us, chd, 'apod', a);
    
    % init display
    figure;
    him = imagesc(us.scan, zeros(us.scan.size), [-50 0]);
    colorbar;
    colormap gray;
    
    % save
    chd0 = chd;
    disp("done!")
end

% process data
if TOGGLE_RFDataProc
    % load data
    L = chd0.T*chd0.M;
    chd0.data(:) = RData(1:L,:);
    
    % process
    chd = (filter(hilbert(chd0),D));
    if isHadamard, chd.data = pagemtimes(chd.data, hadamard(us.xdc.numel)); end
    b   = DAS(us, chd, 'apod', a);
    
    % display
    bim = (mod2db(sum(b(:,:,:),3))); % TODO: handle frames better
    him.CData(:) = gather(bim - max(bim,[],'all')); % normalize to peak
end
return
end
