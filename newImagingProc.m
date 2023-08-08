function chd = newImagingProc(RData, TW, Receive)
global TOGGLE_RFDataProc;
persistent us chd0 D; % UltrasoundSystem, ChannelData, filter
persistent him a; % image handle, apodization
persistent b k PRE_ARGS POST_ARGS; % DAS arguments

% default to false, just in case
if isempty(TOGGLE_RFDataProc), TOGGLE_RFDataProc = false; end

% load data and pre-allocate outputs the first time
if TOGGLE_RFDataProc && (isempty(us) || isempty(chd0) || isempty(him) || ~isvalid(him))
    disp("Initializing ... ");
    
    % init system config and data
    dat = load("qups-conf.mat", "us", "chd"); % load classes
    [us, chd] = deal(dat.us, dat.chd); % assign

    chd.data(:, :, :, :, 1) = chd.data;
    
    % downsampling
    us.xdc.pitch = us.xdc.pitch / 2;
    us.xdc.numel = us.xdc.numel / 2;
    chd.data = chd.data(:, 1:2:end, 1:2:end, :, :);
    chd.t0 = chd.t0(:, 1:2:end);
    us.seq.numPulse = us.seq.numPulse / 2;


    % HACK: fix the start time for the pulse
    if nargin < 2, TW = evalin('base', 'TW'); end
%     [~,TW.peak] = computeTWWaveform(TW);
%     chd.t0 = chd.t0 - TW.peak ./ us.xdc.fc;
    
    % HACK: fix the buffer length
    if nargin < 3, Receive = evalin('base', 'Receive'); end
    sz = size(chd.data);
%     sz(1) = mode([Receive.endSample] - [Receive.startSample] + 1);
    chd.data = zeros(sz, 'like', chd.data);
    
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
    RData = reshape(RData(1:2:end, 1:2:end),[], 64, 64); % reshape rdata
%     RData = permute(RData, [3, 1, 2]); % rearrange dimensions
    chd0.data(:) = RData(1:3200, :);
    
    % process
    chd = (filter(hilbert(chd0),D));
    b   = DAS(us, chd, 'apod', a);
    
    % display
    bim = (mod2db(sum(b(:,:,:),3))); % TODO: handle frames better
    him.CData(:) = gather(bim - max(bim,[],'all')); % normalize to peak
end
return
end
