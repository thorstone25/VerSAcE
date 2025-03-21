function bim = RFDataImg(RData, conf, kwargs)
arguments
    RData {mustBeNumeric}
    conf struct {mustBeScalarOrEmpty} = struct.empty
    kwargs.display (1,1) logical = true;
    kwargs.verbose (1,1) logical = true;
end

global TOGGLE_RFDataImg; %#ok<GVMIS> 
persistent us chd0 D; % UltrasoundSystem, ChannelData, filter
persistent a; % apodization
% persistent b k PRE_ARGS POST_ARGS; % DAS arguments
persistent isHadamard; % is a hadamard encoded transmit?

% default to false, just in case
if isempty(TOGGLE_RFDataImg), TOGGLE_RFDataImg = false; end

% load data and pre-allocate outputs the first time
if ~isempty(conf) || ( true ... either post-processing or triggered and
    && (isempty(us) || isempty(chd0))) % if conf or image not initialized ...
    if kwargs.verbose, disp("Initializing ... "); end
    
    % init system config and data
    if all(isfield(conf, ["us", "chd"]))
        dat = conf; % use input
    else
        dat = load("qups-conf.mat", "us", "chd"); % load classes
    end
    [us, chd] = deal(dat.us, dat.chd); % assign

    % HACK: check if hadamard transmit
    isHadamard = isequal(us.seq.apodization(us.xdc), hadamard(us.xdc.numel));

    % allocate data
    chd = (singleT(chd)); % typing
    chd.data(:,:,:,:,1) = 0; % pre-allocation hack
    
    % filtering
    D = chd.getPassbandFilter(us.xdc.bw,51);
    
    % init beamforming args
    a = us.apApertureGrowth(1); % apodization
    
    % save
    chd0 = chd;
    if kwargs.verbose, disp("done!"); end
end

% process data
if ~isempty(conf) || TOGGLE_RFDataImg
    % load data
    L = chd0.T*chd0.M;
    chd0.data(:) = RData(1:L,:);
    
    % process
    chd = gpuArray(filter(hilbert(chd0),D));
    if isHadamard, chd.data = pagemtimes(chd.data, hadamard(us.xdc.numel)); end
    b   = DAS(us, chd, 'apod', a);
    
    % display
    bim = gather(double(abs(sum(b(:,:,:),3)))); % TODO: handle frames better
else
    bim = zeros(us.scan.size);
end
end
