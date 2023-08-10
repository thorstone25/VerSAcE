function b = imagingProc(RData, conf, kwargs)
arguments
    RData {mustBeNumeric}
    conf struct {mustBeScalarOrEmpty} = struct.empty
    kwargs.display (1,1) logical = true;
    kwargs.verbose (1,1) logical = true;
end

global TOGGLE_imagingProc;
persistent us chd0 D; % UltrasoundSystem, ChannelData, filter
persistent a isHadamard; % apodization, hadamard flag
% persistent b k PRE_ARGS POST_ARGS; % DAS arguments

% default to false, just in case
if isempty(TOGGLE_imagingProc), TOGGLE_imagingProc = false; end

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

    % downsampling
    us.xdc.pitch = us.xdc.pitch / 2;
    us.xdc.numel = us.xdc.numel / 2;
    % chd.data = chd.data(:, 1:2:end, 1:2:end, :, :);
    % chd.t0 = chd.t0(:, 1:2:end);
    % us.seq.numPulse = us.seq.numPulse / 2;
 
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
if TOGGLE_imagingProc
    % load data
    L = chd0.T*chd0.M;
    chd0.data(:) = RData(1:L,:);
    
    % move to gpu
    chd = copy(chd0);
    chd = gpuArray(chd);
    
    % downsample rx
    chd.data = sub(chd.data,1:2:chd.N, chd.ndim);

    % process
    chd = (filter(hilbert(chd),D));
    b   = DAS(us, chd, 'apod', a);
    
    % display
    % bim = (mod2db(sum(b(:,:,:),3))); % TODO: handle frames better
    % him.CData(:) = gather(bim - max(bim,[],'all')); % normalize to peak
    
    % casting
    b = gather(double(abs(b)));
else
    b = zeros(us.scan.size);
end
return
end
