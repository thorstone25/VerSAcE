function b = bfQUPS(RData)

% return a b mode image to VSX
persistent vs us chd0 D prms;

global TOGGLE_bfQUPS;
global QUPS_BF_PARAMS; 

% on trigger
if isempty(TOGGLE_bfQUPS), TOGGLE_bfQUPS = false; end

% init
if isempty(us) || isempty(chd0)
    load('qups-conf', 'us', 'chd'); % load
    chd0 = chd;
    try chd0 = gpuArray(chd0); end %#ok<TRYNC>
    prms = QUPS_BF_PARAMS; % load params into persistent memory
    vs = update_vstruct(); % load configuration
        
    % params
    D = chd0.getPassbandFilter(us.xdc.bw);
    
    % pre-process loading indices
    evs = vs.Event(QUPS_BF_PARAMS.evi(QUPS_BF_PARAMS.bfi)); % Events
    prms.rxs = [evs.rcv]; % corresponding Receive indices
    rcv = vs.Receive(prms.rxs);
    
    prms.i = cell2mat(cellfun(@colon, {rcv.startSample}, {rcv.endSample}, 'UniformOutput', false));
    if isfield(rcv, 'aperture')
        prms.k = [rcv.aperture] - 1; % amount to shift data
    else
        prms.k = zeros(1,numel(rcv));
    end
    prms.aps = vs.Trans.HVMux.ApertureES;
    
    % channel data input/output size pre/post aperture concatenation
    if prms.rx_multi || prms.tx_multi
        % de-multiplex the data
        Mx = (1 + prms.rx_multi) * (1 + prms.tx_multi); % multiplexing size
        sz = size(chd0.data,1:4); % og size
        sz = [sz(1:chd0.mdim-1), Mx, chd0.M/Mx, sz(chd0.mdim+1:end)]; % multi in mdim (txs)
        prms.isz = sz; % input data size
        prms.osz = [sz(1:chd.mdim-1),sz(chd.mdim+1:end)]; % output data size
    end    
end

if ~TOGGLE_bfQUPS
    b = zeros(us.scan.size);
else
    % settings
    rx_multi = prms.rx_multi; % rx multiplexing
    tx_multi = prms.tx_multi; % tx_multiplexing

    % load
    chd = copy(chd0);
    tic;
    switch 3 % different ways of loading
        case 1 % safest - parse and load everything, then cast to single
            chdo = ChannelData.Verasonics({RData}, vs.Receive(prms.rxs), vs.Trans);
            chd.data = single(chdo.data);
        case 2 % two-step - parse and load everything, then input to implicit cast
            chdo = ChannelData.Verasonics({RData}, vs.Receive(prms.rxs), vs.Trans);
            chd.data(:,:,:,:,1) = chdo.data; % allocation trick
        case 3 % optimal - load relevant data into pre-allocated array, then shift to parse
            % load data 
            x = reshape(RData(prms.i,:,:), chd0.T, chd0.M, 128, []); % (T x M x N' x F)
            
            % send to ChannelData object
            for j = unique(prms.k) % each unique tx aperture
                k = prms.k == j; % matching tx indices
                a = prms.aps(:,j+1); % aperture indices
                i = logical(a); % active rx elements
                chd.data(:,k,i,:,1) = x(:,k,a(i),:); % send data to matching receive elements
            end
    end    
    
    % demultiplex by summing left/right apertures
    if rx_multi || tx_multi
        chd.data = reshape(sum(reshape(chd.data, prms.isz), chd.mdim), prms.osz);
    end
    disp("Data loaded in " + toc() + " seconds."); % loading time
    
    % pre-process
    chd = hilbert(filter(chd, D));
    
    % get image
    b = DAS(us, chd);
    b = double(gather(abs(b)));
end


