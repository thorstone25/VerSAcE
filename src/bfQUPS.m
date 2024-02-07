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
    prms.k = [rcv.aperture] - 1; % amount to shift data
    
    % remove pilot pulse t0s
    sz = size(chd0.data);
    R = sz(2) / size(chd0.t0,2); % ratio
    chd0.t0 = chd0.t0(~ismember(1:size(chd0.t0,chd0.mdim), QUPS_BF_PARAMS.ppi));
    sz(2) = R * size(chd0.t0,2); % new size (same ratio)
    chd0.data = zeros(sz, 'like', chd0.data);
    
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
        case 1
            chdo = ChannelData.Verasonics({RData}, vs.Receive(prms.rxs), vs.Trans);
            chd.data(:,:,:,:,1) = chdo.data; % allocation trick
        case 2
            chdo = ChannelData.Verasonics({RData}, vs.Receive(prms.rxs), vs.Trans);
            chd.data = single(chdo.data);
        case 3
            chd.data(:,:,1:128,:,1) = reshape(RData(prms.i,:,:), chd0.T, chd0.M, 128, []); % (T x M x N' x F)
            % x = cat(3,  chd.data, zeros(size(x,1:4)+[0,0,chd.N-2*128,0])); %  (T x M x N x F) pad to full size in elements
            
            for j = unique(prms.k) % each unique tx aperture
                if j % skip 0
                    chd.data(:,prms.k == j,:,:) = circshift(chd.data(:,prms.k == j,:,:), j, 3); % rx shift data over in the matching aperture
                end
            end
            % chd.data = cast(x, 'like', chd.data); % load/cast
    end    
    
    % demultiplex
    if rx_multi || tx_multi
        % de-multiplex the time vector
        if tx_multi
            [t0, dim] = deal(chd.t0, chd.mdim);
            i = 1:2:size(t0,dim); % de-multiplexed indices
            assert(isalmostn(sub(t0,i+0,dim),sub(t0,i+1,dim)));
            chd.t0 = sub(t0,i,dim);
        end
        
        % de-multiplex the data
        Mx = (1 + rx_multi) * (1 + tx_multi); % multiplexing size
        sz = size(chd.data); % og size
        sz = [sz(1:chd.mdim-1), Mx, chd.M/Mx, sz(chd.mdim+1:end)]; % mulit in mdim (txs)
        chd.data = reshape(sum(reshape(chd.data, sz), chd.mdim), [sz(1:chd.mdim-1),sz(chd.mdim+1:end)]);
    end
    toc;
    
    % pre-process
    chd = hilbert(filter(chd, D));
    
    % get image
    b = DAS(us, chd);
    b = double(gather(abs(b)));
end


