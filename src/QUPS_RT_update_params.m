function [prms, chd0] = QUPS_RT_update_params(prms, vs, chd0)
% pre-process loading indices
evs = vs.Event(prms.evi(prms.bfi)); % Events
prms.rxs = [evs.rcv]; % corresponding Receive indices
prms.txs = [evs.tx];  % corresponding Transmit indices
rcv = vs.Receive(prms.rxs);

prms.i = cell2mat(cellfun(@colon, {rcv.startSample}, {rcv.endSample}, 'UniformOutput', false));
if isfield(rcv, 'aperture')
    prms.k = [rcv.aperture] - 1; % amount to shift data
else
    prms.k = zeros(1,numel(rcv));
end
if isfield(vs.Trans, 'HVMux')
    prms.aps = vs.Trans.HVMux.ApertureES;
else
    prms.aps = vs.Trans.Connector;
end

% peak delay in wavelengths (periods)
t0 = [vs.TW([vs.TX(prms.txs).waveform]).peak] + 2*vs.Trans.lensCorrection;
t0 = swapdim(t0, 2, chd0.mdim); % move to tx dimension

% channel data input/output size pre/post aperture concatenation
if prms.rx_multi || prms.tx_multi
    % de-multiplex the data and delays
    Mx = (1 + prms.rx_multi) * (1 + prms.tx_multi); % multiplexing size
    sz = size(chd0.data,1:4); % og size
    sz = [sz(1:chd0.mdim-1), Mx, chd0.M/Mx, sz(chd0.mdim+1:end)]; % multi in mdim (txs)
    prms.isz = sz; % input data size
    prms.osz = [sz(1:chd0.mdim-1),sz(chd0.mdim+1:end)]; % output data size
    t0 = t0(1:Mx:end); % delays
end

% add transmit peak delay
chd0.t0 = chd0.t0 - t0 ./ (1e6*vs.Trans.frequency); % wavelengths

% scalarize if possible
for d = 1:ndims(t0) % each dimension
    t0s = sub(chd0.t0,1,d); % slice
    if size(chd0.t0,d) ~= 1 && all(t0s == chd0.t0)
        chd0.t0 = t0s;
    end
end