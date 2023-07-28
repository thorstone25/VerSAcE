function [Recon, ReconInfo] = setVSXLUT(Recon, ReconInfo, tau_rx, tau_tx, fc, apod_tx)
arguments
    Recon (1,1) {mustBeA(Recon, ["struct", "VSXRecon"])}
    ReconInfo (:,:) {mustBeA(ReconInfo, ["struct", "VSXReconInfo"])}
    tau_rx (:,:,:,:,1) double % delays: pixels x rx
    tau_tx (:,:,:,1,:) double % delays: pixels x tx
    fc double = 1; % transmit frequency
    apod_tx (:,:,:,1,:) double {mustBeInRange(apod_tx, -8388608, 8388608)} = ones(size(tau_tx)) % weights: pixels x tx % NOTE: this can clip
end

M   = size(tau_tx,5); % get sizing
N   = size(tau_rx,4);
IM  = prod(size(tau_tx,1:3)); 
IN  = prod(size(tau_rx,1:3));
I   = IM;
assert(IM == IN, "Transmit and receive delay table pixel sizing mismatch.")
assert(isequal(size(apod_tx,1:5) , size(tau_tx,1:5)), "Apodization array and transmit delay table must have identical sizes.");
if any(apod_tx > (double(intmax('int32'))+0.5)/256, 'all'), warning("Positive weights at maximum representable value will be clipped."); end

% tx side look-up table: [address; weight; delays] x pixels per tx
% TODO: disable pixels with 0 weight and/or coordinate with PData.Regions
addr    = int32(0:I-1)'; % all pixels for now 
wght    = int32(round(double(reshape(apod_tx, [I, M]) * 256))); % cast weight
tau_tx  = max(0,int32(round(double(reshape(tau_tx , [I, M]) * fc) * 16))); % convert and cast delays

% store for each tx
for m = 1:M
    i = [ReconInfo.txnum] == m; % find all matching txs
    lut = [addr, wght(:,m), tau_tx(:,m)]'; % make LUT
    [ReconInfo(i).LUT] = deal(lut); % assign
end

% rx side look-up table: pixel x elements in units of wavelengths
Recon.rcvLUT = uint16(round(double(reshape(tau_rx, [I, N]) * fc) * 16));

