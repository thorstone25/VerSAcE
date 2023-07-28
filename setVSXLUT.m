function [Recon, ReconInfo] = setVSXLUT(Recon, ReconInfo, PData, tau_rx, tau_tx, fc, apod_tx)
arguments
    Recon (1,1) {mustBeA(Recon, ["struct", "VSXRecon"])}
    ReconInfo (:,:) {mustBeA(ReconInfo, ["struct", "VSXReconInfo"])}
    PData (1,1) {mustBeA(PData, ["struct", "VSXPData"])}
    tau_rx (:,:,:,:,1) double % delays: elements x pixels (!!!!!)
    tau_tx (:,:,:,1,:) double % delays: pixels x tx
    fc double = 1; % transmit frequency
    apod_tx (:,:,:,1,:) double {mustBeInRange(apod_tx, -8388608, 8388608)} = ones(size(tau_tx)) % weights: pixels x tx % NOTE: this can clip
end

M   =      size(tau_tx,5); % get sizing
N   =      size(tau_rx,4);
IM  = prod(size(tau_tx,1:3)); 
IN  = prod(size(tau_rx,1:3));
I   = IM;
assert(numel(ReconInfo) == numel(PData.Region));
assert(IM == IN, "Transmit and receive delay table pixel sizing mismatch.")
assert(isequal(size(apod_tx,1:5) , size(tau_tx,1:5)), "Apodization array and transmit delay table must have identical sizes.");
if any(apod_tx > (double(intmax('int32'))+0.5)/256, 'all'), warning("Positive weights at maximum representable value will be clipped."); end

wght    =       round(double(reshape(apod_tx, [I, M]) * 256)); % cast weight
tau_tx  = max(0,round(double(reshape(tau_tx , [I, M]) * fc) * 16)); % convert and cast delays

% tx side look-up table: [address; weight; delays] x pixels per ReconInfo (tx)
for i = 1:numel(ReconInfo) % for each region
    addr = PData.Region(i).PixelsLA'; % all pixels for now
    ReconInfo(i).LUT = int32([addr, wght(:,i), tau_tx(:,i)]'); % make LUT
end

% rx side look-up table: elements x pixels in units of wavelengths
Recon.rcvLUT = uint16(round(double(reshape(tau_rx, [I, N])' * fc) * 16));

