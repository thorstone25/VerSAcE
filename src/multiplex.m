function seq = multiplex(seq, xdc, N)
% MULTIPLEX - Multiplex transmits
% 
% seq = multiplex(seq, xdc, N) creates a Sequence seq for the Transducer
% xdc by replicating and apodizing each transmit aperture to use at most N
% contiguous elements.
% 
% Example:
% xdc = TransducerArray.L12_5v();
% seq = SequenceRadial('type', 'PW', 'angles', -20 : 2 : 20);
% seq = multiplex(seq, xdc, 128); % Vantage UTA-260D -> 128 channels
%
% 

arguments
    seq (1,1) Sequence
    xdc (1,1) Transducer
    N (1,1) {mustBePositive, mustBeInteger} = 128
end

% pass-by-value semantics
seq = copy(seq); 

% sizing
Mx = ceil(xdc.numel / N ); % multiplexing factor ( e.g. 1, 2, 3)
M  =      xdc.numel / Mx ; % active aperture size (e.g. 128, 96, 64)

% active aperture mask (N x Mx)
msk  = cell2mat(arrayfun(@(i){circshift((1:xdc.numel)'<=M,i*M)}, 0:Mx-1));

% multiplexer function
mplx = @(x) repmat(swapdim(x,2,3), [1,Mx,1]); % multiplexed phase/amp (N x S) -> (N x Mx x S)

% set apodization - duplicate and mask
if ~isempty(seq.apodizationf_) % function - leave it alone
    % TODO: warning?
else % value - modify it
    apd = seq.apodization(xdc);
    seq.apodization_ = reshape(msk.*mplx(apd),xdc.numel,[]); % (N x [Mx x S])
end

% set delays - duplicate only
if ~isempty(seq.delaysf_) % function - leave it alone
    % TODO: warning?
elseif ~isempty(seq.delaysv_) % custom value - modify it
    seq.delays_ = reshape(mplx(seq.delays_),xdc.numel,[]); % (N x [Mx x S])
else % default computation - leave it alone
end

% replicate the transmit foci
seq.focus = reshape(repmat(seq.focus, [Mx,1]), 3, []); % replicate foci
if seq.type == "FSA" && isfinite(seq.numPulse), seq.numPulse = Mx * seq.numPulse; end

end



