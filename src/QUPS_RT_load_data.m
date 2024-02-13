function chd = QUPS_RT_load_data(RData, chd0, vs, prms)
% QUPS_RT_load_data - Load data in real-time
%
% chd = QUPS_RT_load_data(RData, chd0, vs, prms) loads the data RData using
% the template ChannelData chd0, the Verasonics param struct vs with VSX
% structs (e.g. vs.Receive, vs.Trans), and other parameters in prms,
% including  prms.rx_multi, prms.tx_multi, prms.rxs, etc.
%
% This function is only partially documented, and subject to change.
%

% load
chd = copy(chd0);
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
if prms.rx_multi || prms.tx_multi
    chd.data = reshape(sum(reshape(chd.data, prms.isz), chd.mdim), prms.osz);
end
