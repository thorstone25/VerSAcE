function [vEvent, vPData] = addReconVSX(scan, vTX, vRcv, vResource, vSeq, kwargs)
arguments
    scan Scan
    vTX VSXTX
    vRcv VSXReceive
    vResource VSXResource {mustBeScalarOrEmpty} = VSXResource.empty
    vSeq (1,:) VSXSeqControl = VSXSeqControl('command', 'returnToMatlab')
    kwargs.numFrames (1,1) double = 1
    kwargs.multipage (1,1) logical = false
    kwargs.display (1,1) logical = true
    kwargs.apod (:,:,:,1,:,:) {mustBeNumericOrLogical} = 1
end

%% Validate sizing
Tx = prod(size(vRcv,2));
assert(all(any(size(kwargs.apod,1:6) == [[scan.size,1,Tx,kwargs.numFrames]; ones(1,6)],1)), ...
    "The size of the apodization(" +...
    join(string(size(kwargs.apod,1:6)),", ") +...
    ") must be compatible with the scan (" + ...
    join(string(scan.size),", ") +...
    ") and the number of transmits (" + Tx + ")." ...
);

%% ReconInfo
% We need 1 ReconInfo structures for each transmit
vReconInfo = copy(repmat(VSXReconInfo('mode', 'accumIQ'), size(vRcv,1:2)));  % default is to accumulate IQ data.

% - Set specific ReconInfo attributes.
for k = 1%:size(vReconInfo,3)
for i = 1:size(vReconInfo,1) % for each rx of the reconinfo object (multiplexing)
for j = 1:size(vReconInfo,2) % for each tx of the reconinfo object
    [vReconInfo(i,j,k).txnum    ] = deal(vTX(   j)); % tx
    [vReconInfo(i,j,k).rcvnum   ] = deal(vRcv(i,j)); % rx
    [vReconInfo(i,j,k).regionnum] = deal(       j );% region index TODO: VSXRegion
    if kwargs.multipage
        [vReconInfo(i,j,k).pagenum] = deal(sub2ind(size(vReconInfo),i,j,k)); % page number
    end
end
end
end

vReconInfo( 1 ).Pre  = "clearInterBuf"; % clear accumulator on init
vReconInfo(end).Post = "IQ2IntensityImageBuf"; % copy result 

% vReconInfo( 1 ).mode = 'replaceIQ'; % on first tx, replace IQ data
% vReconInfo(end).mode = 'accumIQ_replaceIntensity'; % on last tx, accum and "detect"

% threadsync conflict? Overlapping regions are processed in parallel, so
% this introduces a race conditions if they overlap
ts = (any( (sum(kwargs.apod,5:6) > 1) .* kwargs.apod , 1:3)); % Tx x F
[vReconInfo.threadSync] = deal(any(ts,'all')); %%% TODO: can we sort into independent groups, shuffle to trick the FPGA, and compute this ind. per region?

% if isscalar(vReconInfo), vReconInfo.mode = 'replaceIntensity'; end % 1 tx

%% PData
% create the pixel data
vPData = VSXPData.QUPS(scan);

% get the apodization in it's full size
ap = logical(kwargs.apod); % mask
Psz = [scan.size, 1, Tx, 1+0*size(vRcv,3)]; % full sizs (I x Tx x [1*|F])
ap = repmat(ap, Psz ./ size(ap,1:6)); % explicit brodcast

% convert to address / count format
aps = num2cell(ap, 1:3); % pack pixels per cell
cnt  = cellfun(@nnz,   aps,  "UniformOutput", false); % number of pixels
addr = cellfun(@find,  aps,  "UniformOutput", false); % linear address (1-based)
addr = cellfun(@int32, addr, "UniformOutput", false); % to int
addr = cellfun(@(x)x-1,addr, "UniformOutput", false); % to 0-based

%  sanity check
assert(all([cnt{:}] <= scan.nPix), "Detected index attempt higher than the number of pixels - this is a bug!");

% create a default region for each transmit
% TODO: VSXRegion
vPData.Region = cellfun(@(LA, NP) struct('Shape', struct('Name','Custom'), 'PixelsLA', LA, 'numPixels', NP), addr(1:Tx), cnt(1:Tx));

% fill out the regions
% TODO: VSXRegion
% vPData.Region = computeRegions(struct(vPData));

%% Make Buffers
% if inter(mediate) buffers needed, create one
if true || any(contains([vReconInfo.mode], "IQ"))
    % get required number of pages
    if any(contains([vReconInfo.mode], "pages"))
        P = vResource.Parameters.numRcvChannels; % for page acquisitions, each page is a receive channel
    else 
        P = max([1, vReconInfo.pagenum]); % use max page number in reconinfo
    end
    
    % make a buffer
    vbuf_inter = VSXInterBuffer.fromPData(vPData, 'numFrames', kwargs.numFrames, 'pagesPerFrame', P); 
else
    vbuf_inter = VSXInterBuffer.empty; % no buffer
end

% if image buffers needed, create one
if any(contains([[vReconInfo.mode], [vReconInfo.Post]], "Intensity"))
    vbuf_im = VSXImageBuffer.fromPData(vPData, 'numFrames', kwargs.numFrames);
else
    vbuf_im = VSXImageBuffer.empty; % no buffer
end

%% Recon
vRecon = VSXRecon('pdatanum', vPData, 'IntBufDest', vbuf_inter, 'ImgBufDest', vbuf_im, "RINums", vReconInfo(:));

%% Display window
if kwargs.display
    vDisplayWindow = VSXDisplayWindow.QUPS(scan, ...
        'Title', 'VSX Beamformer', ...
        'numFrames', kwargs.numFrames, ...
        'AxesUnits', 'mm', ...
        'Colormap', gray(256) ...
        );

    %% Process
    display_image_process = VSXProcess('classname', 'Image', 'method', 'imageDisplay');
    display_image_process.Parameters = {
        'imgbufnum', vbuf_im,...   % number of buffer to process.
        'framenum',-1,...   % (-1 => lastFrame)
        'pdatanum', vPData,...    % PData structure to use
        'pgain',1.0,...            % pgain is image processing gain
        'reject',2,...      % reject level
        'persistMethod','simple',...
        'persistLevel',20,...
        'interpMethod','4pt',...
        'grainRemoval','none',...
        'processMethod','none',...
        'averageMethod','none',...
        'compressMethod','power',...
        'compressFactor',40,...
        'mappingMethod','full',...
        'display',1,...      % display image after processing
        'displayWindow', vDisplayWindow, ...
        };

    %% Event
    vEvent = VSXEvent(...
        'info', 'VSX Recon', ...
        'recon', vRecon, ...
        'process', display_image_process, ...
        'seqControl', vSeq ...
        );
else
    vEvent = VSXEvent.empty;
end

%% Add to required Resource buffer
if ~isempty(vbuf_inter    ), vResource.InterBuffer(end+1)    = vbuf_inter    ; end
if ~isempty(vbuf_im       ), vResource.ImageBuffer(end+1)    = vbuf_im       ; end
if ~isempty(vDisplayWindow), vResource.DisplayWindow(end+1)  = vDisplayWindow; end
