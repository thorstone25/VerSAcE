function [vBlock, vTrans, vPData, vTW, vTX, vRcv, vRecon, display_image_process, vSeqControl, vTGC, vReconInfo, t0] = QUPS2VSX(us, vResource, xdc_name, kwargs)
    arguments
        us (1,1) UltrasoundSystem
        vResource (1,1) VSXResource
        xdc_name (1,1) string % transducer name
        kwargs.units (1,1) string {mustBeMember(kwargs.units, ["mm", "wavelengths"])} = "mm"
        kwargs.vTW (1,1) VSXTW = VSXTW('type', 'parametric', 'Parameters', [us.xdc.fc/1e6, 0.67, 1, 1]); %-
        kwargs.vTGC (1,1) VSXTGC = VSXTGC('CntrlPts', [0,297,424,515,627,764,871,1000],...
                                          'rangeMax', hypot(us.scan.zb(2), us.scan.xb(2)) ./ us.lambda); %-
    end

    % squash obj to struct warning
    warning_state = warning('off', 'MATLAB:structOnObject');
    
    %% create a scan from us.scan
    scan = us.scan; % should be ScanCartesian
    assert(isa(scan, 'ScanCartesian'));
    
        %% Trans
    if isscalar(xdc_name)
            vTrans.name = char(xdc_name);
            vTrans.units = char(kwargs.units);
            vTrans = computeTrans(vTrans);
    else
        vTrans.spacing          = us.xdc.pitch;
        vTrans.Bandwidth        = us.xdc.bw;
        vTrans.elementWidth     = us.xdc.width;
        vTrans.elementHeight    = us.xdc.height;
        vTrans.numelements      = us.xdc.numel;
        vTrans.frequency        = us.xdc.fc; 
        vTrans.elevationFocusMm = us.xdc.el_focus;

    end
    
    fs_decim = 20.8333; % decimation frequency - this is the supported frequency that verasonics will choose based on the ...
    ...transducer frequency, sample mode, and frequency support list see the Vantage 4.3.0 programming manual, page 106
    c0 = vResource.Parameters.speedOfSound; % I think?
    lambda = c0 / (1e6*vTrans.frequency); % wavelengths
    
    %% ---------------------------- Sequence params -------------------------- %
    % Parameters for Focused Transmit Imaging Sequence
    P.startDepth = 2;   %            acquisition start (wavelengths)
    P.endDepth = 256;   % (requested) acquisition end (wavelengths)
    P.numTx = 26;       % Number of active transmit elements in the aperture.
    P.numRays = 192;    % Number of ray lines in frame.
    P.txFocus = 130;    % transmit focal pt in wavelengths
        
    % set the maximum acquisition depth
    maxAcqLength = ceil(sqrt(P.endDepth^2 + ((vTrans.numelements-1) * vTrans.spacing)^2));
    spw = fs_decim / vTrans.frequency; % samples per wave
    bufLen = 256 * ceil((maxAcqLength - P.startDepth) * spw / 256); % only modulus 256 sample buffer length supported
    
    %% Allocate buffers and Set Parameters
    toggle = true;
    F = 1; % number of frames
    T = 256 * 20; %% bufLen;
    
    na = 7; %%us.seq.numPulse; % Set na = number of angles.
    if (na > 1)
        dtheta = (36*pi/180)/(na-1);
        startAngle = -36*pi/180/2;
    else
        dtheta = 0;
        startAngle=0;
    end % set dtheta to range over +/- 18 degrees.

    vResource.Parameters.speedOfSound = us.seq.c0; 
    vResource.Parameters.verbose = 2; %%
    vResource.RcvBuffer(end+1)   = VSXRcvBuffer('rowsPerFrame', T * us.seq.numPulse,...  %-
                                                'colsPerFrame', vResource.Parameters.numRcvChannels,...
                                                'numFrames', F);
    vResource.InterBuffer(end+1) = VSXInterBuffer('numFrames', F);
    vResource.ImageBuffer(end+1) = VSXImageBuffer('numFrames', F);
  
    %rectangular aperture only
    
   %% PData
    scan = scale(us.scan, 'dist', 1./lambda);
    vPData = VSXPData();
%     vPData.PDelta = [0.5, 0, 0.5]; 
    vPData.PDelta = [scan.dx, 0, scan.dz]; 
    vPData.Size(1) = scan.nz; %-
    vPData.Size(2) = scan.nx; %-
    vPData.Size(3) = scan.ny; %-
    vPData.Origin = [scan.xb(1), scan.yb(1), scan.zb(1)];
    
% %     vPData.Region = computeRegions(struct(vPData));

% %     vResource.DisplayWindow(end+1) = VSXDisplayWindow('ReferencePt', vPData.Origin);
    vResource.DisplayWindow(1).Title = 'L11-5vFlashAngles';
    vResource.DisplayWindow(1).pdelta = scan.dx; % 0.35;
%     vResource.DisplayWindow(1).Position = [250,89.500000000000000,542,535];
    vResource.DisplayWindow(1).Position = [250,89.500000000000000,scan.nx,scan.nz];
    vResource.DisplayWindow(1).ReferencePt = [vPData(1).Origin(1),0,vPData(1).Origin(3)];   % 2D imaging is in the X,Z plane
    vResource.DisplayWindow(1).numFrames = F;
    vResource.DisplayWindow(1).AxesUnits = 'mm';
    vResource.DisplayWindow(1).Colormap = gray(256);

     %% TW
    vTW = kwargs.vTW;

    %% TX
    vTX.Origin = [0.0,0.0,0.0]; %%
    vTX.focus = 0.0; %%
    vTX.Delay = zeros(1,vTrans.numelements); %%
    vTX = repmat(VSXTX(), [1, us.seq.numPulse]);
    vTX = sort(copy(vTX)); 
    [vTX.waveform] = deal(vTW);
    
    % - Set event specific TX attributes.
    if fix(na/2) == na/2       % if na even
        P.startAngle = (-(fix(na/2) - 1) - 0.5)*dtheta;
    else
        P.startAngle = -fix(na/2)*dtheta;
    end
    
    delay  = - us.seq.delays(us.xdc) * us.xdc.fc;
    apod = us.seq.apodization(us.xdc);
    t0 = min(delay, [], 1); % min over elements
    delay = delay - t0;
    for i = 1:us.seq.numPulse
        if us.seq.type == "VS"
            vTX(i).FocalPt = us.seq.focus(:,i)/lambda;
%             vTX(i).Origin(1) = us.seq.focus(1,i)/lambda;
%             vTX(i).focus     = us.seq.focus(3,i)/lambda;
        elseif us.seq.type == "PW"
            vTX(i).Steer = deg2rad([us.seq.angles(i), 0]);
        elseif us.seq.type == "FSA"
            %do nothing
        end
        for j = 1:us.tx.numel
            vTX(i).Apod(j) = apod(j,i);
            vTX(i).Delay(j) = delay(j,i); 
        end
    end
    for n = 1:na   % na transmit events
        vTX(n).Steer = [(P.startAngle+(n-1)*dtheta),0.0];
%         vTX(n).Delay = computeTXDelays(struct(vTX(n))); % requires base
        % workspace variables
    end
    
    for i = 1:us.seq.numPulse
%         vTX(i).TXPD = computeTXPD(struct(vTX(i)), struct(vPData));
    end

    %% TGC
    vTGC = kwargs.vTGC;
    vTGC.Waveform = computeTGCWaveform(vTGC, 1e6*vTrans.frequency);

    %% Rcv 
    % get temporal sampling region
    dnear = 2 * scan.zb(1);   
    dfar  = 2 * hypot(range(scan.xb), scan.zb(end));    
    
    vRcv = VSXReceive();
%     vRcv.aperture = 1;
    vRcv.Apod = ones(1, [vResource.Parameters.numRcvChannels]);
    vRcv.startDepth = floor(dnear); %% 2
    vRcv.endDepth = ceil(dfar); %% 256
    vRcv.TGC = kwargs.vTGC(1);
    vRcv.bufnum = vResource.RcvBuffer;
    vRcv.framenum = 1;
    vRcv.acqNum = 1; 
    vRcv.callMediaFunc = 0;
    
    vRcv = repmat(vRcv,1,us.seq.numPulse);
    vRcv = sort(copy(vRcv));
    
    vRcv(1).callMediaFunc = 1; %%
    % - Set event specific Receive attributes.
    for i = 1:vResource.RcvBuffer(1).numFrames 
        for j = 1:na
            vRcv(j).framenum = i;
            vRcv(j).acqNum = j;
        end
    end

    %% Recon
    vRecon = VSXRecon();
    vRecon.senscutoff = 0.6;
    vRecon.pdatanum = length(vPData);
    vRecon.rcvBufFrame = -1;
    vRecon.IntBufDest = [1,1];
    vRecon.ImgBufDest = [1,-1];
    
    % Define ReconInfo structures.
    % We need na ReconInfo structures for na steering angles.
    vReconInfo = repmat(VSXReconInfo('mode', 'accumIQ', ...  % default is to accumulate IQ data.
                       'txnum', 1, ...
                       'rcvnum', 1, ...
                       'regionnum', 1), 1, na);
    vReconInfo = sort(copy(vReconInfo));
    % - Set specific ReconInfo attributes.
    if na>1
        vReconInfo(1).mode = 'replaceIQ'; % replace IQ data
        for j = 1:na  % For each row in the column
            vReconInfo(j).txnum = j;
            vReconInfo(j).rcvnum = j;
        end
        vReconInfo(na).mode = 'accumIQ_replaceIntensity'; % accum and detect
    else
        vReconInfo(1).mode = 'replaceIntensity';
    end
    
    % associate the reconinfo
    vRecon.RINums = vReconInfo; %%

    %% Process
    display_image_process = VSXProcess();
    display_image_process.classname = 'Image';
    display_image_process.method = 'imageDisplay';
    display_image_process.Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
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
                         'displayWindow',1};
    
    save_rf_data = VSXProcess();
    save_rf_data.classname = 'External';
    save_rf_data.method = 'RFDataStore';
    save_rf_data.Parameters = {'srcbuffer','receive',...
                         'srcbufnum',1,... %TODO: handle in linking stage.
                         'srcframenum',0,...
                         'dstbuffer','none'};


    %% SeqControl
    vSeqControl = VSXSeqControl();
    vSeqControl = repmat(vSeqControl,1,6);
    vSeqControl = sort(copy(vSeqControl));
    
    vSeqControl(1) = setfield(setfield(vSeqControl(1), 'command', 'jump'), 'argument', 1);
    vSeqControl(2) = setfield(setfield(vSeqControl(2), 'command', 'timeToNextAcq'), 'argument', 160);
    vSeqControl(3) = setfield(setfield(vSeqControl(3), 'command', 'timeToNextAcq'), 'argument', 19040);
    vSeqControl(4) = setfield(vSeqControl(4), 'command', 'returnToMatlab');
    vSeqControl(5) = setfield(vSeqControl(5), 'command', 'transferToHost');
    vSeqControl(6) = setfield(setfield(vSeqControl(6), 'command', 'noop'), 'argument', 100/0.2);% 'condition', 'Hw&Sw');
    
    jump_to_image_start = vSeqControl(1);
    wait_for_tx_pulse = vSeqControl(2);
    wait_for_pulse_sequence = vSeqControl(3);
    return_to_matlab = vSeqControl(4);
    transfer_to_host = vSeqControl(5);
    no_operation = vSeqControl(6);
    
    
    

    %% Event

    % loop through all events
    % ---------- Events ------------- %
    vEvent = copy(repmat(VSXEvent('seqControl', wait_for_tx_pulse), [1 us.seq.numPulse]));
    vEvent = sort(vEvent);
    
    for i = 1:us.seq.numPulse % each transmit
        vEvent(i).info = 'Full aperture.';
        vEvent(i).tx  = vTX(i);
        vEvent(i).rcv = vRcv(i);
        vEvent(i).rcv.acqNum = i;
    end
    
    % transfer data to host using the last event
    vEvent(i).seqControl = [wait_for_pulse_sequence, transfer_to_host]; % modify last acquisition vEvent's seqControl
    
    % post-processing events and return to MATLAB   
    vEvent(end+1) = VSXEvent(...
        'info', 'recon and process', ...
        'recon', vRecon, ...
        'process', display_image_process, ...        
        'seqControl', return_to_matlab ...        
    );
    
    % save RF Data
    vEvent(end+1) = VSXEvent(...
        'info', 'Save RF Data', ...
        'process', save_rf_data,...
        'seqControl', no_operation...
        );

    % return to start of block
    vEvent(end+1) = VSXEvent(...
        'info', 'Jump back',...
        'seqControl', [jump_to_image_start]);
    % ------------ Events ------------ %
    jump_to_image_start.argument = vEvent(1);
   
    %% ADDED UI
    vUI = VSXUI();
    vUI.Control =  {'UserB1','Style','VsPushButton','Label', 'SAVE RFData', 'Callback', @doRFDataStore};
    
    %% Block
    vBlock = VSXBlock();
    vBlock.vsxevent = vEvent;
        
    %% added External Functions/Callback 
    % restore warning state
    warning(warning_state);

end