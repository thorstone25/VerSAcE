classdef VSXBlock < matlab.mixin.Copyable
    properties
        capture VSXEvent % data capture events
        post (1,:) VSXEvent = VSXEvent.empty % post-processing events
        next VSXEvent {mustBeScalarOrEmpty} = VSXEvent.empty % event to jump to after loop
        vUI VSXUI = VSXUI.empty % UI events
        vTPC VSXTPC = VSXTPC("name","Default"); % Default TPC (no arguments)
    end
    methods
        function obj = VSXBlock(kwargs)
            arguments, kwargs.?VSXBlock, end
            for f = string(fieldnames(kwargs))'
                obj.(f) = kwargs.(f);
            end
        end
        function vStruct = link(vblock, vResource, Trans, kwargs)
            % link - Index and pre-process the VSXBlock array
            %
            % vStruct = link(vblock, vResource) returns a struct vStruct
            % that can be saved to a mat-file and used with Vantage
            % software given the VSXBlock bvlock and the VSXResource
            % vResource. 
            %
            % vStruct = link(vblock, vResource, Trans) additionally appends
            % the Trans property to vStruct and adds an `aperture` field to
            % `TX` and `Receive` if appropriate.
            %
            % vStruct = link(..., 'TXPD', true) additionally runs
            % `computeTXPD` to compute the transmit power density. The
            % default is true if a `VSXRecon` is present.
            % 
            % Example:
            % % Choose a Transducer
            % Trans = struct('name', 'L12-3v', 'units', 'mm');
            % Trans = computeTrans(Trans);
            % 
            % % Create a system definition in QUPS
            % xdc = Transducer.Verasonics(Trans);
            % seq = SequenceRadial('type','PW','angles',-10:10,'c0',1500);
            % scan = ScanCartesian('x',(-20:1/8:20)*1e-3,'z',(0:1/8:40)*1e-3);
            % us = UltrasoundSystem('xdc', xdc, 'seq', seq, 'scan', scan);
            % 
            % % Convert to a VSXBlock
            % vres = VSXResource(); % system-wide resource
            % vBlock = QUPS2VSX(us, Trans, vres); % convert to VSX objects
            %
            % % Link and pre-process
            % vStruct = vBlock.link(vres, Trans);
            %
            % % Save and run
            % filename = 'qups-vsx.mat';
            % save(filename, '-struct', 'vStruct');
            % VSX;
            %
            % See also: QUPS2VSX

            arguments
            vblock VSXBlock
            vResource (1,1) VSXResource = VSXResource(); % default resource
            Trans {mustBeScalarOrEmpty} = struct.empty % transducer
            kwargs.TXPD (1,1) logical = ~isempty([[vblock.capture.recon], [vblock.post.recon]]) % whether to parse TXPD
        end
            
            % identify which block each "next" value belongs to if any
                % no matching 'other' block - this is independent!
            blkid = num2cell(zeros(1,numel(vblock)));
            for i = numel(vblock):-1:1
                for j = numel(vblock):-1:1
                    if ismember(vblock(i).next, vblock(j).capture)
                        blkid{i} = j; break;
                    end
                end
            end

            % get the sequence of events, in order
            vEvent = arrayfun(@(vb, i) {[vb.capture(:); vb.post(:); jump2event(vb.next, i{1})]'}, vblock, blkid); % order as capture, post, jump2next(next_event)
            vEvent = unique([vEvent{:}], 'stable');

            %% get arrays of all properties
            [vTx,  itx]     = unique([vEvent.tx]        , 'stable'); 
            [vRcv, irx]     = unique([vEvent.rcv]       , 'stable');
            vRecon          = unique([vEvent.recon]     , 'stable');
            vReconInfo      = unique([vRecon.RINums]    , 'stable');
            vProcess        = unique([vEvent.process]   , 'stable');
            vSeqControl     = unique([vEvent.seqControl], 'stable');
            vTw             = unique([vTx.waveform]     , 'stable');
            vTGC            = unique([vRcv.TGC]         , 'stable');
            vUI             = unique([vblock.vUI]       , 'stable'); %#ok<PROPLC> 
            vTPC            = unique([vblock.vTPC]      , 'stable'); %#ok<PROPLC> 
            
            % Processes TODO: extract displays, buffers from here.
            iImg = (arrayfun(@(v)v.classname == "Image", vProcess)); % image processes
            iDop = (arrayfun(@(v)v.classname == "Doppler", vProcess)); % Doppler processes
            iExt = (arrayfun(@(v)v.classname == "External", vProcess)); % Doppler processes

            % get from the resource list
            vDisplayWindow  = unique([vResource.DisplayWindow], 'stable');
            vRcvBuffer      = unique([vResource.RcvBuffer    ], 'stable');
            vInterBuffer    = unique([vResource.InterBuffer  ], 'stable');
            vImageBuffer    = unique([vResource.ImageBuffer  ], 'stable');
            vParameters     = unique([vResource.Parameters   ], 'stable');

            % PData can be referenced in the Parameters field of a Process 
            % (Image or Doppler specifically, if it exists) or in the Recon
            % structure
            vpd = cellfun(@(p) p(2*find("pdatanum" == p(1:2:end))), {vProcess(iImg | iExt).Parameters}, 'UniformOutput', false);
            if ~isempty(vpd), vpd = [vpd{:}]; end
            vPData = unique([vRecon.pdatanum, vpd{:}], 'stable');

            % HACK: Vantage utilities expect `Receive`s to be sorted by buffer, frame, and acquisitions(?) 
            % even if referred to out of order in terms of `Event`s
            [~, i] = ismember([vRcv.bufnum], vRcvBuffer); % buffer indexing
            [~, i] = sortrows([i; [vRcv.framenum]; [vRcv.acqNum];]'); % acquisition/frame/buffer sort order
            vRcv = vRcv(i); % re-order

            %% assign aperture indices if possible
            if ~isempty(Trans) && isfield(Trans, 'HVMux')
                % set TX/RX apertures from available apertures
                % TODO: include generation of new (compliant) apertures
                aps = Trans.HVMux.ApertureES; % list of apertures

                % get tx aperture indices
                out = cellfun(@(apd) {apSelect(apd, aps, "buffer")}, {vTx.Apod}); % choose first match active aperture
                eind = cellfun(@isempty, out); % empty indices
                if any(eind), error("VSXBlock:apertureUnavailable","Unable to find any aperture to satisfy transmit " ...
                        + join(cellstr({vEvent(itx(eind)).info}), ", ") + "."); end
                [vTx.aperture] = deal(out{:});

                % get rx aperture indices
                out = cellfun(@(apd) {apSelect(apd, aps, "buffer")}, {vRcv.Apod}); % choose first match active aperture
                eind = cellfun(@isempty, out); % empty indices
                if any(eind), error("VSXBlock:apertureUnavailble","Unable to find any aperture to satisfy receive " ...
                        + join(string(vEvent(irx(eind)).info), ", ") + "."); end
                [vRcv.aperture] = deal(out{:});
            end

            %% convert to Verasonics definition
            % squash obj to struct warning
            warning_state = warning('off', 'MATLAB:structOnObject');

            % convert to structs
            Event       = arrayfun(@struct, vEvent);
            TX          = arrayfun(@struct, vTx);
            Receive     = arrayfun(@struct, vRcv);
            Recon       = arrayfun(@struct, vRecon);
            ReconInfo   = arrayfun(@struct, vReconInfo);
            Process     = arrayfun(@struct, vProcess);
            SeqControl  = arrayfun(@struct, vSeqControl);     
            TW          = arrayfun(@struct, vTw);
            Resource    = arrayfun(@struct, vResource);
            PData       = arrayfun(@struct, vPData);
            TGC         = arrayfun(@struct, vTGC);
            TPC         = arrayfun(@struct, vTPC); %#ok<PROPLC> 
            UI          = arrayfun(@struct, vUI); %#ok<PROPLC> 
            
            % 
            Resource.DisplayWindow  = arrayfun(@struct, vDisplayWindow);
            Resource.RcvBuffer      = arrayfun(@struct, vRcvBuffer);
            Resource.InterBuffer    = arrayfun(@struct, vInterBuffer);
            Resource.ImageBuffer    = arrayfun(@struct, vImageBuffer);
            Resource.Parameters     = arrayfun(@struct, vParameters);            
            
            % restore warning state
            warning(warning_state);

            % remove empty properties: Receive.aperture, TX.aperture, 
            % TX.FocalPt ReconInfo.LUT, Recon.rcvLUT
            for f = ["Origin","Steer","focus","FocalPt","TXPD","aperture"]
                TX =  rmifempty(TX, f);
            end
            Receive     = rmifempty(Receive, 'aperture');
            Recon       = rmifempty(Recon, 'rcvLUT');
            ReconInfo   = rmifempty(ReconInfo, 'LUT');
            ReconInfo   = rmifempty(ReconInfo, 'pagenum');
            Resource    = rmifempty(Resource, 'RcvBuffer');
            Resource    = rmifempty(Resource, 'InterBuffer');
            Resource    = rmifempty(Resource, 'ImageBuffer');
            Resource    = rmifempty(Resource, 'DisplayWindow');
            
            % remove all empty properties of TPC
            for f = string(fieldnames(TPC))', TPC = rmifempty(TPC, f); end

            %% assign indices
            % assign indices for Event
            for i = 1:numel(Event)
                % find indices and set empty indices to 0
                Event(i).tx         = safeIsMember([vEvent(i).tx]        , vTx);
                Event(i).rcv        = safeIsMember([vEvent(i).rcv]       , vRcv);
                Event(i).recon      = safeIsMember([vEvent(i).recon]     , vRecon);
                Event(i).process    = safeIsMember([vEvent(i).process]   , vProcess);
                Event(i).seqControl = safeIsMember([vEvent(i).seqControl], vSeqControl);
            end
            
            % assign indices for TX
            for i = 1 : numel(TX)
                TX(i).waveform = safeIsMember([vTx(i).waveform], vTw);
            end
            
            % assign indices for Receive
            for i = 1 : numel(Receive)
                Receive(i).bufnum = safeIsMember([vRcv(i).bufnum], vRcvBuffer);
                Receive(i).TGC    = safeIsMember([vRcv(i).TGC   ], vTGC);
            end
            
            % assign indices for Recon
            for i = 1 : numel(Recon)
                Recon(i).RINums   = safeIsMember([vRecon(i).RINums], vReconInfo);
                Recon(i).pdatanum = safeIsMember([vRecon(i).pdatanum], vPData);
                Recon(i).IntBufDest = [safeIsMember([vRecon(i).IntBufDest], vInterBuffer), vRecon(i).IntBufDestFrm];
                Recon(i).ImgBufDest = [safeIsMember([vRecon(i).ImgBufDest], vImageBuffer), vRecon(i).ImgBufDestFrm];
                if ~Recon(i).IntBufDest(1), Recon(i).IntBufDest = []; end % remove if no buffer specified
                if ~Recon(i).ImgBufDest(1), Recon(i).ImgBufDest = []; end % remove if no buffer specified
            end

            % remove extra *Frm argument
            if ~isempty(Recon), Recon = rmfield(Recon, ["IntBufDestFrm", "ImgBufDestFrm"]); end

            % assign indices for ReconInfo
            for i = 1:numel(vReconInfo)
                ReconInfo(i).txnum  = safeIsMember([vReconInfo(i).txnum ], vTx );
                ReconInfo(i).rcvnum = safeIsMember([vReconInfo(i).rcvnum], vRcv);
                % TODO: add regionnum property vs. VSXRegion
                if isempty(ReconInfo(i).regionnum)
                    ReconInfo(i).regionnum = 1; % if no regions specified, use 1st region
                end
            end

            % For UI: allow an 'auto' option for 'button positions' arguments
            pos = "User" + ["A","B","C"] + flip(1:8)';
            prf = ["UserB"+(1:7), "UserC"+(6:8), "UserA"+(7:8), "UserC"+(4:8), "UserB8"];
            pos = union(prf, pos, 'stable');
            for i = 1:numel(UI)
                upos = arrayfun(@(UI) string(UI.Control{1}), UI); % positions in use
                fpos = setdiff(pos, upos, 'stable');
                switch upos(i)
                    case "auto", UI(i).Control{1} = char(fpos(1)); % first free position
                end
            end

            %% Parse special arguments
            % TW - Remove unused properties
            if all(cellfun(@isempty, {TW.VDASArbwave}))
                TW = rmfield(TW, 'VDASArbwave');
                for i = 1 : numel(TW)
                    if TW(i).type ~= "envelope"
                       [TW(i).envNumCycles, TW(i).envFrequency, TW(i).envPulseWidth]  = deal([]); % remove envelope properties
                    end
                    if TW(i).type ~= "parametric"
                        [TW(i).Parameters, TW(i).equalize] = deal([]); % remove parametric properties
                    end
                    if TW(i).type ~= "sampled"
                        TW(i).Waveform = [];
                    end
                    if TW(i).type ~= "sampled" && TW(i).type ~= "function"
                        [TW(i).frequency, TW(i).Descriptors] = deal([]); % remove function/sampled properties
                    end
                end
            end

            % SeqControl
            for i = 1:numel(SeqControl)
                % link the index for jump commands
                if  SeqControl(i).command == "jump"
                    SeqControl(i).argument = safeIsMember([SeqControl(i).argument], vEvent);
                end

                % link the index for setTPCProfile commands
                if  SeqControl(i).command == "setTPCProfile"
                    SeqControl(i).argument = safeIsMember([SeqControl(i).argument], vTPC); %#ok<PROPLC>
                end

                % 0 arguments -> empty TODO: use nan to mark empty arguments
                if isnan(SeqControl(i).argument)
                    SeqControl(i).argument = [];
                end
            end

            % Parse Process arguments
            for i = 1:numel(Process)
                % turn the params into a struct for easier processing
                prms = struct(Process(i).Parameters{:});
                
                % for all referenced buffers or display windows
                for f = ["displayWindow", "pdatanum", "imgbufnum", "srcbufnum", "dstbufnum"]
                    % if it's specified
                    if isfield(prms, f)
                        % get the corresponding property
                        switch f
                            case "displayWindow", arr = vDisplayWindow;
                            case "pdatanum"     , arr = vPData;
                            case "imgbufnum"    , arr = vImageBuffer;
                            case {"srcbufnum", "dstbufnum"}
                                % get the property name of the type of buffer
                                if     f == "srcbufnum", sb = "srcbuffer";
                                elseif f == "dstbufnum", sb = "dstbuffer";
                                else, error("Unrecognized buffer type.");
                                end

                                % get the type of buffer
                                switch prms.(sb)
                                    case "receive",             arr = vRcvBuffer;
                                    case "inter",               arr = vInterBuffer;
                                    case {'image', 'imageP'},   arr = vImageBuffer;
                                    case 'none',                continue;
                                    otherwise
                                        error("Unrecognized " + sb + " value for Process " + i + " of class " + Process(i).classname + ".");
                                end
                        end
                        % get argument index
                        j = 2*find(f == string(Process(i).Parameters(1:2:end)));

                        % link
                        if ~isempty(j), Process(i).Parameters{j} = safeIsMember(Process(i).Parameters{j}, arr); end
                    end
                end
            end

            %% Casting
            % convert to a single struct
            nms  = {'Event', 'TX', 'Receive', 'Recon', 'ReconInfo', 'Process', 'SeqControl', 'TW', 'Resource', 'PData', 'TGC', 'TPC', 'UI'};
            vals = { Event ,  TX ,  Receive ,  Recon ,  ReconInfo ,  Process ,  SeqControl ,  TW ,  Resource ,  PData ,  TGC ,  TPC ,  UI };
            args = [nms; vals];
            vStruct = struct(args{:});

            % add Trans
            if ~isempty(Trans), vStruct.Trans = Trans; end

            % convert some logicals to double
            [vStruct.Receive.callMediaFunc] = dealfun(@double, vStruct.Receive.callMediaFunc);
            if ~isempty(ReconInfo), [vStruct.ReconInfo.threadSync ] = dealfun(@double, vStruct.ReconInfo.threadSync ); end
            [vStruct.TW.sysExtendBL       ] = dealfun(@double, vStruct.TW.sysExtendBL       );
            [vStruct.TW.perChWvfm         ] = dealfun(@double, vStruct.TW.perChWvfm         );
            for i = 1:numel(vStruct.Process)
                j = vStruct.Process(i).Parameters(1:2:end) == "display";
                k = 2*find(j);
                if any(j), vStruct.Process(i).Parameters{k} = double(vStruct.Process(i).Parameters{k}); end
            end

            % remove empty structures, or structures of all empty
            % properties
            for f = string(fieldnames(vStruct))'
                if isempty(vStruct.(f))
                    vStruct = rmfield(vStruct, f);
                end

            end

            % convert all strings to char
            vStruct = recursiveString2Char(vStruct);

            % clear nan values from TX struct
            for f = ["focus", "Steer", "FocalPt"]
                if isfield(vStruct.TX, f)
                    for i = 1:numel(vStruct.TX)
                        if any(isnan(vStruct.TX(i).(f)))
                            vStruct.TX(i).(f) = []; % replace with empty
                        end
                    end
                end
            end

            %% Compute
            % ... once we have all native structs and chars (not classes
            % and strings), we can run native Vantage functions, although
            % we'll need to hijack the base workspacae to stay compliant.
            
            % Compute transmit power density for corresponding Recons
            if kwargs.TXPD
                % assert the local vars exist
                % nms = ["TW", "Trans", "Resource"]; % vars silently required?
                nms = string(nms(1:end-2)); % is literally everything required???
                exg = find(arrayfun(@(n) evalin('base', "exist('"+n+"','var');"), nms)); % exists in global 
                exl = arrayfun(@(n) isfield(vStruct, n), nms); % exist in local
                if ~all(exl)
                    error("Fields " + join(nms, " and ") ...
                        + " are required for TXPD computation."); 
                end

                % swap local and global (base) workspace
                vls = struct(); for f = nms(exg), vls.(f) = evalin('base', f); end % store global vars in local
                for f = nms(exl), assignin("base", f, vStruct.(f)); end % store local vars in global

                % compute transmit power density in base workspace
                for i = 1:numel(Recon) % for each reconstruction process
                    for j = unique([ReconInfo(Recon(i).RINums).txnum]) % for each TX index
                        assignin('base', 'j', j); assignin('base', 'i', i);
                        evalin('base', "TX(j).TXPD = computeTXPD(TX(j), PData(Recon(i).pdatanum));") % requires global vars
                    end
                end

                % restore swap local/global workspace vars
                for f = nms(exl), vStruct.(f) = evalin("base", f); end % copy global to local 
                evalin("base", join(["clearvars", nms(exl)], " ")); % delete all local vars temp. copied to base
                for f = nms(exg), assignin("base", f, vls.(f)); end % return local vars to global
            end
        end
    end
    
    % validation
    methods(Static)
        function [TX, PData, Trans] = validateTXPD(TX, PData, Trans)
            
            % showTXPD visualize transmit TXPD data with cutoffs
            %   Call with no arguments to process all TX structures with first PData struct.
            %   If a range of TX values is desired, set the first varargin to the range, e.g. [1:10].
            %   If a PData struct other than one is to be specified, set the range of TX values and
            %      then specify the PData number in the 2nd argument.
            %
            % Last Update: 1-10-2017
                        
            PDataBase = evalin('base','PData');
            PDsldrMax = length(PDataBase);
            PDsldrMin = 0.99;
            if length(PDataBase) > 1
                PDsldrSS  = 1/max(1,length(PDataBase)-1);
            else
                PDsldrSS = 1;
            end
            
            Trans = evalin('base','Trans');
            
            % default value (indices)
            TXIn = 1:size(TX,2);
            PDIn = 1;
            
            pdelta=[]; pdeltaX=[]; pdeltaY=[]; pdeltaZ=[]; pdeltaR=[]; pdeltaT=[];
            %% Get the default value, PDIn is pdata num, txn is tx num
            txn = TXIn(1);
            PData = PDataBase(PDIn);
            
            % Check for presence of the TXPD data, and compute if not already computed.
            if ~isfield(TX,'TXPD')
                h = waitbar(0,'Calculate TPXD, please wait!');
                for i=1:size(TXIn,2)
                    TX(TXIn(i)).TXPD = computeTXPD(TX(TXIn(i)),PData);
                    waitbar(i/size(TXIn,2))
                end
                [TX.peakCutOff] = deal(1.0);
                [TX.peakBLMax] = deal(4.0);
                close(h)
            else
                for i=1:size(TXIn,2)
                    if isempty(TX(TXIn(i)).TXPD)
                        TX(TXIn(i)).TXPD = computeTXPD(TX(TXIn(i)),PData);
                        TX(TXIn(i)).peakCutOff = 1.0;
                        TX(TXIn(i)).peakBLMax = 4.0;
                    end
                end
            end
            
            % update the pdelta setting for showTXPD plot
            Delta = updatePData(PData);
            deltaFieldsName = fieldnames(Delta);
            for i = 1:length(deltaFieldsName)
                eval(['pdelta',deltaFieldsName{i},'= Delta.',deltaFieldsName{i},';']);
            end
            Nrow = PData.Size(1); Ncol = PData.Size(2); Nsec = PData.Size(3);
            
            function Delta = updatePData(PData)
                
                if isfield(PData,'PDelta') && ~isempty(PData.PDelta)
                    if isfield(PData,'Coord')
                        switch PData.Coord
                            case 'rectangular'
                                if size(PData.PDelta,2) ~= 3
                                    error('computeRegions: PData(x).PDelta number of columns must be 3.');
                                end
                                Delta.X = PData.PDelta(1);
                                Delta.Y = PData.PDelta(2);
                                Delta.Z = PData.PDelta(3);
                            case 'polar'
                                if size(PData.PDelta,2) ~= 3
                                    error('computeRegions: PData(x).PDelta number of columns must be 3.');
                                end
                                Delta.T = PData.PDelta(1);
                                Delta.R = PData.PDelta(2);
                                Delta.Z = PData.PDelta(3);
                            case 'cylindrical'
                                error('computeRegions: PData(x).Coor = ''cylindrical'' is not yet supported.');
                            otherwise
                                error('computeRegions: Unrecognized PData(x).Coor string.');
                        end
                    else
                        % Default to rectangular coordinates if no Coord attribute
                        if size(PData.PDelta,2) ~= 3
                            error('computeRegions: PData(x).PDelta number of columns must be 3.');
                        end
                        Delta.X = PData.PDelta(1);
                        Delta.Y = PData.PDelta(2);
                        Delta.Z = PData.PDelta(3);
                    end
                elseif isfield(PData,'pdelta') && ~isempty(PData.pdelta)
                    % if no PDelta and PData.pdelta is set, it implies all non-specified pdelta? are equal to its value.
                    Delta.pdelta = PData.pdelta;
                    Delta.X = pdelta;
                    Delta.Y = pdelta;
                    Delta.Z = pdelta;
                    Delta.R = pdelta;
                    Delta.T = pdelta/64;  % pdeltaT is in radians, for example if pdelta=0.5, pdeltaT = .0078 or
                    % about 200 radial lines in 90 degrees
                else
                    % if PDelta and pdelta are missing, at the least, pdeltaX and pdeltaZ must be specified.
                    if ~isfield(PData,'pdeltaX') || isempty(PData.pdeltaX)
                        error('computeRegions: PData(x).pdeltaX must be specified if PData(x).PDelta or pdelta missing.');
                    elseif ~isfield(PData,'pdeltaZ') || isempty(PData.pdeltaZ)
                        error('computeRegions: PData(x).pdeltaZ must be specified if PData(x).PDelta or pdelta missing.');
                    end
                    % Read specified pdeltaX,Y,Z,R or T.
                    Delta.X = PData.pdeltaX;
                    if isfield(PData,'pdeltaY')&&(~isempty(PData.pdeltaY)), Delta.Y = PData.pdeltaY; end
                    Delta.Z = PData.pdeltaZ;
                    if isfield(PData,'pdeltaR')&&(~isempty(PData.pdeltaR)), Delta.R = PData.pdeltaR; end
                    if isfield(PData,'pdeltaT')&&(~isempty(PData.pdeltaT)), Delta.T = PData.pdeltaT; end
                end
            end
            
            %%
            % Determine limits for axes. Default is for 2D display, will be changed by
            % 3D display
            if isfield(PData,'Coord')&&(strcmp(PData.Coord,'polar'))
                XLim = [pdeltaT*(1-PData.Size(2))/2, pdeltaT*(PData.Size(2)-1)/2];
                YLim = [0, pdeltaR * PData.Size(1)];
            else
                XLim = [PData.Origin(1), PData.Origin(1) + pdeltaX * PData.Size(2)];
                YLim = [PData.Origin(3), PData.Origin(3) + pdeltaZ * PData.Size(1)];
            end
            
            %% UIbutton group for plot selection.            
            % TX peak cutoff slider
            if ~isfield(TX,'peakCutOff')||isempty(TX(txn).peakCutOff), TX(txn).peakCutOff = 0; end
            
            if ~isfield(TX,'peakBLMax')||isempty(TX(txn).peakBLMax), TX(txn).peakBLMax = BLMax; end
            
            % Region number slider
            if ~isfield(PData,'Region')||isempty(PData.Region)
                evalin('base',['[PData(',int2str(PDIn),').Region] = computeRegions(PData(',int2str(PDIn),'));']);
                PData = evalin('base',['PData(',int2str(PDIn),')']);
            end
            regNum = 1;
           
            AxesUnit = 'wavelength';
            
            speedOfSound = 1540;
            if evalin('base','exist(''Resource'',''var'')')
                Resource = evalin('base','Resource');
                if isfield(Resource,'Parameters')
                    if isfield(Resource.Parameters,'speedOfSound')
                        speedOfSound = Resource.Parameters.speedOfSound/1000; % speed of sound in mm/usec
                    end
                end
                
                if isfield(Resource,'DisplayWindow')
                    if isfield(Resource.DisplayWindow,'Title')
                        hf.Name = ['showTXPD - ',Resource.DisplayWindow.Title];
                    end
                end
                
            else
                disp('Speed of sound is not defined, use 1540 m/s as the default value');
            end
                        
            % Calculate the maximum intensity of the TXPD data.
            maxIntensity = 0;
            for i = 1:size(TXIn,2)
                for j = 1:Nsec
                    TXPD = double(TX(TXIn(i)).TXPD(:,:,j,1));
                    peak = max(max(TXPD));
                    if peak>maxIntensity, maxIntensity = peak; end
                    clear TXPD
                end
            end
            maxIntensity = maxIntensity/256;
            plotCode = 'intensity';  % default plotCode
            [Data,clims] = updateImage;
            
            %% Nested functions
            function [Data,clims] = updateImage()
                Nrow = PData.Size(1); Ncol = PData.Size(2); Nsec = PData.Size(3);
                setYDir = 0;
                if Nsec > 1
                    switch shapeOrientation
                        case 'xz'
                            TXPD = zeros(Nsec,Ncol,3);
                            TXPD(:,:,1) = permute(double(TX(txn).TXPD(sectNum,:,:,1)),[3 2 1])/256;
                            TXPD(:,:,2) = permute(double(TX(txn).TXPD(sectNum,:,:,2)),[3 2 1])/16;
                            TXPD(:,:,3) = permute(double(TX(txn).TXPD(sectNum,:,:,3)),[3 2 1])/16;
                            XLim = [-PData.Origin(2), -PData.Origin(2) + pdeltaY * PData.Size(1)];
                            YLim = [PData.Origin(3), PData.Origin(3) + pdeltaZ * PData.Size(3)];
                            sliceLA = sub2ind(PData.Size, sectNum*(ones(Nsec,Ncol)),ones(Nsec,1)*(1:1:Ncol),(ones(Ncol,1)*(1:1:Nsec))');
                        case 'yz'
                            TXPD = zeros(Nsec,Nrow,3);
                            TXPD(:,:,1) = permute(double(TX(txn).TXPD(:,sectNum,:,1)),[3 1 2])/256;
                            TXPD(:,:,2) = permute(double(TX(txn).TXPD(:,sectNum,:,2)),[3 1 2])/16;
                            TXPD(:,:,3) = permute(double(TX(txn).TXPD(:,sectNum,:,3)),[3 1 2])/16;
                            XLim = [PData.Origin(1), PData.Origin(1) + pdeltaX * PData.Size(2)];
                            YLim = [PData.Origin(3), PData.Origin(3) + pdeltaZ * PData.Size(3)];
                            sliceLA = sub2ind(PData.Size, ones(Nsec,1)*(1:1:Nrow),sectNum*ones(Nsec,Nrow),(ones(Nrow,1)*(1:1:Nsec))');
                        case 'xy'
                            TXPD = zeros(Nrow,Ncol,3);
                            TXPD(:,:,1) = (double(TX(txn).TXPD(:,:,sectNum,1)))/256;
                            TXPD(:,:,2) = (double(TX(txn).TXPD(:,:,sectNum,2)))/16;
                            TXPD(:,:,3) = (double(TX(txn).TXPD(:,:,sectNum,3)))/16;
                            XLim = [PData.Origin(1), PData.Origin(1) + pdeltaX * PData.Size(2)];
                            YLim = [PData.Origin(2), PData.Origin(2) - pdeltaY * PData.Size(1)];
                            sliceLA = sub2ind(PData.Size, (ones(Ncol,1)*(1:1:Nrow))',ones(Nrow,1)*(1:1:Ncol),sectNum*ones(Nrow,Ncol));
                            setYDir = 1;
                    end
                else
                    TXPD = double(TX(txn).TXPD);
                    TXPD(:,:,1) = TXPD(:,:,1)/256;
                    TXPD(:,:,2:3) = TXPD(:,:,2:3)/16;
                    if isfield(PData,'Coord')&&(strcmp(PData.Coord,'polar'))
                        XLim = [pdeltaT*(1-PData.Size(2))/2, pdeltaT*(PData.Size(2)-1)/2];
                        YLim = [0, pdeltaR * PData.Size(1)];
                    else
                        XLim = [PData.Origin(1), PData.Origin(1) + pdeltaX * PData.Size(2)];
                        YLim = [PData.Origin(3), PData.Origin(3) + pdeltaZ * PData.Size(1)];
                    end
                    sliceLA = reshape(find(ones(Nrow,Ncol)),Nrow,Ncol);
                end
                
                P = (TXPD(:,:,1)>TX(txn).peakCutOff)&(TXPD(:,:,3)<TX(txn).peakBLMax);
                P = +P;  % convert logical to double.
                switch plotCode
                    case 'intensity'
                        Data = TXPD(:,:,1) .* P;
                        clims = [0 maxIntensity];
                    case 'peakTime'
                        Data = TXPD(:,:,2) .* P;
                        clims = [0 max(max(TXPD(:,:,2)))];
                    case 'burstLength'
                        Data = TXPD(:,:,3) .* P;
                        clims = [0 max(0.5,max(max(Data)))]; % don't allow max to go to 0
                end
                
                if isequal(clims(1),clims(2)), clims(2) = clims(1)+0.01; end % prevent error message from caxis
                
                if (get(showReg,'Value') == 1)&&(regNum~=0)
                    regNum = round(get(nregSldr,'Value'));
                    if ~isfield(PData,'Region')||isempty(PData.Region), return, end
                    if ~isfield(PData.Region(regNum),'numPixels')||isempty(PData.Region(regNum).numPixels), return, end
                    ci = clims(2)/4; % intensity of Region pixels
                    [row,col] = find(ismember(sliceLA,PData.Region(regNum).PixelsLA+1));
                    Data(sub2ind(size(Data),row(iseven(row+col)),col(iseven(row+col)))) = ci;
                end
                
                axes(ha),
                if strcmp(AxesUnit,'mm')
                    XLim = XLim * speedOfSound/Trans.frequency;
                    YLim = YLim * speedOfSound/Trans.frequency;
                end
                
                imagesc(XLim,YLim,Data,clims); if setYDir, set(ha,'YDir','normal'); end
                if isfield(PData,'Coord')&&(strcmp(PData.Coord,'polar'))
                    xlabel('Angle in Radians');
                    ylabel('Depth in Wavelengts');
                    set(changeUnit,'Visible','off');
                end
                
                % plot Transducer elements
                if Trans.type ~= 2
                    scaleToWvl = Trans.frequency/speedOfSound; % conversion factor from mm to wavelengths
                    
                    if strcmp(Trans.units, 'mm') && strcmp(AxesUnit,'wavelength')
                        ElementPos = Trans.ElementPos * scaleToWvl;
                    elseif strcmp(Trans.units, 'wavelengths') && strcmp(AxesUnit,'mm')
                        ElementPos = Trans.ElementPos / scaleToWvl;
                    else
                        ElementPos = Trans.ElementPos;
                    end
                    if isfield(PData,'Coord')&&(strcmp(PData.Coord,'polar'))
                        ElementPos(:,1) = linspace(XLim(1),XLim(2),Trans.numelements);
                        switch Trans.type
                            case 0 % phase array
                                if isfield(PData.Region(regNum).Shape,'Position')
                                    z = PData.Region(regNum).Shape.Position(3);
                                else
                                    z = PData.Origin(3); % entire PData
                                end
                                ElementPos(:,3) = abs(z./cos(ElementPos(:,1)));
                            case 1 % curved array
                                ElementPos(:,3) = Trans.radius;
                        end
                    end
                    hold on;
                    if isfield(Trans,'HVMux')
                        if isfield(TX(1,txn),'aperture') == 1
                            aperture = TX(1,txn).aperture;  %Import aperture Data
                        else
                            aperture = 1;
                        end
                        if isfield(Trans.HVMux,'ApertureES')
                            ActEle = Trans.HVMux.ApertureES(:,aperture);
                        else
                            ActEle = Trans.HVMux.Aperture(:,aperture);
                        end
                        ActEle(ActEle ~= 0) = 1;
                        ActEle(ActEle == 0) = NaN;
                    else
                        ActEle = ones(1,Trans.numelements);
                    end
                    plot(ElementPos(ActEle==1,1),ElementPos(ActEle==1,3),'rs');
                    hold off
                end
            end
                        
        end
        
    end
end

function X = safeIsMember(A, B)
   [~, X] = ismember(A, B);
   if isempty(X)
       X = 0;
   end
end

function y = recursiveString2Char(x)
y = x; % init
if isstruct(x) % struct container
    for i = 1:numel(x) % for all elements
        for f = string(fieldnames(x))' % for all fields
            y(i).(f) = recursiveString2Char(x(i).(f)); % do recursive string2char on all props of x
        end
    end
elseif iscell(x) % cell container
    for i = 1:numel(x) % for all elements
         y{i} = recursiveString2Char(x{i}); % do recursive string2char on all props of x
    end
elseif isstring(x)
    y = char(x);
    if xor(isempty(x), isempty(y))
        error("VSXOOD:VSXBlock:string2charEmptiness", "String and char emptiness not identical - use string.empty not """" for an empty string.");
    end
else
    % do nothing
end
end

function vEvent = jump2event(vEvent, i)
arguments, vEvent {mustBeScalarOrEmpty}, i double {mustBeScalarOrEmpty} = []; end
if isscalar(vEvent)
    if i ~= 0 % 0 -> no matching block -> vEvent is the next Event!
        vEvent = VSXEvent('info', "Jump to " + vEvent.info + (" of block " + i), 'seqControl', ...
            VSXSeqControl('command', 'jump', 'argument', vEvent) ...
            );
    end
end
end

function s = rmifempty(s, f)
% remove field if empty
%
% s = rmifempty(s, f) removes field f from struct s if f is a field of s
% and s(i).f is empty for all i
% 
arguments
    s struct
    f (1,1) string
end

% remove s.(f) property if unneeded / unspecified
if isfield(s, f) && all(cellfun(@isempty, {s.(f)}))
    s = rmfield(s, f);
end


end

function i = apSelect(apd, aps, mode)
arguments
    apd (:,1) {mustBeNumericOrLogical} % elems
    aps (:,:) {mustBeNumericOrLogical} % elems x apertures
    mode (1,1) string {mustBeMember(mode, ["first", "last", "median", "buffer"])} = "buffer"
end
val = all((~apd) | aps, 1); % elems off or aperture active -> valid aperture
if isempty(val), i = []; return; end % short-circuit on invalid apertures
switch mode
    case {'first', 'last'}
        i = find(val, 1, mode);
    case {'median'}
        i = find(val);
        i = i(ceil(numel(i)/2)); % center aperture
    case "buffer"
        caps = num2cell(aps(:,val),1); % place each aperture in a cell
        capd = num2cell(aps(:,val) & apd,1); % same with overlap

        % get start and end of aperture/apodization
        saps = cellfun(@(x)find(x,1,'first'), caps);
        eaps = cellfun(@(x)find(x,1,'last' ), caps);
        sapd = cellfun(@(x)find(x,1,'first'), capd);
        eapd = cellfun(@(x)find(x,1,'last' ), capd);

        % get max buffer index
        i = argmax(min(abs(saps - sapd),abs(eaps - eapd))); % max the min dist from aperture start/end to apodization start/end
        j = find(val); 
        i = j(i);
end
end
