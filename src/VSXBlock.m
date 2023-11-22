classdef VSXBlock < matlab.mixin.Copyable
    properties
        capture VSXEvent % data capture events
        post (1,:) VSXEvent = VSXEvent.empty % post-processing events
        next {mustBeScalarOrEmpty} = VSXEvent.empty % event to jump to after loop
        vUI VSXUI = VSXUI.empty % UI events
    end
    methods
        function obj = VSXBlock(kwargs)
            arguments, kwargs.?VSXBlock, end
            for f = string(fieldnames(kwargs))'
                obj.(f) = kwargs.(f);
            end
        end
        function vStruct = link(vblock, vResource)
            arguments
            vblock VSXBlock
            vResource (1,1) VSXResource = VSXResource(); % default resource
        end
            
            % squash obj to struct warning
            warning_state = warning('off', 'MATLAB:structOnObject');

            % identify which block each "next" value belongs to if any
            for i = numel(vblock):-1:1
                for j = numel(vblock):-1:1
                    if ismember(vblock(i).next, vblock(j).capture)
                        blkid{i} = j; break;
                    end
                end
            end

            % get the sequence of events, in order
            vEvent = arrayfun(@(vb, i) {[vb.capture(:); vb.post(:); jump2event(vb.next, i{1})]'}, vblock, blkid); % order as capture, post, jump2next(next_event)
            vEvent = [vEvent{:}]; 

            %% get arrays of all properties
            vTx             = unique([vEvent.tx]        , 'stable'); 
            vRcv            = unique([vEvent.rcv]       , 'stable');
            vRecon          = unique([vEvent.recon]     , 'stable');
            vReconInfo      = unique([vRecon.RINums]    , 'stable');
            vProcess        = unique([vEvent.process]   , 'stable');
            vSeqControl     = unique([vEvent.seqControl], 'stable');
            vTw             = unique([vTx.waveform]     , 'stable');
            vTGC            = unique([vRcv.TGC]         , 'stable');
            vUI             = unique([vblock.vUI]       , 'stable'); %#ok<PROPLC> 
            
            vDisplayWindow  = unique([vResource.DisplayWindow], 'stable');
            vRcvBuffer      = unique([vResource.RcvBuffer    ], 'stable');
            vInterBuffer    = unique([vResource.InterBuffer  ], 'stable');
            vImageBuffer    = unique([vResource.ImageBuffer  ], 'stable');
            vParameters     = unique([vResource.Parameters   ], 'stable');

            % PData can be referenced in the Parameters field of a Process 
            % (Image or Doppler specifically, if it exists) or in the Recon
            % structure
            vpd = cellfun(@(p) p(2*find("pdatanum" == p(1:2:end))), {vProcess.Parameters}, 'UniformOutput', false);
            vpd = [vpd{:}];
            vPData = unique([vRecon.pdatanum, vpd{:}], 'stable');
            
            %% convert to Verasonics definition
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
            UI          = arrayfun(@struct, vUI); %#ok<PROPLC> 
            
            % 
            Resource.DisplayWindow  = arrayfun(@struct, vDisplayWindow);
            Resource.RcvBuffer      = arrayfun(@struct, vRcvBuffer);
            Resource.InterBuffer    = arrayfun(@struct, vInterBuffer);
            Resource.ImageBuffer    = arrayfun(@struct, vImageBuffer);
            Resource.Parameters     = arrayfun(@struct, vParameters);            
            
            % remove TX.aperture property if unneeded / unspecified
            if all(cellfun(@isempty, {TX.aperture}))
                TX = rmfield(TX, 'aperture'); 
            end
            if all(cellfun(@(x)all(isnan(x)), {TX.FocalPt}))
                TX = rmfield(TX, 'FocalPt'); 
            end
            
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

            %% Parse special arguments
            % TW - Remove unused properties
            if all(cellfun(@isempty, {TW.VDASArbwave}))
                TW = rmfield(TW, 'VDASArbwave');
            end

            % SeqControl
            for i = 1:numel(SeqControl)
                % link the index for jump commands
                if SeqControl(i).command == "jump"
                   SeqControl(i).argument = safeIsMember([SeqControl(i).argument], vEvent);  
                end
                
                % 0 arguments -> empty TODO: use nan to mark empty arguments
                if SeqControl(i).argument == 0
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
            nms  = {'Event', 'TX', 'Receive', 'Recon', 'ReconInfo', 'Process', 'SeqControl', 'TW', 'Resource', 'PData', 'TGC', 'UI'};
            vals = {Event, TX, Receive, Recon, ReconInfo, Process, SeqControl, TW, Resource, PData, TGC, UI};
            args = [nms; vals];
            vStruct = struct(args{:});
           
            % convert some logicals to double
	    [vStruct.Receive.callMediaFunc] = dealfun(@double, vStruct.Receive.callMediaFunc);
	    if ~isempty(ReconInfo), [vStruct.ReconInfo.threadSync ] = dealfun(@double, vStruct.ReconInfo.threadSync ); end
	    [vStruct.TW.sysExtendBL       ] = dealfun(@double, vStruct.TW.sysExtendBL       );
	    [vStruct.TW.perChWvfm         ] = dealfun(@double, vStruct.TW.perChWvfm         );

	    % remove empty structures
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

            % restore warning state
            warning(warning_state);
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
arguments, vEvent {mustBeScalarOrEmpty}, i (1,:) double = []; end
if isscalar(vEvent)
vEvent = VSXEvent('info', "Jump to " + vEvent.info + (" of block " + i), 'seqControl', ...
    VSXSeqControl('command', 'jump', 'argument', vEvent) ...
    );
end
end
