classdef VSXBlock < matlab.mixin.Copyable
    properties
        vsxevent (1,:) VSXEvent
    end
    methods
        function vStruct = link(self, vResource, vPData, vUI)
        arguments
            self VSXBlock
            vResource
            vPData
            vUI VSXUI = VSXUI.empty
        end
            
            % squash obj to struct warning
            warning_state = warning('off', 'MATLAB:structOnObject');
            

            %% convert string to char
 %{
           for i = 1:numel(Event)
                for f1 = string(fieldnames(Event))'
                    for f2 = string(fieldnames(Event.(f1)))'
                        Event(i).(f1).(f2) = stringToChar(Event(i).(f1).(f2));
                    end
                end
            end
            for i = 1:numel(vPData)
                for f = string(fieldnames(vPData))'
                    vPData(i).(f)= stringToChar(vPData(i).(f));
                end
            end
            for i = 1:numel(vDisplayWindow)
                for f = string(fieldnames(vDisplayWindow))'
                    vDisplayWindow(i).(f)= stringToChar(vDisplayWindow(i).(f));
                end
            end
            for i = 1:numel(vRcvBuffer)
                for f = string(fieldnames(vRcvBuffer))'
                    vRcvBuffer(i).(f)= stringToChar(vRcvBuffer(i).(f));
                end
            end
            for i = 1:numel(vInterBuffer)
                for f = string(fieldnames(vInterBuffer))'
                    vInterBuffer(i).(f)= stringToChar(vInterBuffer(i).(f));
                end
            end
            for i = 1:numel(vImageBuffer)
                for f = string(fieldnames(vImageBuffer))'
                    vImageBuffer(i).(f)= stringToChar(vImageBuffer(i).(f));
                end
            end
            for i = 1:numel(vParameters)
                for f = string(fieldnames(vParameters))'
                    vParameters(i).(f)= stringToChar(vParameters(i).(f));
                end
            end
            for i = 1:numel(vRcv)
                for f = string(fieldnames(vRcv))'
                    vRcv(i).(f)= stringToChar(vRcv(i).(f));
                end
            end
%}
            %% get arrays of all properties
            vEvent              = self.vsxevent; % all events
            vTx                 = unique([vEvent.tx]); % array of all transmits
            vRcv                = unique([vEvent.rcv]);
            vRecon              = unique([vEvent.recon]);
            vReconInfo          = unique([vRecon.RINums]);
            vProcess            = unique([vEvent.process]);
            vSeqControl         = unique([vEvent.seqControl]);
            vTw                 = unique([vTx.waveform]);
            vTGC                = unique([vRcv.TGC]);
            
            % is this necessary?
            vDisplayWindow = unique([vResource.DisplayWindow]);
            vRcvBuffer     = unique([vResource.RcvBuffer    ]);
            vInterBuffer   = unique([vResource.InterBuffer  ]);
            vImageBuffer   = unique([vResource.ImageBuffer  ]);
            vParameters    = unique([vResource.Parameters   ]);
            
            %% convert to Verasonics definition
            % convert to structs
            Event = arrayfun(@struct, vEvent);
            TX = arrayfun(@struct, vTx);
            Receive = arrayfun(@struct, vRcv);
            Recon = arrayfun(@struct, vRecon);
            ReconInfo = arrayfun(@struct, vReconInfo);
            Process = arrayfun(@struct, vProcess);
            SeqControl = arrayfun(@struct, vSeqControl);     
            TW = arrayfun(@struct, vTw);
            Resource = arrayfun(@struct, vResource);
            PData = arrayfun(@struct, vPData);
            
            % 
            Resource.DisplayWindow = arrayfun(@struct, vDisplayWindow);
            Resource.RcvBuffer = arrayfun(@struct, vRcvBuffer);
            Resource.InterBuffer = arrayfun(@struct, vInterBuffer);
            Resource.ImageBuffer = arrayfun(@struct, vImageBuffer);
            Resource.Parameters = arrayfun(@struct, vParameters);
            
            
            TGC = arrayfun(@struct, vTGC);
            UI = arrayfun(@struct, vUI);
            
            % remove TX.aperture property is unneeded / unspecified
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
                Event(i).tx = safeIsMember([vEvent(i).tx], vTx);
                Event(i).rcv = safeIsMember([vEvent(i).rcv], vRcv);
                Event(i).recon = safeIsMember([vEvent(i).recon], vRecon);
                Event(i).process = safeIsMember([vEvent(i).process], vProcess);
                Event(i).seqControl = safeIsMember([vEvent(i).seqControl], vSeqControl);
            end
            
            % assign indices for TX
            for i = 1 : numel(TX)
                TX(i).waveform = safeIsMember([vTx(i).waveform], vTw);
            end
            
            % assign indices for Receive
            for i = 1 : numel(Receive)
                Receive(i).bufnum = safeIsMember([vRcv(i).bufnum], vRcvBuffer);
                Receive(i).TGC = safeIsMember([vRcv(i).TGC], vTGC);
            end
            
            % assign indices for Recon
            for i = 1 : numel(Recon)
                Recon(i).RINums = safeIsMember([vRecon(i).RINums], vReconInfo); 
            end
            
            % checking UI for empty values
            for i = 1 : numel(UI)
                for f = string(fieldnames(UI))'
                    if isempty(UI(i).(f))
                        UI(i).(f) = 0;
                    end
                end
            end

            % ----seqcontrol
            for i = 1:numel(SeqControl)
                if SeqControl(i).command == "jump"
                   SeqControl(i).argument = safeIsMember([SeqControl(i).argument], vEvent);  
                end
                if SeqControl(i).argument == 0
                    SeqControl(i).argument = [];
                end
            end
            
            nms  = {'Event', 'TX', 'Receive', 'Recon', 'ReconInfo', 'Process', 'SeqControl', 'TW', 'Resource', 'PData', 'TGC', 'UI'};
            vals = {Event, TX, Receive, Recon, ReconInfo, Process, SeqControl, TW, Resource, PData, TGC, UI};
            args = [nms; vals];
            vStruct = struct(args{:});
            
            % remove empty structures
            for f = string(fieldnames(vStruct))'
                if isempty(vStruct.(f))
                    vStruct = rmfield(vStruct, f);
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

    if isstruct(x)
        recursiveString2Char(x) % do recursive string2char on all props of x
    elseif isstring(x)
        y = char(x);
    else
        % do nothing
    end

end

