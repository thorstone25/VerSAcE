classdef VSXParameters < matlab.mixin.Copyable
    properties
        % Connector (1,:) double = 1
        numTransmit (1,1) double = 128
        numRcvChannels (1,1) double = 128
        speedOfSound (1,1) double = 1540
% %         speedCorrectionFactor (1,1) double = 1.0
        % startEvent (1,1) double = 1
        verbose (1,1) double {mustBeLessThanOrEqual(verbose,3)} = 2
        initializeOnly (1,1) double = 0
        simulateMode (1,1) double = 0
        waitForProcessing (1,1) double = 0
        % UpdateFunction (1,:) char = 'VsUpdate'
        % ProbeConnectorLED double
        % ProbeThermistor
        % SystemLED
        % numLogDataRecs (1,1) double
        % sizeApod (1,1) double
        % GUI (1,:) char = 'vsx_gui' 
    end
    methods
        function obj = VSXParameters(kwargs)
            arguments, kwargs.?VSXParameters, end
            for f = string(fieldnames(kwargs))', obj.(f) = kwargs.(f); end
        end
    end
end

        