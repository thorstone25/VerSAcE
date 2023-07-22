classdef VSXTX < matlab.mixin.Copyable
    properties
        waveform VSXTW {mustBeScalarOrEmpty} = VSXTW.empty
        Origin (1,3) double
        aperture double {mustBeScalarOrEmpty} = []
        Apod (1,:) double 
        focus (1,1) double = nan
        Steer (1,2) double = nan
        FocalPt (1,3) double = nan
        % FocalPtMm (1,3) double
        Delay (1,:) double
        TXPD (:,:,:,:) uint16
        % %         peakCutOff (1,1) double
        % %         peakBLMax (1,1) double
        %         VDASApod (1,:) double
        % %         VDASDelay (1,:) double = 0
    end
    methods
        function obj = VSXTX(kwargs)
            arguments
                kwargs.?VSXTX
            end
            for f = string(fieldnames(kwargs))'
                obj.(f) = kwargs.(f); 
            end
%             if length(obj.TXPD)>length(VSXPData)
%                 error('transmit pixel data out of range.')
%             end
        end
    end
end
