classdef VSXTX < matlab.mixin.Copyable
    properties
        waveform VSXTW {mustBeScalarOrEmpty} = VSXTW.empty
        Origin (1,3) double = 0
        aperture double {mustBeScalarOrEmpty, mustBeInteger, mustBePositive} = []
        Apod (1,:) double  = []
        focus (1,1) double = 0
        Steer (1,2) double = 0
        FocalPt (1,3) double = nan
        % FocalPtMm (1,3) double
        Delay (1,:) double = []
        TXPD (:,:,:,:) uint16 = [] % size(PData) x 3
        peakCutOff double {mustBeScalarOrEmpty} = 1.0 % [] % defaults from showTXPD
        peakBLMax  double {mustBeScalarOrEmpty} = 4.0 % [] % defaults from showTXPD
        %         VDASApod (1,:) double
        %         VDASDelay (1,:) double = 0
    end
    methods
        function obj = VSXTX(kwargs)
            arguments
                kwargs.?VSXTX
            end
            for f = string(fieldnames(kwargs))'
                obj.(f) = kwargs.(f); 
            end
        end
    end
end
