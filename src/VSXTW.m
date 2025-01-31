classdef VSXTW < matlab.mixin.Copyable
    properties
        % all
        type (1,1) string {mustBeMember(type, [ ...
            "parametric",...
            "envelope",...
            "pulseCode",... deprecated - will be spat back out
            "states",...
            "function",...
            "sampled"])} = "parametric";
        States (:,2) double
        TriLvlWvfm (:,1) double {mustBeMember(TriLvlWvfm, [-1, 0, 1])}
        Wvfm1Wy (:,1) double
        Wvfm2Wy (:,1) double
        peak (1,1) double
        numsamples (1,1) double 
        estimatedAvgFreq (1,1) double
        Bdur (1,1) double
        sysExtendBL (1,1) double % true if more than 25 cycle pulse
        CumOnTime (1,1) double
        Numpulses (1,1) double {mustBeInteger}
        integralPkUsec (1,2) double
        fluxHVlimit double {mustBeScalarOrEmpty} 
        VDASArbwave uint16
        
        % added?
        perChWvfm (1,1) logical = false
        refchnum (1,1) double
        refRelPulseWidth (1,1) double
        maxPulseusec (1,1) double
        Zload (1,1) double
        Zsource (1,1) double
        chIpk1V (1,1) double
        wvfmIntegral (1,1) double
        simChNum double {mustBeInteger, mustBeScalarOrEmpty} = []
        %TODO: use subclasses

        % parametric
        Parameters (:,4) double = []
        equalize (1,1) {mustBeMember(equalize, [0 1 2])} = 1

        % envelope
        envNumCycles  (1,:) double {mustBeInteger, mustBeInRange(envNumCycles, 1, 10000)} = []
        envFrequency  (1,:) double {mustBeInRange(envFrequency, 0.5, 32)} = []
        envPulseWidth (1,:) double {mustBeInRange(envPulseWidth, -1, 1)} = []

        % function/sampled
        frequency (1,1) = NaN
        Descriptors (1,4) = NaN
        Waveform (1,:) = []        
    end
       
    methods
        function obj = VSXTW(kwargs)
            arguments, kwargs.?VSXTW, end
            for f = string(fieldnames(kwargs))'
                obj.(f) = kwargs.(f);
            end
        end
    end
end
