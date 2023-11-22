classdef VSXTGC < matlab.mixin.Copyable
    properties
        CntrlPts (1,8) double
        rangeMax (1,1) double
        Waveform (1,512) double
    end
    methods
        function obj = VSXTGC(kwargs)
            arguments
                kwargs.?VSXTGC
            end
            for f = string(fieldnames(kwargs))'
                obj.(f) = kwargs.(f); 
            end
        end
        function wv = computeTGCWaveform(TGC, fc)
            arguments
                TGC (1,1) VSXTGC
                fc (1,1) double
            end
            Trans = struct('frequency', fc/1e6); %#ok<NASGU> used in computeTGCWaveform
            wv = computeTGCWaveform(struct(TGC));
        end
    end
end
