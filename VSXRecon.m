classdef VSXRecon < matlab.mixin.Copyable
    properties
        senscutoff (1,1) double
        pdatanum (1,1) double
        rcvBufFrame (1,1) double
        IntBufDest (1,2) double
        ImgBufDest (1,2) double 
        RINums (1,:) VSXReconInfo
    end
    methods
        function obj = VSXRecon(kwargs)
            arguments
                kwargs.?VSXRecon
            end
            for f = string(fieldnames(kwargs))'
                obj.(f) = kwargs.(f); 
            end
        end
    end
end
