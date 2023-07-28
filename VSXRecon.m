classdef VSXRecon < matlab.mixin.Copyable
    properties
        senscutoff (1,1) double = 0.6
        pdatanum (1,1) VSXPData
        rcvBufFrame (1,1) double = -1
        IntBufDest (1,1) VSXInterBuffer
        IntBufDestFrm (1,1) double = 1
        ImgBufDest (1,1) VSXImageBuffer
        ImgBufDestFrm (1,1) double = -1
        RINums (1,:) VSXReconInfo = VSXReconInfo.empty
        newFrameTimeout (1,1) double = 1000 % timeout in ms
        rcvLUT (:,:) uint16 = []
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
